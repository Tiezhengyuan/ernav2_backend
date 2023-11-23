'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
import os
from django.conf import settings
from django.utils import timezone

from rna_seq.models import Project, Task, TaskTree, TaskExecution,\
    ExecutionTree, MethodTool, Tool, Method, SampleProject, Genome
from .align import Align
from utils.dir import Dir


class ExecuteTasks:
    def __init__(self, project_id:str, task_id:str=None, chain:bool=None):
        self.pool = [(project_id, task_id, None)]
        self.chain = chain if chain else False


    def __call__(self) -> bool:
        while self.pool:
            project_id, task_id, parent_params = self.pool.pop(0)
            params = self.retreive_metadata(project_id, task_id)
            # stop execution if the task is occupied.
            if self.skip_task(params):
                return False
            # try to run the task
            params['parent_params'] = parent_params
            self.init_task(params)
            self.run_task(params)
            self.end_task(params)

            # add the next task into self.pool
            if self.chain and params.get('children'):
                for child in params['children']:
                    self.pool.append(
                        (child.project.project_id, child.task_id, params)
                    ) 
        return True

    def skip_task(self, params:dict) -> bool:
        '''
        confirm if the task is locked due to execution
        '''
        if params['task'].task_execution and \
            params['task'].task_execution.status == 'run' and \
            params['task'].task_execution.command:
            print(f'''
                Skipp the task  because that is running.
                task id={params['task'].id}
                execution id={params['task'].task_execution.id},
                CMD={params['task'].task_execution.command}.
                ''')
            return True
        return False

    def retreive_metadata(self, project_id:str, task_id:str):
        '''
        retrieve all data from DB related project+task
        '''
        params = {}
        # Project
        project = Project.objects.get(project_id=project_id)
        params['project'] = project

        # Task
        task = Task.objects.get(project=project, task_id=task_id) \
            if task_id else Task.objects.filter(project=project).first()
        params['task'] = task

        # parent/children
        parents = TaskTree.objects.filter(child=task)
        params['parents'] = [t.task for t in parents]
        children = TaskTree.objects.filter(task=task)
        params['children'] = [t.child for t in children]

        # Method and Tool
        method_tool = MethodTool.objects.get(pk=task.method_tool.pk)
        params['method'] = Method.objects.get(pk=method_tool.method.pk)
        params['tool'] = Tool.objects.get(pk=method_tool.tool.pk)

        # Samples
        project_samples = SampleProject.objects.filter(project=project)
        sample_files = [ obj.sample_file for obj in project_samples]
        params['sample_files'] = sample_files

        # Genome
        params['genome'] = Genome.objects.get(pk=project.genome.pk)

        # for k in params:
        #     print(f"{k}\t\t{params[k]}")
        return params


    def init_task(self, params:dict):
        '''
        '''
        # variables
        params['cmd'] = None
        params['output'] = []
        params['output_dir'] = os.path.join(
            settings.RESULTS_DIR,
            params['project'].project_id,
            params['method'].method_name
        )
        Dir(params['output_dir']).init_dir()
        
        #update TaskExecution
        params['task_execution'] = TaskExecution.objects.run(params['task'])
        # update Task
        params['task'].task_execution = params['task_execution']
        params['task'].save()
        #update ExecutionTree
        if params.get('parent_params'):
            ExecutionTree.objects.create(
                execution = params['parent_params']['task_execution'],
                child_execution = params['task_execution'],
            )
        
    def run_task(self, params:dict):
        match params['method'].method_name:
            case 'build_index':
                return Align(params).build_index()
            case 'align_transcriptome':
                return Align(params).align_transcriptome()
        return None
    
    def end_task(self, params:dict):
        print(params['task_execution'].id, params['output'])
        params['task_execution'].end_execution(
            params.get('cmd'), params.get('output')
        )

    
