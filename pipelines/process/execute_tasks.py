'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
import os
from django.conf import settings

from rna_seq.models import Project, Task, TaskTree, TaskExecution,\
    ExecutionTree, MethodTool, Tool, Method, SampleProject, Genome
from rna_seq.models.constants import METHODS
from utils.dir import Dir
from .align import Align
from .convert_format import ConvertFormat
from .assemble import Assemble



class ExecuteTasks:
    def __init__(self, project_id:str, task_id:str=None, chain:bool=None, force:bool=None):
        self.pool = [(project_id, task_id, None)]
        self.chain = chain if chain else False
        self.force = force if force else True


    def __call__(self) -> bool:
        while self.pool:
            project_id, task_id, parent_params = self.pool.pop(0)
            params = self.retreive_metadata(project_id, task_id)
            params['parent_params'] = parent_params
            
            # stop execution if the task is occupied.
            if self.skip_task(params):
                return False
            
            # try to run the task
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

    def run_task(self, params:dict):
        match params['method'].method_name:
            case 'build_index':
                return Align(params).build_index()
            case 'align_transcriptome':
                return Align(params).align_transcriptome()
            case 'convert_format':
                return ConvertFormat(params).sam_to_bam()
            case 'assemble_transcripts':
                return Assemble(params).align_transcriptome()


    def skip_task(self, params:dict) -> bool:
        '''
        confirm if the task is locked due to execution
        '''
        methods = [i['method_name'] for i in METHODS]
        if params['method'].method_name not in methods:
            print(f"wrong method name. {params['method'].method_name}.")
            return False
        
        status = params['task'].task_execution.status if \
            params['task'].task_execution else None
        if status == 'run' and not self.force:
            print(f'''
                Skipp the task  because that is running.
                task id={params['task'].id}
                execution id={params['task'].task_execution.id},
                command={params['task'].task_execution.command}.
            ''')
            return True
        return False

    def retreive_metadata(self, project_id:str, task_id:str) -> dict:
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


    def init_task(self, params:dict) -> None:
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
        return None
        

    
    def end_task(self, params:dict) -> None:
        print('#end task:', params['task_execution'].id, params['output'])
        params['task_execution'].end_execution(params.get('output'))
        return None

    
