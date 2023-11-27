'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
from copy import deepcopy
import os
from django.conf import settings

from rna_seq.models import Project, Task, TaskTree, TaskExecution,\
    ExecutionTree, MethodTool, Tool, Method, SampleProject, Genome, \
    Annotation
from rna_seq.models.constants import METHODS
from utils.dir import Dir
from .align import Align
from .assemble import Assemble



class ExecuteTasks:
    def __init__(self, project_id:str, task_id:str=None, chain:bool=None, force:bool=None):
        self.pool = []
        self.chain = True
        if task_id:
            self.pool = [(project_id, task_id, None)]
            self.chain = eval(chain) if chain else False
        else:
            for head_task_id in self.detect_head_tasks(project_id):
                self.pool.append((project_id, head_task_id, None))
        self.force = force if force else True

    def detect_head_tasks(self, project_id:str):
        head_task_ids = []
        tasks = Task.objects.filter(project_id=project_id)
        for task in tasks:
            obj = TaskTree.objects.filter(child=task)
            if not obj:
                head_task_ids.append(task.task_id)
        return head_task_ids

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

    def run_task(self, params:dict)->None:
        match params['method'].method_name:
            case 'build_index':
                return Align(params).build_index()
            case 'build_genome_index':
                return Align(params).build_genome_index()
            case 'align_transcriptome':
                return Align(params).align_transcriptome()
            case 'convert_format':
                from .convert_format import ConvertFormat
                return ConvertFormat(params).sam_to_bam()
            case 'assemble_transcripts':
                return Assemble(params).assemble_transcripts()
            case 'merge_transcripts':
                return Assemble(params).merge_transcripts()
            case 'count_reads':
                from .collect import Collect
                return Collect(params).count_reads()
            case 'quality_control':
                from .quality_control import QualityControl
                return QualityControl(params)()

        return None


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
        parents = [t.task for t in parents]
        params['parents'] = parents
        params['parent_outputs'] = self.combine_parents_output(parents)
        children = TaskTree.objects.filter(task=task)
        params['children'] = [t.child for t in children]

        # Method and Tool
        if task.method_tool:
            method_tool = MethodTool.objects.get(pk=task.method_tool.pk)
            params['method'] = Method.objects.get(pk=method_tool.method.pk)
            params['tool'] = Tool.objects.get(pk=method_tool.tool.pk) \
                if method_tool.tool else None

        # Samples
        project_samples = SampleProject.objects.filter(project=project)
        sample_files = [ obj.sample_file for obj in project_samples]
        params['sample_files'] = sample_files

        # Genome
        if project.genome:
            params['genome'] = Genome.objects.get(pk=project.genome.pk)
            params['annotations'] = Annotation.objects.filter(genome=params['genome'])
        # for k in params:
        #     print(f"{k}\t\t{params[k]}")
        return params

    def combine_parents_output(self, parents:list) -> list:
        '''
        suppose that a task has one or multiple parents.
        '''
        if not parents:
            return []
        parent_output = deepcopy(parents[0].task_execution.get_output())
        if len(parents) > 1:
            for parent in parents[1:]:
                another = deepcopy(parent.task_execution.get_output())
                for a,b in zip(parent_output, another):
                    a.update(b)
        return parent_output

    def init_task(self, params:dict) -> None:
        '''
        '''
        # variables
        params['cmd'] = None
        params['force_run'] = True
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
        print('#end task:', params['task_execution'].id)
        params['task_execution'].end_execution(params.get('output'))
        return None

    
