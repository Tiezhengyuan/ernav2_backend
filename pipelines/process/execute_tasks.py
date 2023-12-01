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
        self.pool = [(project_id, task_id, None),]
        self.chain = True if (chain or task_id is None) else False
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

    def run_task(self, params:dict)->None:
        match params['method'].method_name:
            case 'import_data':
                # task T00
                from .collect import Collect
                return Collect(params).import_data()
            case 'trim_sequences':
                from .trim_adapter import TrimAdapter
                return TrimAdapter(params)()
            case 'build_index':
                return Align(params).build_index()
            case 'build_genome_index':
                return Align(params).build_genome_index()
            case 'align_transcriptome':
                return Align(params).align()
            case 'align_short_reads':
                return Align(params).align()
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
            case 'convert_format':
                from .convert_format import ConvertFormat
                return ConvertFormat(params).sam_to_bam()

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
        task = Task.objects.get(project=project, task_id=task_id) if \
            task_id else Task.objects.get_head_task(project=project)
        params['task'] = task

        # parent/children
        parents = TaskTree.objects.filter(child=task)
        parents = [t.task for t in parents]
        params['parents'] = parents
        params['parent_outputs'] = self.combine_parents_output(parents)
        print(params['parents'])
        print('###', params['parent_outputs'])
        children = TaskTree.objects.filter(task=task)
        params['children'] = [t.child for t in children]

        # Method and Tool
        if task.method_tool:
            method_tool = MethodTool.objects.get(pk=task.method_tool.pk)
            params['method'] = Method.objects.get(pk=method_tool.method.pk)
            params['tool'] = Tool.objects.get(pk=method_tool.tool.pk) \
                if method_tool.tool else None

        # Genome
        if project.genome:
            params['genome'] = Genome.objects.get(pk=project.genome.pk)
            params['annotations'] = Annotation.objects.filter(genome=params['genome'])
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
            f"{params['task'].task_id}_{params['method'].method_name}",
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

    
