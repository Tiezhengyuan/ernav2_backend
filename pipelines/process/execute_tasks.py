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
from .collect import Collect
from .count import Count

class ExecuteTasks:
    def __init__(self, project_id:str, task_id:str=None, chain:bool=None, force:bool=None):
        self.pool = []
        self.chain = True if (chain or task_id is None) else False
        if not task_id:
            task_id = 'T00'
        self.pool = [(project_id, task_id, None),]
        self.force = force if force else True

    def __call__(self) -> bool:
        i = 0
        while i < len(self.pool):
            project_id, task_id, parent_params = self.pool[i]
            params = self.retreive_metadata(project_id, task_id)
            params['parent_params'] = parent_params
            
            # stop execution if the task is occupied.
            if self.skip_task(params):
                return False
            
            # postpone the task if parent tasks are not ready
            if 0 < i < len(self.pool) - 1:
                if self.postpone_task(params, i):
                    continue

            # try to run the task
            self.init_task(params)
            self.run_task(params)
            self.end_task(params)

            # add the next task into self.pool
            if self.chain and params.get('children'):
                for child in params['children']:
                    tag = 1
                    for item in self.pool:
                        if item[1] == child.task_id:
                            tag = 0
                    if tag == 1:
                        next_task = (project_id, child.task_id, params)
                        self.pool.append(next_task)
            i += 1
        return True

    def postpone_task(self, params:dict, i:int) -> bool:
        parent_tasks = [item[1] for item in self.pool[:i+1]]
        # if no parent is detected, just run the task
        if not parent_tasks:
            return False
        # parents exist
        for parent in params['parents']:
            # at least one parent is not finished, then postpone the task
            if parent.task_id not in parent_tasks:
                if i < len(self.pool) -1:
                    self.pool[i:] = self.pool[i+1:] + [self.pool[i],]
                return True
        return False

    def run_task(self, params:dict)->None:
        match params['method'].method_name:
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

            case 'import_data':
                # task T00
                return Collect(params).import_data()
            case 'merge_transcripts':
                return Collect(params).merge_transcripts()
            
            case 'count_reads':
                return Count(params).count_reads()
            case 'merge_read_counts':
                return Count(params).merge_read_counts()

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
        task = Task.objects.get(project=project, task_id=task_id)
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

        # Genome
        if project.genome:
            genome = Genome.objects.get(pk=project.genome.pk)
            params['genome'] = genome
            params['annot_genomic_dna'] = Annotation.objects.genome_annot(genome, 'fna')
            params['annot_genomic_gtf'] = Annotation.objects.genome_annot(genome, 'gtf')
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
        print(f"Try to run task={params['task'].task_id}")
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
        task_id = params['task'].task_id
        exec_id = params['task_execution'].id
        print(f"\nThe task {task_id} ends. execution={exec_id}\n\n")
        params['task_execution'].end_execution(params.get('output'))
        return None

    
