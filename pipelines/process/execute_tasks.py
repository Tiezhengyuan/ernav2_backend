'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
from copy import deepcopy
import json
import os
from django.conf import settings
from rnaseqdata import load_seqdata


from rna_seq.models import Project, Task, TaskTree, TaskExecution,\
    ExecutionTree, MethodTool, Tool, Method, Genome, Annotation
from rna_seq.constants import METHODS
from pipelines.utils.dir import Dir
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
            params = self.init_params(project_id, task_id)
            # postpone the task if parent tasks are not ready
            if self.postpone_task(params, i):
                continue
            
            # update params if the task seems being ready.
            params = self.retreive_metadata(params)
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
                    tag = 1
                    for item in self.pool:
                        if item[1] == child.task_id:
                            tag = 0
                    if tag == 1:
                        next_task = (project_id, child.task_id, params)
                        self.pool.append(next_task)
            i += 1
        return True

    def init_params(self, project_id:str, task_id:str) -> dict:
        '''
        get project and task instance 
        '''
        project = Project.objects.get(project_id=project_id)
        task = Task.objects.get(project=project, task_id=task_id)
        parents = TaskTree.objects.filter(child=task)
        # load seqdata
        seqdata_file = os.path.join(settings.RESULTS_DIR,
            project.project_id, "T00_import_data", "seqdata.obj")
        # load seqdata
        seqdata_file = os.path.join(settings.RESULTS_DIR,
            project.project_id, "T00_import_data", "seqdata.obj")
        params = {
            'project': project,
            'task': task,
            'parents': [t.task for t in parents],
            'seqdata': load_seqdata(seqdata_file),
            'seqdata_path': seqdata_file,
        }
        return params

    
    def postpone_task(self, params:dict, i:int) -> bool:
        # Task T00 should be always executed
        if params['task'].task_id == 'T00':
            return False
        # the last task
        if i == len(self.pool) - 1:
            return False

        # parents exist
        if params['parents']:
            # get task_ids of all finished tasks
            finished_task_ids = [item[1] for item in self.pool[:i+1]]
            for parent in params['parents']:
                # at least one parent is not finished, then postpone the task
                if parent.task_id not in finished_task_ids:
                    if i < len(self.pool) -1:
                        self.pool[i:] = self.pool[i+1:] + [self.pool[i],]
                        return True
        return False

    def retreive_metadata(self, params:dict) -> dict:
        '''
        retrieve all data from DB related project+task
        '''
        # parent/children
        params['parent_outputs'] = self.combine_parents_output(params['parents'])
        children = TaskTree.objects.filter(task=params['task'])
        params['children'] = [t.child for t in children]

        # Method and Tool
        if params['task'].method_tool:
            method_tool = MethodTool.objects.get(pk=params['task'].method_tool.pk)
            params['method'] = Method.objects.get(pk=method_tool.method.pk)
            params['tool'] = Tool.objects.get(pk=method_tool.tool.pk) \
                if method_tool.tool else None

        # Genome
        if params['project'].genome:
            genome = Genome.objects.get(pk=params['project'].genome.pk)
            params['genome'] = genome
            params['annot_genomic_dna'] = Annotation.objects.genome_annot(genome, 'fna')
            #Note: GFF works, but GTF doesn't work in StringTie (bug)
            params['annot_genomic_gtf'] = Annotation.objects.genome_annot(genome, 'gff')
        return params

    def combine_parents_output(self, parents:list) -> list:
        '''
        suppose that a task has one or multiple parents.
        '''
        if not parents or not parents[0].task_execution:
            return []
        parent_output = deepcopy(parents[0].task_execution.get_output())
        if len(parents) > 1:
            for parent in parents[1:]:
                another = deepcopy(parent.task_execution.get_output())
                for a,b in zip(parent_output, another):
                    a.update(b)
        return parent_output

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

    def init_task(self, params:dict) -> None:
        '''
        '''
        print(f"Try to run task={params['task'].task_id}.")
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
            case 'merge_stringtie_read_counts':
                return Count(params).merge_stringtie_read_counts()

            case 'quality_control':
                from .quality_control import QualityControl
                return QualityControl(params)()

            case 'convert_format':
                from .convert_format import ConvertFormat
                return ConvertFormat(params).sam_to_bam()
        return None
    
    def end_task(self, params:dict) -> None:
        task_id = params['task'].task_id
        exec_id = params['task_execution'].id
        print(f"\nThe task {task_id} ends. execution={exec_id}\n\n")
        # update db.TaskExecution
        params['task_execution'].end_execution(params.get('output'))
        # save meta data into json
        outfile = os.path.join(params['output_dir'], 'output.json')
        with open(outfile, 'w') as f:
            json.dump(params['output'], f, indent=4)
        return None
  
