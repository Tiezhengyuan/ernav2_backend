'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
import json
import os
from django.conf import settings
from rna_seq.models import *
from .align import Align
from utils.dir import Dir

RESULTS_DIR = settings.RESULTS_DIR


class ExecuteTask:
    def __init__(self, project_id:str, task_id:str):
        self.project_id = project_id
        self.task_id = task_id
        self.params = {}

    def __call__(self):
        self.retreive_metadata()
        self.init_task()
        self.run_task()
        self.end_task()
        return None

    def retreive_metadata(self):
        '''
        retrieve all data from DB related project+task
        '''
        # Project
        project = Project.objects.get(project_id=self.project_id)
        self.params['project'] = project

        # Task
        task = Task.objects.get(project=project, task_id=self.task_id)
        self.params['task'] = task
        self.params['task'].params = json.loads(task.params)
        # parent/children
        parents = TaskTree.objects.filter(child=task)
        self.params['parents'] = [t.task for t in parents]
        children = TaskTree.objects.filter(task=task)
        self.params['children'] = [t.child for t in children]

        # Method and Tool
        method_tool = MethodTool.objects.get(pk=task.method_tool.pk)
        self.params['method'] = Method.objects.get(pk=method_tool.method.pk)
        self.params['tool'] = Tool.objects.get(pk=method_tool.tool.pk)

        # Samples
        project_samples = SampleProject.objects.filter(project=project)
        sample_files = [ obj.sample_file for obj in project_samples]
        self.params['sample_files'] = sample_files

        # Genome
        self.params['genome'] = Genome.objects.get(pk=project.genome.pk)

        # for k in self.params:
        #     print(f"{k}\t\t{self.params[k]}")


    def init_task(self):
        '''
        '''
        # output dir
        self.params['output_dir'] = os.path.join(
            RESULTS_DIR,
            self.params['project'].project_id,
            self.params['method'].method_name
        )
        Dir(self.params['output_dir']).init_dir()
        
        #update TaskExecution
        self.params['task_execution'] = TaskExecution.objects\
            .run(self.params['task'])

    def end_task(self):
        self.params['task_execution'] = TaskExecution.objects\
            .finish(self.params['task_execution'].id)
        
    def run_task(self):
        '''
        '''
        match self.params['method'].method_name:
            case 'build_index':
                return Align(self.params).build_index()
            case 'align_transcriptome':
                return Align(self.params).align_transcriptome()
        return None