'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
import json
from rna_seq.models import *
from .align import Align


class ExecuteTask:
    def __init__(self, project_id:str, task_id:str):
        self.project_id = project_id
        self.task_id = task_id
        self.params = {}

    def __call__(self):
        self.retreive_metadata()
        # for k in self.params:
        #     print(f"{k}\t\t{self.params[k]}")
        return self.run_task()

    def retreive_metadata(self):
        '''
        retrieve all data from DB related project+task
        '''
        # Project
        project = Project.objects.get(project_id=self.project_id)
        self.params['project'] = vars(project)
        # Genome
        self.params['genome'] = Genome.objects.get(pk=project.genome.pk)
        # Task
        task = Task.objects.get(project=project, task_id=self.task_id)
        self.params['task'] = vars(task)
        self.params['task']['params'] = json.loads(task.params)

        # Method and Tool
        method_tool = MethodTool.objects.get(pk=task.method_tool.pk)
        self.params['method'] = Method.objects.get(pk=method_tool.method.pk)
        self.params['tool'] = Tool.objects.get(pk=method_tool.tool.pk)
        print(self.params['tool'].exe_name)

        # Samples
        project_samples = SampleProject.objects.filter(project=project)
        sample_files = [ obj.sample_file for obj in project_samples]
        self.params['sample_files'] = sample_files
    


    def run_task(self):
        '''
        '''
        match self.params['method'].method_name:
            case 'build_index':
                return Align(self.params).build_index()
            case 'align_transcriptome':
                return Align(self.params).align_transcriptome()

        return None