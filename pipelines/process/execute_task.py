'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
import json
import os
from django.conf import settings
from django.utils import timezone

from rna_seq.models import Project, Task, TaskTree, TaskExecution,\
    MethodTool, Tool, Method, SampleProject, Genome
from .align import Align
from utils.dir import Dir



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
            settings.RESULTS_DIR,
            self.params['project'].project_id,
            self.params['method'].method_name
        )
        Dir(self.params['output_dir']).init_dir()
        # 
        self.params['cmd'] = None
        self.params['output'] = []
        
        #update TaskExecution
        self.params['task_execution'] = TaskExecution.objects\
            .run(self.params['task'])

    def end_task(self):
        print(self.params['task_execution'].id, self.params['output'])
        # update TaskExecution
        if self.params.get('cmd'):
            TaskExecution.objects.update_command(
                self.params['task_execution'],
                self.params['cmd']
            )
            self.params['task_execution'].status = 'finish'
        else:
            self.params['task_execution'].status = 'skip'
        if self.params.get('output'):
            self.params['task_execution'].output = json.dumps(
                self.params['output'])
        self.params['task_execution'].end_time = timezone.now()
        self.params['task_execution'].save()
    
        
    def run_task(self):
        '''
        '''
        match self.params['method'].method_name:
            case 'build_index':
                return Align(self.params).build_index()
            case 'align_transcriptome':
                return Align(self.params).align_transcriptome()
        return None