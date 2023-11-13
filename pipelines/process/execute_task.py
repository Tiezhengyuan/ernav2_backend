'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''
import json
from rna_seq.models import Project, Task, MethodTool, \
    Method, Tool, Genome
from .align import Align


class ExecuteTask:
    def __init__(self, project_id:str, task_id:str):
        self.project_id = project_id
        self.task_id = task_id
        self.params = {}

    def __call__(self):
        self.retreive_metadata()
        print(self.params)
        # self.run_task()

    def retreive_metadata(self):
        '''
        retrieve all data from DB related project+task
        '''
        project = Project.objects.get(project_id=self.project_id)
        self.params['project'] = vars(project)
        genome = Genome.objects.get(pk=project.genome.pk)
        self.params['genome'] = vars(genome) if genome else {}
        task = Task.objects.get(project=project, task_id=self.task_id)
        self.params['task'] = vars(task)
        self.params['task']['params'] = json.loads(task.params)
        method_tool = MethodTool.objects.get(pk=task.method_tool.pk)
        method = Method.objects.get(pk=method_tool.method.pk)
        self.params['method'] = vars(method)
        tool = Tool.objects.get(pk=method_tool.tool.pk)
        self.params['tool'] = vars(tool)



    def run_task(self):
        '''
        '''
        match self.params['method']['method_name']:
            case 'build_index':
                return Align(self.params).build_index()
            case 'align_transcriptome':
                return Align(self.params).align_transcriptome()

        return None