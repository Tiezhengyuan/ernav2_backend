import json
from django.db import models

from .project import Project
from .method import Method
from .method_tool import MethodTool
from .genome import Genome
from .annotation import Annotation

class TaskManager(models.Manager):

    def next_task_id(self, project_id):
        '''
        first confirm project_id, then consecutive on task_id
        maximum 99 tasks per project
        '''
        project = Project.objects.get(project_id=project_id)
        last = self.filter(project=project).last()
        if last:
            next_id = str(int(last.task_id[1:]) + 1)
            return f"T{next_id.zfill(2)}"
        return 'T01'
    
    def get_head_task(self, project):
        return self.get(project=project, task_id='T00')

    def add_head_task(self, project):
        '''
        ID of head task is always T00
        that is automatically created
        '''
        print(Method.objects.head_method())
        method_tool = MethodTool.objects.get(
            method=Method.objects.head_method()
        )
        defaults = {
            'method_tool': method_tool,
            'task_name': 'head_task',
            'params': json.dumps({}),
            'is_ready': False,
        }
        task = self.update_or_create(
            project=project,
            task_id='T00',
            defaults=defaults
        )
        return task


    def add_task(self, project, task_id, data:dict):
        '''
        task_id could be generated automatically
        The task could be newly created or updated
        '''
        method_tool = None
        if 'method_tool_id' in data:
            method_tool = MethodTool.objects.get(pk=data['method_tool_id'])
        elif 'method_name' in data:
            if 'tool' in data:
                method_tool = MethodTool.objects.get_method_tool(
                    data['method_name'],
                    data['tool']['exe_name'],
                    data['tool']['version']
                )
            else:
                method_tool = MethodTool.objects.get_method_tool(
                    data['method_name']
                )

        annot = None
        if 'genome' in data:
            genome = Genome.objects.get(**data['genome'])
            annot = Annotation.objects.get(
                genome=genome,
                file_format=data['annotation']['file_format'],
                annot_type=data['annotation']['annot_type']
            )

        defaults = {
            'task_name': data.get('task_name'),
            'params': json.dumps(data.get('params', {})),
            'is_ready': True if data.get('is_ready') else False,
            'method_tool': method_tool,
            'annotation': annot,
        }
        print(defaults)
        task = self.update_or_create(project=project,
            task_id=task_id, defaults=defaults)
        return task

    def load_tasks(self, project_id:str, tasks_data:list):
        '''
        load multiple tasks given a existing project
        '''
        # get Project
        project = Project.objects.get(project_id=project_id)

        # first add head task
        self.add_head_task(project=project)

        # add tasks
        for data in tasks_data:
            task_id = data.get('task_id')
            task = self.add_task(project, task_id, data)
        else:
            # update Project
            project.status = 'ready'
            project.save()
        return self.filter(project=project_id)

    def delete_tasks(self, project_id:str=None, task_id:str=None):
        '''
        delete all/some tasks, or single task
        '''
        if project_id:
            # single task
            if task_id:
                tasks = self.filter(project_id=project_id, task_id=task_id)
                return tasks.delete()
            # project tasks
            else:
                tasks = self.filter(project_id=project_id)
                return tasks.delete()
        else:
            # all tasks
            if task_id is None:
                return self.all().delete()
        return []


class Task(models.Model):
    # project_id + task_id = pk
    project = models.ForeignKey(
        'rna_seq.Project',
        related_name = 'tasks',
        on_delete=models.CASCADE,
        verbose_name="Project ID"
    )
    task_id = models.CharField(
        max_length=10,
        verbose_name="Task ID"
    )
    # some conditions should be met
    # at least "method_name" should be defined
    is_ready = models.BooleanField(default=False)
    method_tool = models.ForeignKey(
        'rna_seq.MethodTool',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        verbose_name = "Method and Tool",
    )
    annotation = models.ForeignKey(
        'rna_seq.Annotation',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    task_name = models.CharField(
        max_length=100,
        null=True,
        blank=True, 
        verbose_name="Task Name",
    )
    params = models.CharField(
        max_length=1256,
        null=True,
        blank=True, 
        verbose_name="Parameters (in json string format)",
    )
    task_execution = models.OneToOneField(
        'rna_seq.TaskExecution',
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='task_execution'
    )
    create_time = models.DateTimeField(auto_now_add=True)
    modified_time = models.DateTimeField(auto_now=True)

    objects = TaskManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('project', 'task_id',)
        unique_together = ('project', 'task_id')

    def update_params(self, input:dict) -> dict:
        params = json.loads(self.params) if self.params else {}
        if input:
            params.update(input)
            self.params = json.dumps(params)
            self.save()
        return self.params
    
    def get_params(self) -> dict:
        try:
            return json.loads(self.params)
        except Exception as e:
            print(e)
        return self.params

