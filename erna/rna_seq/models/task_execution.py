from django.utils import timezone
from django.db import models
from .task import Task

class TaskExecutionManager(models.Manager):

    def initialize(self, project_id, task_id, command):
        '''
        1. new task is created and settings are saved.
        2. change parameters of commands
        In this case, the task is ready for execution
        '''
        task = Task.objects.get_task(project_id, task_id)
        return self.model.objects.create(task=task, \
                status='R', command=command)
    
    def run(self, id:int):
        '''
        run/execute a task
        '''
        return self.model.objects.filter(id=id).update(status='E', \
            start_time=timezone.now(), end_time=None)
    
    def pause(self, id):
        return self.model.objects.filter(id=id)\
            .update(status='P', end_time=timezone.now())

    def stop(self, id):
        return self.model.objects.filter(id=id)\
            .update(status='D', end_time=timezone.now())

    def fail(self, id):
        return self.model.objects.filter(id=id)\
            .update(status='F', end_time=timezone.now())

    def get_lastest_task(self, project_id, task_id):
        task = Task.objects.get_task(project_id, task_id)
        return self.model.objects.filter(task=task).last()
    
    def get_project_status(self, project_id):
        '''
        Returns the status of all lasted tasks given a project
        '''
        status = {}
        tasks = Task.objects.get_project_tasks(project_id)
        for task in tasks:
            last = self.model.objects.filter(task=task).last()
            status[last.id] = last.status
        return status

    def get_command(self, id):
        return self.model.objects.get(id=id).command


class TaskExecution(models.Model):
    '''
    one task may be executed 0-many times.
    only store lastest execution given a task
    '''
    task = models.ForeignKey(
        Task,
        related_name = 'task_executions',
        on_delete=models.CASCADE
    )
    status = models.CharField(
        max_length=10,
        default = 'suspend',
        choices=[
            #newly added, parameters is added
            ('suspend', 'suspend'), 
            # relying tasks are finished. being ready for run
            ('ready', 'ready'),
            #pause this task if task is not launched.
            ('pause', 'pause'), 
            ('run', 'run'),
            ('finish', 'finish'),
            ('fail', 'fail'),
        ],
    )
    # string type converted from json format
    command = models.CharField(
        max_length=5000,
        verbose_name="command for tool launching"
    )
    start_time = models.DateTimeField(blank=True, null=True)
    end_time = models.DateTimeField(blank=True, null=True)

    objects = TaskExecutionManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('task', 'id', 'start_time')
