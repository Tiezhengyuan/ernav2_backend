import json
from django.utils import timezone
from django.db import models
from .task import Task

class TaskExecutionManager(models.Manager):

    def run(self, task):
        '''
        run/execute a task
        '''
        return self.model.objects.create(task=task,
            status='run', start_time=timezone.now())

    def finish(self, id):
        '''
        task is executed successfully
        '''
        return self.model.objects.filter(id=id)\
            .update(status='finish', end_time=timezone.now())

    def skip(self, id):
        return self.model.objects.filter(id=id)\
            .update(status='skip', end_time=timezone.now())

    def pause(self, id):
        return self.model.objects.filter(id=id)\
            .update(status='pause', end_time=timezone.now())


    def fail(self, id):
        return self.model.objects.filter(id=id)\
            .update(status='fail', end_time=timezone.now())

    def lastest_execution(self, project_id, task_id):
        task = Task.objects.get_task(project_id, task_id)
        return self.filter(task=task).last()
    
    # def get_project_status(self, project_id):
    #     '''
    #     Returns the status of all lasted tasks given a project
    #     '''
    #     status = {}
    #     tasks = Task.objects.get_project_tasks(project_id)
    #     for task in tasks:
    #         last = self.model.objects.filter(task=task).last()
    #         status[last.id] = last.status
    #     return status

    def update_command(self, task_execution, new_command:list):
        '''
        run/execute a task
        '''
        command = json.loads(task_execution.command) if \
            task_execution.command else []
        command.append(' '.join(new_command))
        task_execution.command = json.dumps(command)
        task_execution.save()
        return task_execution

        
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
            # 
            ('run', 'run'),
            # 
            ('finish', 'finish'),
            ('skip', 'skip'),
            ('fail', 'fail'),
        ],
    )
    start_time = models.DateTimeField(blank=True, null=True)
    end_time = models.DateTimeField(blank=True, null=True)
    # string type converted from json format
    command = models.CharField(
        max_length=1025,
        null=True,
        blank=True,
        verbose_name="command for tool launching"
    )
    output = models.CharField(
        max_length=5000,
        null=True,
        blank=True,
        verbose_name="Path of output files"
    )

    objects = TaskExecutionManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('task', 'id', 'start_time')
