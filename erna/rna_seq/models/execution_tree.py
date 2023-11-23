from django.db import models

from .task_execution import TaskExecution

class ExecutionTreeManager(models.Manager):
    pass

class ExecutionTree(models.Model):
    execution = models.ForeignKey(
        TaskExecution,
        on_delete=models.CASCADE,
        verbose_name='executed task',
    )
    child_execution = models.ForeignKey(
        TaskExecution,
        on_delete=models.CASCADE,
        related_name='child_executions',
        verbose_name='Child of executed task',
    )

    objects = ExecutionTreeManager()

    class Meta:
        app_label = 'rna_seq'
