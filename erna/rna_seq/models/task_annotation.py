from django.db import models


from .task import Task
from .annotation import Annotation


class TaskAnnotationManager(models.Manager):
    pass

class TaskAnnotation(models.Model):
    annotation = models.ForeignKey(
        Annotation,
        on_delete=models.CASCADE
    )
    task = models.ForeignKey(
        Task,
        on_delete=models.CASCADE
    )

    objects = TaskAnnotationManager()

    class Meta:
        app_label = 'rna_seq'
        unique_together = ['annotation', 'task']
        ordering = ['annotation', 'task']
