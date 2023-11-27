from django.db import models

from .method_tool import MethodTool

class PipelineManager(models.Manager):

    def load_pipeline(self, pipeline_name:str, steps:list):
        res = []
        child = None
        while steps:
            method_name, exe_name, tool_version = steps.pop()
            method_tool = MethodTool.objects.get_method_tool(
                method_name, exe_name, tool_version
            )
            data = {
                'pipeline_name': pipeline_name,
                'step': method_tool
            }
            if child is not None:
                data['next_step'] = child
            child = method_tool
            obj = self.create(**data)
            res.append(obj)
        return res


class Pipeline(models.Model):
    pipeline_name = models.CharField(max_length=128)
    step = models.ForeignKey(
        MethodTool,
        related_name = 'steps',
        on_delete=models.CASCADE
    )
    next_step = models.ForeignKey(
        MethodTool,
        related_name = 'next_steps',
        on_delete=models.CASCADE,
        null=True,
        blank=True
    )

    objects = PipelineManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ['pipeline_name', 'step']