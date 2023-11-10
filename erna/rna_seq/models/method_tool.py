'''

'''
from django.db import models
from .tool import Tool
from .constants import METHODS
from .method import Method

class MethodToolManager(models.Manager):
  def refresh(self):
    self.all().delete()
    res = []
    for method in METHODS:
      method_obj = Method.objects.get(method_name=method['method_name'])
      if method['tool_name']:
        for tool_name in method['tool_name']:
          tools = Tool.objects.filter(tool_name=tool_name)
          for tool in tools:
            obj = self.create(method=method_obj, tool = tool)
            res.append(obj)
      else:
        obj = self.create(method=method_obj)
        res.append(obj)
    return res

  def get_method_tool(self, method_name:str, tool_name:str, version:str=None):
    method = Method.objects.get(method_name=method_name)
    tool = Tool.objects.get_tool(tool_name=tool_name, version=version)
    return self.get(method=method, tool=tool)

class MethodTool(models.Model):
  method = models.ForeignKey(
    Method,
    related_name = 'tools',
    on_delete=models.CASCADE
  )
  tool = models.ForeignKey(
    Tool,
    related_name = 'methods',
    on_delete=models.CASCADE,
    null=True,
    blank=True
  )

  objects = MethodToolManager()

  class Meta:
    app_label = 'rna_seq'
    unique_together = ['method', 'tool']
    ordering = ['method', 'tool']