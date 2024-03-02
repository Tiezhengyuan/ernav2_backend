'''

'''
import json
from django.db import models

from .tool import Tool
from .method import Method
from rna_seq.constants import METHODS

# model manager
class MethodToolManager(models.Manager):
  def refresh(self):
    self.all().delete()
    res = []
    for method in METHODS:
      method_obj = Method.objects.get(method_name=method['method_name'])
      if method.get('exe_name'):
        for exe_name in method['exe_name']:
          tools = Tool.objects.filter(exe_name=exe_name)
          for tool in tools:
            obj = self.create(method=method_obj, tool = tool)
            res.append(obj)
      else:
        obj = self.create(method=method_obj)
        res.append(obj)
    return res

  def get_method_tool(self, method_name:str, exe_name:str=None, version:str=None):
    method = Method.objects.get(method_name=method_name)
    if exe_name or version:
      tool = Tool.objects.get_tool(exe_name=exe_name, version=version)
      if tool:
        return self.get(method=method, tool=tool)
    return self.get(method=method)

  def get_method_tools(self):
    '''
    export method_tools
    '''
    res = {}
    for obj in self.all():
      method_name = obj.method.method_name
      if method_name not in res:
          res[method_name] = []
      item = {}
      value = {
          'method_tool_id': obj.id,
          'method_name': obj.method.method_name,
      }
      if obj.tool:
          default_params = json.loads(obj.tool.default_params) if \
              obj.tool.default_params else {}
          value.update({
              'tool_name': obj.tool.tool_name,
              'exe_name': obj.tool.exe_name,
              'version': obj.tool.version,
              'params': default_params,
          })
          item = {
              'text': f"{obj.tool.tool_name}_{obj.tool.version}",
              'value': value,
          }
      else:
          if obj.method.default_params:
              value['params'] = json.loads(obj.method.default_params) \
                  if obj.method.default_params else {}
          item = {
              'text': 'built-in',
              'value': value,
          }
      # include empty methods of which no tools are defined.
      res[method_name].append(item)
    return res

# model
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
