'''
model: Method
'''
import json
from django.db import models

from rna_seq.constants import METHODS, ROOT_METHOD

# manager
class MethodManager(models.Manager):

  def refresh(self):
    '''
    load methods defined in METHODS
    '''
    self.all().delete()
    res = []
    for method in METHODS:
      data = {
        'method_name': method['method_name'],
        'description': method.get('description'),
        'default_params': json.dumps(method['default_params']) if \
          'default_params' in method else None,
      }
      obj = self.create(**data)
      res.append(obj)
    return res

  def head_method(self):
    return self.get(method_name=ROOT_METHOD['method_name'])

  def method_names(self):
    res = []
    for method in self.all():
      obj = {
        'method_name': method.method_name,
        'label': method.method_name.title().replace('_', ' '),
        'component': method.method_name.title().replace('_', ''),
      }
      res.append(obj)
    return res

# model
class Method(models.Model):
  method_name = models.CharField(
    max_length=96,
    primary_key=True
  )
  description = models.CharField(
    max_length=1028,
    null=True,
    blank=True
  )
  # optional
  default_params = models.CharField(
      max_length=1028,
      blank=True,
      null=True,
      verbose_name='Default parameters'
  )
  
  objects = MethodManager()

  class Meta:
    app_label = 'rna_seq'
    ordering = ['method_name',]
  
