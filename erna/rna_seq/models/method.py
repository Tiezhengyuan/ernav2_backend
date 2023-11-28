'''

'''
from django.db import models
from .constants import METHODS


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
      }
      obj = self.create(**data)
      res.append(obj)
    return res

  def head_method(self):
    return self.get(method_name='import_data')

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

  objects = MethodManager()

  class Meta:
    app_label = 'rna_seq'
    ordering = ['method_name',]