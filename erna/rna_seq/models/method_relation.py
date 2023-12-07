'''

'''
from django.db import models
from .method import Method
from rna_seq.constants import METHODS


class MethodRelationManager(models.Manager):
  def refresh(self):
    '''
    load constant METHODS into db
    '''
    self.all().delete()
    res = []
    for method in METHODS:
      name = method['method_name']
      for child in method.get('child_method', []):
        obj = self.create(
          method=Method.objects.get(method_name=name),
          child=Method.objects.get(method_name=child)
        )
        res.append(obj)
    return res


class MethodRelation(models.Model):
  method = models.ForeignKey(
    Method,
    on_delete=models.CASCADE
  )
  child = models.ForeignKey(
    Method,
    related_name = 'children',
    on_delete=models.CASCADE
  )

  objects = MethodRelationManager()

  class Meta:
    app_label = 'rna_seq'
    ordering = ['method', 'child']