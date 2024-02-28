'''

'''
from django.db import models
from .method import Method
from rna_seq.constants import METHODS, ROOT_METHOD


class MethodRelationManager(models.Manager):
  def refresh(self):
    '''
    load constant METHODS into db
    '''
    self.all().delete()
    root = Method.objects.get(method_name=ROOT_METHOD['method_name'])

    res = []
    for method in METHODS:
      name = method['method_name']
      parent = Method.objects.get(method_name=name)
      # add root-parent
      obj = self.create(method=root, child=parent)
      res.append(obj)

      # add parent-child
      for child in method.get('child_method', []):
        obj = self.create(
          method=parent,
          child=Method.objects.get(method_name=child)
        )
        res.append(obj)
    return res

  def get_children(self):
      methods = Method.objects.all()
      res = dict([(i.method_name, []) for i in methods])
      for obj in MethodRelation.objects.all():
          parent_name = obj.method.method_name
          child_name = obj.child.method_name
          res[parent_name].append(child_name)
      return res
    
  def get_parents(self):
    methods = Method.objects.all()
    res = dict([(i.method_name, []) for i in methods])
    for obj in MethodRelation.objects.all():
        parent_name = obj.method.method_name
        child_name = obj.child.method_name
        res[child_name].append(parent_name)
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