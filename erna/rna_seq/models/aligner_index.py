'''
reference index used by index
'''
import json
import os
from typing import Iterable
from django.db import models
from django.conf import settings
from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType

import rna_seq.models
from .tool import Tool
# from .annotation import Annotation
from pipelines.utils.dir import Dir


class AlignerIndexManager(models.Manager):
  # def refresh_annotation(self):
  #   res = []
  #   annotations = Annotation.objects.filter(file_format='fna', annot_type='genomic')
  #   tools = Tool.objects.filter(exe_name__contains='-build')
  #   for annot in annotations:
  #     fa_path = annot.file_path
  #     index_dir_path = os.path.join(os.path.dirname(fa_path), 'index')
  #     if os.path.isdir(index_dir_path):
  #       for tool in tools:
  #         index_path = os.path.join(index_dir_path, \
  #           f"{tool.tool_name}_{tool.version}_")
  #         if len(os.listdir(index_dir_path)) > 0:
  #           obj = self.load_reference(tool, annot, index_path)
  #           res.append(obj)
  #   return res
  
  def scan_index_dir(self) -> Iterable:
    index_dir = settings.INDEX_DIR
    for name in os.listdir(index_dir):
      infile = os.path.join(index_dir, name, 'info.json')
      if os.path.isfile(infile):
        with open(infile, 'r') as f:
          yield (int(name), json.load(f))

  def refresh(self):
    '''
    update references if index is built given one aligner
    '''
    # delete all
    self.all().delete()

    res = []
    indexes = self.scan_index_dir()
    for digit_name, info in indexes:
      this_model = getattr(rna_seq.models, info['model_name'])
      defaults = {
        'tool': Tool.objects.get(**info['tool']),
        'index_path': info['index_path'],
        'content_object': this_model.objects.get(**info['model']),
      }
      obj = self.update_or_create(id=digit_name, defaults=defaults)
      res.append(obj)
    return res

class AlignerIndex(models.Model):
  tool = models.ForeignKey(
    'rna_seq.Tool',
    on_delete=models.CASCADE,
  )
  index_path = models.CharField(
    max_length=256,
    verbose_name= 'index path used by aligner',
  )
  # related models: Annotation, RNA
  # The model must define a certain fa_path
  content_type = models.ForeignKey(
    ContentType,
    on_delete=models.CASCADE,
  )
  object_id = models.PositiveIntegerField()
  content_object = GenericForeignKey('content_type', 'object_id')

  objects = AlignerIndexManager()

  class Meta:
    app_label = "rna_seq"
    ordering = ["tool", "index_path"]

