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
from pipelines.utils.dir import Dir

INFO_FILE = 'info.json'

class AlignerIndexManager(models.Manager):
  
  def scan_index_dir(self) -> Iterable:
    '''
    scan index dir
    '''
    index_dir = settings.INDEX_DIR
    for name in os.listdir(index_dir):
      infile = os.path.join(index_dir, name, INFO_FILE)
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
        'tool': Tool.objects.get(pk=info['tool_id']),
        'index_path': info['index_path'],
        'content_object': this_model.objects.get(**info['model_query']),
      }
      obj = self.update_or_create(id=digit_name, defaults=defaults)
      res.append(obj)
    return res

  def new_index(self, tool, related_obj):
    '''
    Note: index_path is not defined in record when that is created
    '''
    print(tool, related_obj)
    obj = self.create(tool=tool, content_object=related_obj)
    index_dir_path = os.path.join(settings.INDEX_DIR, str(obj.id))
    Dir(index_dir_path).init_dir()
    fa_path = related_obj.file_path
    file_name, _ = os.path.splitext(os.path.basename(fa_path))
    index_path = os.path.join(index_dir_path, file_name)
    return obj, index_dir_path, index_path


class AlignerIndex(models.Model):
  tool = models.ForeignKey(
    'rna_seq.Tool',
    on_delete=models.CASCADE,
  )
  index_path = models.CharField(
    max_length=256,
    null=True,
    blank=True,
    verbose_name= 'index path used by aligner',
  )
  # related models: Annotation, RNA, MolecularAnnotation
  # The model must define a certain fa_path
  content_type = models.ForeignKey(
    ContentType,
    on_delete=models.CASCADE,
  )
  # id of the object of the related model
  # ids could be identical but come from different models
  object_id = models.PositiveIntegerField()
  content_object = GenericForeignKey('content_type', 'object_id')

  objects = AlignerIndexManager()

  class Meta:
    app_label = "rna_seq"
    unique_together = ('tool', 'content_type')
    ordering = ("tool", "index_path")

  def update_index(self, meta_data:dict) -> None:
    self.index_path = meta_data['index_path']
    self.save()
    # save info.json
    outfile = os.path.join(meta_data['index_dir_path'], INFO_FILE)
    with open(outfile, 'w') as f:
      json.dump(meta_data, f, indent=4)