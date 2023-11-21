'''
reference index used by index
'''
import os
from django.db import models
from django.conf import settings

from .annotation import Annotation
from .genome import Genome
from .tool import Tool

class ReferenceManager(models.Manager):
  def load_reference(self, tool:Tool, annotation:Annotation, index_path):
    data = {'index_path': index_path,}
    res = self.update_or_create(
      tool=tool,
      annotation=annotation,
      defaults = data
    )
    return res

  def refresh(self):
    '''
    update references if index is built given one aligner
    '''
    # delete all
    self.all().delete()

    res = []
    annotations = Annotation.objects.filter(file_format='fna', annot_type='genomic')
    tools = Tool.objects.filter(exe_name__contains='-build')
    for annot in annotations:
      fa_path = annot.file_path
      index_dir_path = os.path.join(os.path.dirname(fa_path), 'index')
      for tool in tools:
        index_path = os.path.join(index_dir_path, \
          f"{tool.tool_name}_{tool.version}_")
        if len(os.listdir(index_dir_path)) > 0:
          obj = self.load_reference(tool, annot, index_path)
          res.append(obj)
    return res
  
  # def refresh(self):
  #   '''
  #   update references if index is built given one aligner
  #   '''
  #   res = []
  #   ref_dir = settings.REFERENCES_DIR
  #   for data_source in os.listdir(ref_dir):
  #     genome_dir = os.path.join(ref_dir, data_source, 'genome')
  #     for specie in os.listdir(genome_dir):
  #       specie_dir = os.path.join(genome_dir, specie)
  #       for version in os.listdir(specie_dir):
  #         genome = Genome.objects.get_genome(specie, version)
  #         index_path = os.path.join(specie_dir, version, 'index')
  #         index_files = os.listdir(index_path)
  #         aligner_names = [i.split('.', 1)[0] for i in index_files]
  #         for aligner in list(set(aligner_names)):
  #           defaults = {
  #             'index_path': index_path,
  #           }
  #           obj = self.update_or_create(
  #             genome=genome,
  #             aligner=aligner,
  #             defaults = defaults)
  #           res.append(obj)
  #   return res

class Reference(models.Model):
  tool = models.ForeignKey(
    Tool,
    on_delete=models.CASCADE
  )
  annotation = models.ForeignKey(
    Annotation,
    on_delete=models.CASCADE,
    related_name = 'references',
  )
  index_path = models.CharField(max_length=256)

  objects = ReferenceManager()

  class Meta:
    app_label = "rna_seq"
    unique_together = ["tool", "annotation"]
    ordering = ["tool", "annotation"]
