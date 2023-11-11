'''
reference index used by index
'''
import os
from django.db import models
from django.conf import settings
from .genome import Genome

class ReferenceManager(models.Manager):
  def load_reference(self, genome, aligner, fa_path, index_path):
    data = {
      'fa_path': fa_path,
      'index_path': index_path,
    }
    res = self.update_or_create(
      genome=genome,
      aligner=aligner,
      defaults = data
    )
    return res
  
  def refresh(self):
    '''
    update references if index is built given one aligner
    '''
    res = []
    ref_dir = settings.REFERENCES_DIR
    for data_source in os.listdir(ref_dir):
      genome_dir = os.path.join(ref_dir, data_source, 'genome')
      for specie in os.listdir(genome_dir):
        specie_dir = os.path.join(genome_dir, specie)
        for version in os.listdir(specie_dir):
          genome = Genome.objects.get_genome(specie, version)
          index_path = os.path.join(specie_dir, version, 'index')
          index_files = os.listdir(index_path)
          aligner_names = [i.split('.', 1)[0] for i in index_files]
          for aligner in list(set(aligner_names)):
            defaults = {
              'index_path': index_path,
            }
            obj = self.update_or_create(
              genome=genome,
              aligner=aligner,
              defaults = defaults)
            res.append(obj)
    return res

class Reference(models.Model):
  genome = models.ForeignKey(
    Genome,
    related_name = 'references',
    on_delete=models.CASCADE
  )
  aligner = models.CharField(max_length= 10)
  index_path = models.CharField(max_length=256)
  fa_path = models.CharField(
    max_length=256,
    null=True,
    blank=True,
  )

  objects = ReferenceManager()

  class Meta:
    app_label = "rna_seq"
    unique_together = ["genome", "aligner"]
    ordering = ["genome", "aligner"]
