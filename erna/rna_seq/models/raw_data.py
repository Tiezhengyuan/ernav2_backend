'''
one file one record.
'''
from django.db import models
from django.conf import settings
import os

STANDER_FORMAT = {
  'FQ': 'FSTQ',
  'FASTQ': 'FASTQ',
  'FA': 'FASTA', 
  'FASTA': 'FASTA',
  'SAM': 'SAM',
  'BAM': 'BAM',
  'GBK': 'GBK',
  'GTF': 'GTF',
  'GFF': 'GFF',
  'GFF3': 'GFF3',
}

class RawDataManager(models.Manager):
  def get_batch_files(self, batch_name:str):
    return self.filter(batch_name=batch_name)

  def detect_file_format(self, file_name:str):
    split_tup = os.path.splitext(file_name)
    file_type = split_tup[1][1:].upper() if len(split_tup) > 1 else ''
    return STANDER_FORMAT.get(file_type, "UN")
  
  def detect_file_type(self, file_name:str, file_format:str=None):
    if file_format == 'FASTQ':
      if 'R2' in file_name:
        return 'R2'
      return 'R1'
    return 'UN'
  
  def add_data(self, batch_name:str, file_path:str, file_name:str):
    file_format = self.detect_file_format(file_name)
    file_type = self.detect_file_type(file_name, file_format)
    res = self.create(batch_name=batch_name, file_path=file_path, \
      file_name=file_name, file_type=file_type, file_format=file_format)
    return res

  def load_data(self, raw_data:dict):
    res = []
    for batch_name in raw_data:
      for file_path, file_name in raw_data[batch_name]:
        obj = self.add_data(batch_name, file_path, file_name)
        res.append(obj)
    return res


class RawData(models.Model):
  file_path = models.CharField(max_length=512)
  file_name = models.CharField(max_length=128)
  file_format = models.CharField(
    max_length=10,
    blank=True,
    null=True,
  )
  file_type = models.CharField(
    max_length=10,
    blank=True,
    null=True
  )
  batch_name = models.CharField(
    max_length=20,
    blank=True,
    null=True
  )
  parsed = models.BooleanField(default=False)

  objects = RawDataManager()

  class Meta:
    app_label = 'rna_seq'
    unique_together = ['batch_name', 'file_path', 'file_name']
    ordering = ['batch_name', 'file_path', 'file_name']
  
  def to_dict(self):
    return {
      'raw_data_id': self.id,
      'batch_name': self.batch_name,
      'file_name': self.file_name,
      'full_path': self.full_file_path,
    }

  @property
  def storage_path(self):
    return settings.DATA_DIR

  @property
  def full_file_path(self):
    return os.path.join(self.file_path, self.file_name)

  def file_exists(self) -> bool:
    full_path = self.full_file_path()
    if os.path.isfile(full_path):
      return True
    return False