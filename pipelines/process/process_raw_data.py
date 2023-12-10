'''
process raw data namely fastq
'''
import os
from typing import Iterable
from django.conf import settings
from django.core import serializers

from rna_seq.models import RawData, Sample, SampleFile, SampleProject
from .process import Process

class ProcessRawData:

  def uncompress_raw_data(self):
    '''
    uncompress *.gz
    '''
    for rawdata_dir in settings.RAW_DATA_DIRS:
      for root, dirs, files in os.walk(rawdata_dir):
        for file in files:
          if file.endswith('.gz'):
            gz_file = os.path.join(root, file)
            Process.uncompress_gz(gz_file)

  def scan_raw_data(self):
    '''
    batch_name ~ (path, filename)
    '''
    raw_data = {'UN': []}
    for rawdata_dir in settings.RAW_DATA_DIRS:
      for batch in os.listdir(rawdata_dir):
        path = os.path.join(rawdata_dir, batch)
        if os.path.isdir(path):
          raw_data[batch] = []
          for root, dirs, files in os.walk(path, topdown=False):
            for file in files:
              raw_data[batch].append((root, file))
        else:
          raw_data['UN'].append((rawdata_dir, batch))
    return raw_data
  
  # TODO: deprecated in the future
  def refresh_raw_data(self):
    '''
    refresh table RawDATA
    '''
    RawData.objects.all().delete()
    raw_data = self.scan_raw_data()
    data = RawData.objects.load_data(raw_data)
    return serializers.serialize('json', data)

  def reset_sample(self):
    '''
    Reset all tables defined in the app sample
    1. Delete all data in RawData, Sample, SampleFile, SampleProject
    2. RawData: refresh raw data
    '''
    res = {
      'deleted': {},
      'created': {},
    }
    for db in (RawData, Sample, SampleFile, SampleProject):
      queryset = db.objects.all()
      res['deleted'][db.__name__] = queryset.count()
      queryset.delete()
    # create
    raw_data = self.scan_raw_data()
    data = RawData.objects.load_data(raw_data)
    res['created'][RawData.__name__] = len(data)
    return res

  def parse_sample_data(self, study_name, prefix=None, postfix=None):
    '''
    one raw data match one sample
    one sample may match to 1-many raw data
    '''
    study_samples = Sample.objects.filter(study_name=study_name)
    unparsed_data = RawData.objects.filter(parsed=False)
    sample_data = []
    for raw_data in unparsed_data:
      file_name = raw_data.file_name
      for sample in study_samples:
        target = sample.sample_name
        if prefix:
          target = prefix + target
        if postfix:
          target += postfix
        if target in file_name:
          pair = {
            'sample_id': sample.id,
            'sample_name': sample.sample_name,
            'raw_data_id': raw_data.id,
            'file_name': raw_data.file_name,
            'batch_name': raw_data.batch_name,
          }
          sample_data.append(pair)
          break
    return sample_data

