'''
process raw data namely fastq
'''
import os
from typing import Iterable
from sample.models import RawData as RawDataModel

class RawData:
  def __init__(self, dir_raw_data):
    self.dir_raw_data = dir_raw_data

  def scan_raw_data(self):
    raw_data = {'UN': []}
    for batch in os.listdir(self.dir_raw_data):
      path = os.path.join(self.dir_raw_data, batch)
      if os.path.isdir(path):
        raw_data[batch] = []
        for root, dirs, files in os.walk(path, topdown=False):
          for file in files:
            raw_data[batch].append((root, file))
      else:
        raw_data['UN'].append((self.dir_raw_data, batch))
    return RawDataModel.objects.load_data(raw_data)
