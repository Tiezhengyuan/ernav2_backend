'''
'''
import gzip
from Bio import SeqIO
from typing import Iterable


class FASTQ:
  def __init__(self, infile:str):
    self.infile = infile
  
  def parse_records(self)->Iterable:
    if self.infile.endswith('.gz'):
      with gzip.open(self.infile, 'rt') as f:
        for rec in SeqIO.parse(f, 'fastq'):
          yield rec
    else:
      with open(self.infile, 'r') as f:
        for rec in SeqIO.parse(f, 'fastq'):
          yield rec

