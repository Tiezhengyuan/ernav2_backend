'''
K-mer
'''
from typing import Iterable

class Kmer:
  def __init__(self, seq:Iterable):
    self.seq = seq
    self.seq_len = len(seq)
  
  def kmer(self, k:int, start:int=None) -> Iterable:
    '''
    forward
    '''
    if k < 2:
      k = 2
    elif k > self.seq_len:
      k = self.seq_len
    if start is None or abs(start) >= self.seq_len:
      start = 0
    elif start < 0:
      start += self.seq_len
    
    # Note: must add 1 in range()
    for i in range(start, self.seq_len - k + 1):
      subseq = self.seq[i:i+k]
      yield subseq
  
