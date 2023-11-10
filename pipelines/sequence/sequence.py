'''
sequence: DNA, RNA, Protein
'''
from Bio.Seq import Seq

class Sequence:
  def __init__(self, seq):
    self.seq = Seq(seq)

  def nt_bitvector(self)->dict:
    '''
    key-value: reversed bitset of positions
    '''
    def bitset(seq, mask, mask_len):
      if not seq:
        return mask
      letter = seq[0]
      tag = 0
      for k in list(mask):
        if k == letter:
          mask[k] += 2 ** mask_len
          tag = 1
          break
      if tag == 0:
        mask[letter] = 2 ** mask_len
      return bitset(seq[1:], mask, mask_len + 1)

    #
    mask_len = len(self.seq)
    mask = bitset(self.seq, {}, 0)
    mask_str = dict([(k, format(v, f'0{mask_len}b')) for k,v in mask.items()])
    return mask, mask_str
