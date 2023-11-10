'''
exact match pattern:
shift-and:
shift-or:
'''
from sequence.sequence import Sequence

class Shift:
  def __init__(self, target_seq:str):
    self.pattern, _ = Sequence(target_seq).nt_bitvector()
    self.pattern_len = len(target_seq)
  
  def format_code(self, code:int):
    return format(code, f"0{self.pattern_len}b")
  
  def shift_and(self, text:str, first_occurrence=False):
    '''
    Shift-and
    '''
    p, res = 0, []
    target = 2 ** (self.pattern_len -1)
    for i, v in enumerate(text):
      v_code = self.pattern.get(v, 0)
      p = (p << 1 | 1) & v_code
      # print(self.format_code(v_code), self.format_code(p))
      if p >= target:
        res.append(i - self.pattern_len + 1)
        if first_occurrence:
          return res
    return res