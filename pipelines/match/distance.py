

class Distance:
  def __init__(self, seq1:str):
    self.seq1 = seq1
  
  def match_5end(self, seq2:str):
    '''
    match from 5'-end
    return distance score
    '''
    dist = 0
    end = len(seq2) if len(seq2) <= len(self.seq1) \
      else len(self.seq1)
    for i in range(0, end):
      if self.seq1[i] == seq2[i]:
        dist += 1
    return dist
