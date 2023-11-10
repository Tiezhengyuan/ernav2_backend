'''
trim adapter for miRNA-seq
'''
from Bio import SeqIO,SeqRecord, Seq


from sequence.trim_seq import TrimSeq
from biofile.fastq import FASTQ

class TrimAdapter:
  def __init__(self, params:dict):
    self.params = params

  def __call__(self):
    trim_end = self.params.get('trim_end', '3end')
    if  trim_end == '3end':
      return self.trim_3end_adapter()
    elif  trim_end == '5end':
      return self.trim_5end_adapter()
    else:
      print('skipp adapter trimming.')

  def trim_3end_adapter(self):
    info = {
      'total_reads': 0,
      'trimmed_reads': 0,
    }
    trimmer = TrimSeq(
      self.params['adapter'],
      self.params.get('min_match'),
      '3end',
      self.params.get('max_err')
    )
    # read fastq
    with open(self.params['output'], 'w') as out_handle:
      fq_iter = FASTQ(self.params['input']).parse_records()
      for rec in fq_iter:
        info['total_reads'] += 1
        trimmed_seq, pos = trimmer.trim_3end(rec.seq)
        if pos > -1:
          rec = rec[:pos]
          info['trimmed_reads'] += 1
        SeqIO.write(rec, out_handle, 'fastq')
    info['trim_percentage'] = round(info['trimmed_reads']\
      /info['total_reads'], 4)
    return info
  
    def trim_5end_adapter(self):
      pass