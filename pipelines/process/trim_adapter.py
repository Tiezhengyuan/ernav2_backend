'''
trim adapter for miRNA-seq
'''
from Bio import SeqIO
import os

from sequence.trim_seq import TrimSeq
from biofile.fastq import FASTQ

class TrimAdapter:
  def __init__(self, params:dict):
    self.params = params

  def __call__(self):
    task_params = self.params['task'].get_params()
    for parent_output in self.params['parent_outputs']:
      file_path = parent_output['file_path']
      outfile = os.path.join(self.params['output_dir'], os.path.basename(file_path))
      if 'adapter_3end' in task_params:
        self.trim_3end_adapter(task_params, file_path, outfile)
      self.params['output'].append({
        'sample_name': parent_output['sample_name'],
        'file_path': outfile,
      })
    print('skipp adapter trimming.')

  def trim_3end_adapter(self, task_params, infile:str, outfile:str):
    info = {'total_reads': 0, 'trimmed_reads': 0,}
    trimmer = TrimSeq(
      task_params['adapter_3end'],
      task_params.get('min_match'),
      '3end',
      task_params.get('max_err')
    )
    # read fastq
    with open(outfile, 'w') as out_handle:
      fq_iter = FASTQ(infile).parse_records()
      for rec in fq_iter:
        info['total_reads'] += 1
        trimmed_seq, pos = trimmer.trim_3end(rec.seq)
        if pos > -1:
          rec = rec[:pos]
          info['trimmed_reads'] += 1
        SeqIO.write(rec, out_handle, 'fastq')
    info['trim_percentage'] = round(info['trimmed_reads']/info['total_reads'], 4)
    print(info)
    return info
  
