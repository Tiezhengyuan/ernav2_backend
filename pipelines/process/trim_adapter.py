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
    for item in self.params['parent_outputs']:
      output = {'sample_name': item['sample_name']}
      for type in [i for i in item if i in ('R1', 'R2')]:
        output[type] = []
        for infile in item[type]:
          outfile = os.path.join(self.params['output_dir'], os.path.basename(infile))
          if 'adapter_3end' in task_params:
            self.trim_3end_adapter(task_params, infile, outfile)
          output[type].append(outfile)
      self.params['output'].append(output)
    print('skipp adapter trimming.')

  def trim_3end_adapter(self, task_params, infile:str, outfile:str):
    info = {'total_reads': 0, 'trimmed_reads': 0, 'short_trimmed_reads': 0}
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
          if len(rec.seq) >= task_params.get('min_len', 12):
            info['trimmed_reads'] += 1
            SeqIO.write(rec, out_handle, 'fastq')
          else:
            # doesn't export to trimmed file
            info['short_trimmed_reads'] += 1
        else:
          # not trimmed
          SeqIO.write(rec, out_handle, 'fastq')

    info['trim_percentage'] = round(info['trimmed_reads']/info['total_reads'], 4)
    print(info)
    return info
  
