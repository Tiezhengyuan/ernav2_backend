'''
trim adapter for miRNA-seq
'''
from Bio import SeqIO
from biofile import FASTQ
from eseq import TrimSeq
import json
import os

from rna_seq.constants import TRIM

class TrimAdapter:
  def __init__(self, params:dict):
    self.params = params

  def __call__(self):
    stat = []
    task_params = self.params['task'].get_params()
    for item in self.params['parent_outputs']:
      # R1
      R1_files = []
      for infile in item.get('R1', []):
        outfile = os.path.join(self.params['output_dir'], os.path.basename(infile))
        R1_files.append(outfile)
        if os.path.isfile(outfile):
          break
        info = {
          'sample_name': item['sample_name'],
          'outfile': outfile,
        }
        res = {}
        if 'adapter_3end' in task_params:
          res = self.trim_3end_adapter(task_params, infile, outfile)
        elif 'keep_5end' in task_params:
          res = self.keep_5end(task_params, infile, outfile)
        info.update(res)
        stat.append(info)

      # output
      output = {'sample_name': item['sample_name'],}
      if R1_files:
        output['R1'] = R1_files
      self.params['output'].append(output)
    with open(os.path.join(self.params['output_dir'], 'stat.json'), 'w') as f:
      json.dump(stat, f, indent=4)

  def trim_3end_adapter(self, task_params, infile:str, outfile:str):
    '''
    trim adapter sequence at 3-end given adapter sequence
    '''
    info = {'total_reads': 0, 'trimmed_reads': 0, 'short_trimmed_reads': 0}
    trimmer = TrimSeq(
      task_params['adapter_3end'],
      task_params.get('min_match'),
      '3end',
      task_params.get('max_err')
    )
    # read fastq
    with open(outfile, 'w') as out_handle:
      fq_iter = FASTQ().parse_records(infile)
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

    total = info['total_reads'] if info['total_reads'] > 0 else 1
    info['trim_percentage'] = round(info['trimmed_reads']/total, 4)
    print(info)
    return info
  
  def keep_5end(self, task_params, infile, outfile):
    '''
    trim fixed length from 3-end
    '''
    info = {'total_reads': 0, 'trimmed_reads': 0}
    pos = TRIM['min_len'] # minimum kept length
    try:
      pos = int(task_params.get('keep_5end'))
    except Exception as e:
      pass
      
    # trim
    with open(outfile, 'w') as out_handle:
      fq_iter = FASTQ().parse_records(infile)
      for rec in fq_iter:
        info['total_reads'] += 1
        if pos >= TRIM['min_len']:
          rec = rec[:pos+1]
          info['trimmed_reads'] += 1
        SeqIO.write(rec, out_handle, 'fastq')
    return info