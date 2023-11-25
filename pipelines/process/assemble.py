'''
assembly genome or transcriptome
'''
import os
from .process import Process

class Assemble:
  def __init__(self, params:dict):
    self.params = params
  
  def assemble_transcripts(self):
    '''
    run method: assemble_transcripts
    '''
    for parent_output in self.params['parent_outputs']:
      # prepare commands
      if self.params['tool'].tool_name == 'stringtie':
        self.cmd_stringtie(parent_output)
      elif self.params['tool'].tool_name == 'cufflinks':
        self.cmd_cufflinks(parent_output)
      Process.run_subprocess(self.params)
    return None

  def cmd_stringtie(self, parent_output:dict):
    '''
    stringtie <in.bam ..>  [options]
    '''
    sample_name = parent_output.get('sample_name', '_')
    output_prefix = os.path.join(self.params['output_dir'], sample_name)
    gtf_file = output_prefix + '.gtf'
    self.params['cmd'] = [
      self.params['tool'].exe_path,
      parent_output['sorted_bam_file'],
      '-o', gtf_file,
    ]
    # append reference annotation file in gtf/gff
    annot_file = self.annotation_file()
    ballgown_file, abundance_file, covered_file = None, None, None
    if annot_file:
      ballgown_file = output_prefix + '.ctab'
      abundance_file = output_prefix + '.abund'
      covered_file = output_prefix + '.cov'
      self.params['cmd'] += [
        '-G', annot_file,
        '-e', '-b', ballgown_file,
        '-A', abundance_file,
        '-C', covered_file,
      ]

    # update output
    self.params['output'].append({
      'annotation_file': annot_file,
      'cmd': ' '.join(self.params['cmd']),
      'sample_name': sample_name,
      'output_prefix': output_prefix,
      'gtf_file': gtf_file,
      'ballgown_file': ballgown_file,
      'abundance_file': abundance_file,
      'covered_file': covered_file,
    })


  def cmd_cufflinks(self, parent_output_item:dict):
    pass

  def annotation_file(self):
    '''
    Use a reference annotation file (in GTF or GFF3 format)
    to guide the assembly process. 
    '''
    # firstly try task.params
    task_params = self.params['task'].get_params()
    if task_params.get('annotation_file'):
      return task_params['annotation_file']
    
    # secondly try Annotation with genome
    for annot in self.params['annotations']:
      if annot.annot_type == 'genomic' and annot.file_format in ('gff', 'gtf'):
        self.params['task'].update_params({
          'annotation_file': annot.file_path,
        })
        return annot.file_path
    return None


  def merge_transcripts(self):
    '''
    run method: merge transcripts
    '''
    if self.params['tool'].tool_name == 'stringtie':
      self.cmd_stringtie_merge()
    return Process.run_subprocess(self.params)

  def cmd_stringtie_merge(self):
    '''
    stringtie --merge [Options] { gtf_list | strg1.gtf ...}
    '''
    outputs = self.params['parent_outputs']
    annotation_file = outputs[0]['annotation_file']
    merged_gtf_file = os.path.join(self.params['output_dir'], 'merged_transcripts.gtf')
    self.params['cmd'] = [
      self.params['tool'].exe_path, '--merge',
      '-G', annotation_file,
      '-o', merged_gtf_file,
      ' '.join([i['gtf_file'] for i in outputs]),
    ]

    # update output
    self.params['output'].append({
      'cmd': ' '.join(self.params['cmd']),
      'annotation_file': annotation_file,
      'merged_transcripts': merged_gtf_file,
    })
    