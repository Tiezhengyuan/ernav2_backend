'''
assembly genome or transcriptome
'''
import os
from .process import Process
from .process_cmd import ProcessCMD

class Assemble:
  def __init__(self, params:dict):
    self.params = params
  
  def assemble_transcripts(self):
    '''
    run method: assemble_transcripts
    '''
    tool = self.params.get('tool')
    if tool is None:
      return None
    
    for parent_output in self.params['parent_outputs']:
      # prepare commands
      sample_name = parent_output.get('sample_name', '_')
      output_prefix = os.path.join(self.params['output_dir'], sample_name)
      input_data = {
        'sample_name': sample_name,
        'sorted_bam_file': parent_output['sorted_bam_file'],
        'output_prefix': output_prefix,
        'annotation_file': self.params['annot_genomic_gtf'].file_path,
      }
      # assemble transcripts
      if self.params['tool'].tool_name == 'stringtie':
        input_data['stringtie_gtf_file'] = output_prefix + '.gtf'
        if not os.path.isfile(input_data['stringtie_gtf_file']):
          self.params['cmd'] = ProcessCMD.stringtie_assemble(tool, input_data)
          Process.run_subprocess(self.params)
      # update output
      self.params['output'].append(input_data)
    return None


  # def annotation_file(self):
  #   '''
  #   Use a reference annotation file (in GTF or GFF3 format)
  #   to guide the assembly process. 
  #   '''
  #   # firstly try task.params
  #   task_params = self.params['task'].get_params()
  #   if task_params.get('annotation_file'):
  #     return task_params['annotation_file']
    
  #   # secondly try Annotation with genome
  #   if self.params['annot_genomic_gtf']:
  #     self.params['task'].update_params({
  #       'annotation_file': self.params['annot_genomic_gtf'].file_path,
  #     })
  #     return self.params['annot_genomic_gtf'].file_path
  #   return None


