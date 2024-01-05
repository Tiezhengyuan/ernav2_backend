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
      }
      # assemble transcripts
      if self.params['tool'].tool_name == 'stringtie':
        # add annotations from gene.json to seqdata
        input_data['annotation_file'] = self.params['genome_annot']['gff'].file_path
        cmd = ProcessCMD.stringtie_assemble(tool, input_data)
        if not os.path.isfile(input_data['stringtie_gtf_file']):
          self.params['cmd'] = cmd
          Process.run_subprocess(self.params)
      # update output
      self.params['output'].append(input_data)
    return None



