'''
retrieve data for further analysis
'''
from copy import deepcopy
import os
import pandas as pd
import pysam
from .process import Process

from rna_seq.models import SampleProject
from pipelines.utils.utils import Utils

class Collect:
  def __init__(self, params:dict):
    self.params = params

  def import_data(self):
    # Samples
    project_samples = SampleProject.objects.filter(
      project=self.params['project']
    )
    sample_files = [obj.sample_file for obj in project_samples]
    self.params['sample_files'] = sample_files
    res = {}
    for sf in sample_files:
      path = os.path.join(sf.raw_data.file_path, sf.raw_data.file_name)
      sample_name = sf.sample.sample_name
      file_type = sf.raw_data.file_type
      Utils.init_dict(res, [sample_name, file_type], [])
      res[sample_name][file_type].append(path)
    for k,v in res.items():
      v['sample_name'] = k 
      self.params['output'].append(v) 

  def merge_transcripts(self):
    '''
    run method: merge transcripts
    '''
    if self.params['tool'].tool_name == 'stringtie':
      self.cmd_stringtie_merge_transcripts()
    return Process.run_subprocess(self.params)

  def cmd_stringtie_merge_transcripts(self):
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
    self.params['force_run'] = False if os.path.isfile(merged_gtf_file) else True

    # update output
    self.params['output'].append({
      'cmd': ' '.join(self.params['cmd']),
      'annotation_file': annotation_file,
      'merged_transcripts': merged_gtf_file,
    })




      
