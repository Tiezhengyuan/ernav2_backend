'''
retrieve data for further analysis
'''
import os
from biofile import GTF, GFF
from biosequtils import Dir, KeyValue

from rna_seq.models import SampleProject
from .process import Process


class Collect:
  def __init__(self, params:dict):
    self.params = params

  def import_data(self):
    '''
    launched by task T00
    '''
    # get same ~ raw data
    self.import_sample_data()

  def import_sample_data(self):
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
      KeyValue.init_dict(res, [sample_name, file_type], [])
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
      ' '.join([i['stringtie_gtf_file'] for i in outputs]),
    ]
    self.params['force_run'] = False if os.path.isfile(merged_gtf_file) else True

    # update output
    self.params['output'].append({
      'cmd': ' '.join(self.params['cmd']),
      'annotation_file': annotation_file,
      'merged_transcripts': merged_gtf_file,
    })


  @staticmethod
  def import_gtf_annotations(annot_file):
    res = {}
    file_type = annot_file[-3:]
    outdir =  os.path.join(os.path.dirname(annot_file), f'{file_type}_features')
    if os.path.isdir(outdir):
      for file_name in os.listdir(outdir):
        if file_name.endswith('json'):
          feature, _ = os.path.splitext(file_name)
          res[feature] = os.path.join(outdir, file_name)
    else:
      Dir(outdir).init_dir()
      res = GTF(annot_file, outdir).split_by_feature()
    return res
  
  @staticmethod
  def import_gff_annotations(annot_file):
    res = {}
    file_type = annot_file[-3:]
    outdir =  os.path.join(os.path.dirname(annot_file), f'{file_type}_features')
    if os.path.isdir(outdir):
      for file_name in os.listdir(outdir):
        if file_name.endswith('json'):
          feature, _ = os.path.splitext(file_name)
          res[feature] = os.path.join(outdir, file_name)
    else:
      Dir(outdir).init_dir()
      res = GFF(annot_file, outdir).split_by_feature()
    return res


      
