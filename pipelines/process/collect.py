'''
retrieve data for further analysis
'''
from copy import deepcopy
import pandas as pd
import os
from rna_seq.models import SampleProject
from utils.utils import Utils

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

  def count_reads(self):
    '''
    reads counting
    '''
    for parent in self.params['parents']:
      print(parent)
      output = parent.task_execution.get_output()
      if parent.method_tool.tool.tool_name == 'stringtie':
        self.stringtie_counting(output)
  

  def stringtie_counting(self, parent_output:dict):
    first = parent_output[0]
    df_tpm, df_fpkm = self.read_abund(first['abundance_file'], first['sample_name'])
    for item in parent_output[1:]:
      df1, df2 = self.read_abund(item['abundance_file'], item['sample_name'])
      df_tpm = pd.merge(df_tpm, df1, how='outer').fillna(0)
      df_fpkm = pd.merge(df_fpkm, df2, how='outer').fillna(0)
    #
    tpm_file = os.path.join(self.params['output_dir'], 'TPM.txt')
    df_tpm.to_csv(tpm_file, index=False, sep='\t')
    fpkm_file = os.path.join(self.params['output_dir'], 'FPKM.txt')
    df_fpkm.to_csv(fpkm_file, index=False, sep='\t')
    self.params['output'].append({
      'TPM': tpm_file,
      'FPKM': fpkm_file,
    })

  def read_abund(self, infile, sample_name):
    df=pd.read_csv(infile, sep='\t')
    # FPKM
    df1=df.loc[:, df.columns != 'TPM']
    df1.columns = df1.columns.str.replace('FPKM', sample_name)
    # TPM
    df2=df.loc[:, df.columns != 'FPKM']
    df2.columns = df2.columns.str.replace('TPM', sample_name)
    return df1, df2
