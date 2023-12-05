'''
sequence alignment
'''
import os
from django.conf import settings

import rna_seq.models
from rna_seq.models import Genome, Reference, Tool
from utils.dir import Dir
from .process import Process
from .process_cmd import ProcessCMD

EXTERNALS_DIR = getattr(settings, 'EXTERNALS_DIR')
REFERENCES_DIR = getattr(settings, 'REFERENCES_DIR')


class Align:
  def __init__(self, params:dict=None):
    self.params = params if params else {}

  '''
  build index
  '''
  def build_index(self):
    '''
    method: build_index
    '''
    tool = self.params.get('tool')
    if tool is None:
      return None
    
    task_params = self.params['task'].get_params()
    this_model = getattr(rna_seq.models, task_params['model'])
    obj = this_model.objects.get(**task_params['query'])
    index_dir_path, index_path = obj.index_path(tool.tool_name, tool.version)
    Dir(index_dir_path).init_dir()

    input_data = {
      'fa_path': obj.fa_path,
      'index_path': index_path,
    } 
    self.params['cmd'], output_data = ProcessCMD.aligner_build_index(tool, input_data)
    if self.no_index_files(index_dir_path):
      Process.run_subprocess(self.params)
    self.params['output'].append(output_data)
    return None

  # used by erna_app.py
  def index_builder(self, specie:str, genome_version:str, \
      aligner:str, aligner_version:str) -> None:
    '''
    update self.params for building index
    if none is passed to self.params
    '''
    #get Genome
    genome = Genome.objects.get_genome(specie, genome_version)
    self.params['genome'] = genome

    # get Tool
    aligner_map = {
      'bowtie': 'bowtie2-build',
      'hisat2': 'hisat2-build',
      'STAR': 'STAR',
    }
    if aligner_map.get(aligner):
      self.params['tool'] = Tool.objects.get(
        tool_name=aligner,
        version=aligner_version,
        exe_name=aligner_map[aligner]
      )
    return None

  def no_index_files(self, path:str) -> bool:
    '''
    check if index files exist
    aligner: bowtie, bowtie2, hisat2
    '''
    tool_name = self.params['tool'].tool_name
    index_files = [i for i in os.listdir(path) if tool_name in i]
    if len(index_files) >= 4:
      return False
    return True

  def get_index_path(self):
    '''
    get index for aligner
    '''
    # Firstly check Task.params
    if self.params.get('task') and self.params['task'].get_params():
      return self.params['task'].get_params().get('index_path')
    
    # secondly check its parent TaskExecution
    if self.params.get('parent_params'):
      parent_output = self.params['parent_params']['output']
      for item in parent_output:
        if 'index_path' in item:
          return item['index_path']

    # Finally check its execution of parent Task
    for parent_output in self.params['parent_outputs']:
      if 'index_path' in parent_output:
        return parent_output['index_path']
    return None


  def build_genome_index(self):
    '''
    build index for genome alignment
    '''
    tool = self.params.get('tool')
    if tool is None:
      return None
    
    # build index given a specific file
    if self.params['annot_genomic_dna']:
      fa_path = self.params['annot_genomic_dna'].file_path
      index_dir_path = os.path.join(os.path.dirname(fa_path), 'index')
      index_path = os.path.join(index_dir_path, f"{tool.tool_name}_{tool.version}_")
      Dir(index_dir_path).init_dir()
      gtf_path = self.params['annot_genomic_gtf'].file_path if \
        self.params.get('annot_genomic_gtf') else ''

      input_data = {
        'index_path': index_path,
        'fa_path': fa_path,
        'gtf_path': gtf_path,
      }
      cmd, output_data = [], {}
      if tool.tool_name == 'star':
        cmd, output_data = ProcessCMD.star_build_index(tool, input_data)
      else:
        cmd, output_data = ProcessCMD.aligner_build_index(tool, input_data)
      self.params['cmd'] = cmd

      # skip building if index files exist
      if self.no_index_files(index_dir_path):
        Process.run_subprocess(self.params)

      # update annot.Reference
      Reference.objects.load_reference(tool, self.params['annot_genomic_dna'], index_path)
      
      # update Task
      self.params['output'].append(output_data)
      for child_task in self.params['children']:
        child_task.update_params(output_data)
    return None


  '''
  sequence alignment
  '''
  def align(self):
    '''
    '''
    tool = self.params.get('tool')
    if tool is None:
      return None
    
    sample_files = self.params['parent_outputs']
    for input_data in sample_files:
      sample_name = input_data['sample_name']
      output_prefix = os.path.join(self.params['output_dir'], sample_name)
      input_data['index_path'] = self.get_index_path()
      input_data['output_prefix'] = output_prefix,

      cmd, output_data = [], {}
      exe_name = self.params['tool'].exe_name
      if exe_name == 'hisat2':
        cmd, output_data = ProcessCMD.hisat2_align(tool, input_data)
      elif exe_name == 'bowtie2':
        cmd, output_data = ProcessCMD.bowtie2_align(tool, input_data)
      elif exe_name == 'STAR':
        cmd, output_data = ProcessCMD.star_align(tool, input_data)
      self.params['cmd'] = cmd

      # run process
      self.params['force_run'] = False if os.path.isfile(output_data['sam_file']) else True
      Process.run_subprocess(self.params)

      # update params['output']
      self.params['output_prefix'] = output_prefix
      self.params['output'].append(output_data)
    return None


