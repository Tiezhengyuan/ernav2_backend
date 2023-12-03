'''
sequence alignment
'''
import os
from .process import Process
from django.conf import settings

import rna_seq.models
from rna_seq.models import Genome, Reference, Tool
from utils.dir import Dir

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
    task_params = self.params['task'].get_params()
    this_model = getattr(rna_seq.models, task_params['model'])
    obj = this_model.objects.get(**task_params['query'])
    index_dir_path, index_path = obj.index_path(
      self.params['tool'].tool_name,
      self.params['tool'].version
    )
    Dir(index_dir_path).init_dir()
    self.params['cmd'] = [
      self.params['tool'].exe_path,
      obj.fa_path,
      index_path
    ]
    if self.no_index_files(index_dir_path):
      Process.run_subprocess(self.params)
    output = {
      'cmd': ' '.join(self.params['cmd']),
      'index_path': index_path
    }
    self.params['output'].append(output)

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

    for annot in self.params['annotations']:
      # build index given a specific file
      if annot.annot_type == 'genomic' and annot.file_format == 'fna':
        fa_path = annot.file_path
        index_dir_path = os.path.join(os.path.dirname(fa_path), 'index')
        index_path = os.path.join(index_dir_path, f"{tool.tool_name}_{tool.version}_")
        Dir(index_dir_path).init_dir()

        # skip building if index files exist
        self.params['cmd'] = [tool.exe_path, fa_path, index_path,]
        if self.no_index_files(index_dir_path):
          Process.run_subprocess(self.params)

        # update annot.Reference
        Reference.objects.load_reference(tool, annot, index_path)
        
        # update Task
        output = {
          'cmd': ' '.join(self.params['cmd']),
          'index_path': index_path
        }
        self.params['output'].append(output)
        for child_task in self.params['children']:
          child_task.update_params(output)
    return None


  '''
  sequence alignment
  '''
  def align(self):
    '''
    '''
    sample_files = self.params['parent_outputs']
    for input_data in sample_files:
      exe_name = self.params['tool'].exe_name
      if exe_name == 'hisat2':
        self.cmd_hisat2(input_data)
      elif exe_name == 'bowtie2':
        self.cmd_bowtie2(input_data)

      # run process
      Process.run_subprocess(self.params)
    return None
 
  def cmd_hisat2(self, input_data:dict):
    '''
    format cmd
    update self.params
    '''
    cmd = [
      self.params['tool'].exe_path,
      '-x', self.get_index_path(),
    ]
    if input_data.get('bam'):
      cmd += ['-b', input_data['bam']]
    else:
      if input_data.get('R1'):
        cmd += ['-1', ','.join(input_data['R1'])]
      if input_data.get('R2'):
        cmd += ['-2', ','.join(input_data['R2'])]

    sample_name = input_data['sample_name']
    output_prefix = os.path.join(self.params['output_dir'], sample_name)
    sam_file = f"{output_prefix}.sam"
    cmd += ['-S', sam_file]
    self.params['cmd'] = cmd
    self.params['output_prefix'] = output_prefix
    self.params['output'].append({
      'cmd': ' '.join(cmd),
      'sample_name': sample_name,
      'output_prefix': output_prefix,
      'sam_file': sam_file,
    })
    self.params['force_run'] = False if os.path.isfile(sam_file) else True
    return cmd


  def cmd_bowtie2(self, input_data:dict):
    cmd = [
      self.params['tool'].exe_path,
      '-x', self.get_index_path(),
    ]
    if input_data.get('R1') and input_data.get('R2'):
      cmd += [
        '-1', ','.join(input_data['R1']),
        '-2', ','.join(input_data['R2']),
      ]
    elif input_data.get('bam'):
      cmd += ['-b', input_data['bam']]
    elif input_data.get('unaligned'):
      cmd += ['-f', input_data['unaligned']]
    else:
      raw_data = input_data.get('R1', []) + input_data.get('R2', [])
      cmd += ['-U', ','.join(raw_data)]

    sample_name = input_data['sample_name']
    output_prefix = os.path.join(self.params['output_dir'], sample_name)
    sam_file = f"{output_prefix}.sam"
    cmd += ['-S', sam_file]
    self.params['cmd'] = cmd
    self.params['output_prefix'] = output_prefix
    self.params['output'].append({
      'cmd': ' '.join(cmd),
      'sample_name': sample_name,
      'output_prefix': output_prefix,
      'sam_file': sam_file,
    })
    self.params['force_run'] = False if os.path.isfile(sam_file) else True
    return cmd