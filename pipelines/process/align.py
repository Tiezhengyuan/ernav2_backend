'''
sequence alignment
'''
import os
import subprocess
from django.conf import settings

from rna_seq.models import Genome, Reference, Tool, Annotation
from utils.dir import Dir

EXTERNALS_DIR = getattr(settings, 'EXTERNALS_DIR')
REFERENCES_DIR = getattr(settings, 'REFERENCES_DIR')


class Align:
  def __init__(self, params:dict=None):
    self.params = params if params else {}

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
    exe_name = None
    if aligner == 'bowtie':
      exe_name='bowtie2-build'
    elif aligner == 'hisat2':
      exe_name='hisat2-build'
    if exe_name:
      self.params['tool'] = Tool.objects.get(
        tool_name=aligner,
        version=aligner_version,
        exe_name=exe_name
      )
    return None

  def no_index_files(self, index_dir_path:str) -> bool:
    '''
    check if index files exist
    '''
    tool_name = self.params['tool'].tool_name
    index_files = [ name for name in os.listdir(index_dir_path) \
      if tool_name in name]
    print(index_files)
    if len(index_files) >= 4:
      return False
    return True

  
  def build_index(self):
    '''
    build index for genome alignment
    '''
    tool = self.params.get('tool')
    if tool is None:
      return None

    annotations = Annotation.objects.filter(genome=self.params['genome'])
    for annot in annotations:
      # build index given a specific file
      if annot.annot_type == 'genomic' and annot.file_format == 'fna':
        fa_path = annot.file_path
        index_dir_path = os.path.join(os.path.dirname(fa_path), 'index')
        index_path = os.path.join(index_dir_path, \
          f"{tool.tool_name}_{tool.version}_")
        Dir(index_dir_path).init_dir()
                
        # skip building if index files exist
        if self.no_index_files(index_dir_path):
          cmd = [tool.exe_path, fa_path, index_path,]
          print(cmd)
          res = subprocess.run(cmd, capture_output=True, text=True)
          print(res.stdout)
          print(res.stderr)

        # update annot.Reference
        return Reference.objects.load_reference(
          tool, annot, index_path)
    return None



    
  def align_transcriptome(self):
    '''
    '''
    sample_files = self.project_sample_files()
    cmd = None
    if self.params['tool']['exe_name'] == 'bowtie2':
      cmd = self.cmd_bowtie()
    elif self.params['tool']['exe_name'] == 'hisat2':
      cmd = self.cmd_hisat2()

    if cmd:
      res = subprocess.run(cmd, capture_output=True, text=True)
      print(res.stdout)
      print(res.stderr)
      return res
    return None

  def project_sample_files(self):
    '''
    retrieve fastq files
    {
      'sample_name': {
        'R1': 'sample_R1.fq',
        'R2': 'sample_R2.fq',
      },
    }
    '''
    res = {}
    for sample_file in self.params['sample_files']:
      sample_name = sample_file.sample.sample_name
      if sample_name not in res:
        res[sample_name] = {}
      file_type = sample_file.raw_data.file_type
      if file_type not in res[sample_name]:
        res[sample_name][file_type] = []
      file = os.path.join(sample_file.raw_data.file_path, \
        sample_file.raw_data.file_name)
      if file not in res[sample_name][file_type]:
        res[sample_name][file_type].append(file)
    return res
    

  def cmd_hisat2(self, input_data:dict):
    '''

    '''
    cmd = [
      self.params['tool']['exe_path'],
      '-x', params['index_path'],
      '-S', os.path.join(input_data['sam_file']),
    ]
    if input_data['R1']:
      cmd.append(f"-1 {input_data['R1']}")
    if input_data['R2']:
      cmd.append(f"-2 {input_data['R2']}")
    if input_data['bam']:
      cmd.append(f"-b {input_data['bam']}")
    return cmd


  def cmd_bowtie(self, params:dict, input_data:dict):
    '''

    '''
    cmd = [
      os.path.join(EXTERNALS_DIR, 'bowtie2'),
      '-x', params['index_path'],
      '-S', os.path.join(input_data['sam_file']),
    ]
    if input_data['R1']:
      cmd.append(f"-1 {input_data['R1']}")
    if input_data['R2']:
      cmd.append(f"-2 {input_data['R2']}")
    if input_data['bam']:
      cmd.append(f"-b {input_data['bam']}")
    return cmd

