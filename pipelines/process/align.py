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
  def __init__(self, params:dict):
    self.parms = params

  def build_index(self, data_source:str, specie:str, version:str):
    '''
    build index for genome alignment
    '''
    genome = Genome.objects.get(data_source=data_source, \
      specie=specie, version=version)
    annotations = Annotation.objects.filter(genome=genome)
    tool = self.index_builder()
    if tool is None:
      return None

    for annot in annotations:
      fa_path = annot.file_path
      index_dir_path = os.path.join(os.path.dirname(fa_path), 'index')
      index_path = os.path.join(index_dir_path, \
        f"{tool.tool_name}_{tool.version}_")
      if len(os.listdir(index_dir_path)) > 0:
        return Reference.objects.load_reference(
          tool, annot, index_path)

      # build index
      if 'genomic.fna' in fa_path and '_from_' not in fa_path:
        Dir(index_dir_path).init_dir()
        cmd = [self.index_builder(), fa_path, index_path,]
        print(cmd)
        res = subprocess.run(cmd, capture_output=True, text=True)
        print(res.stdout)
        print(res.stderr)
        # update annot.Reference
        return Reference.objects.load_reference(
          tool, annot, index_path)
    return None

  def index_builder(self) -> str:
    '''
    update db.Annotation based on db.Genome
    '''
    exe_name = None
    aligner = self.params['tool']['tool_name']
    if aligner == 'bowtie':
      exe_name='bowtie2-build'
    elif aligner == 'hisat2':
      exe_name='hisat2-build'
    if exe_name:
      return Tool.objects.get(
        tool_name=aligner,
        version=self.params['tool']['version'],
        exe_name=exe_name
      )
    return None

    
  def genome_alignment(self, params:dict=None):
    '''
    genome alignment
    example args: params = {
      'project_id': 'P00001',
    }
    '''
    sample_files = self.project_sample_files(params['project_id'])
    cmd = self.cmd_bowtie()
    # res = subprocess.run(cmd, capture_output=True, text=True)
    # print(res.stdout)
    # print(res.stderr)
    # return res

  def project_sample_files(self, project_id:str):
    '''
    retrieve fastq files
    '''
    



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

