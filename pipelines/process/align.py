'''
sequence alignment
'''
import os
import subprocess
from django.conf import settings

from rna_seq.models import Genome, Reference
from utils.dir import Dir

EXTERNALS_DIR = getattr(settings, 'EXTERNALS_DIR')
REFERENCES_DIR = getattr(settings, 'REFERENCES_DIR')


class Align:
  def __init__(self, aligner:str):
    self.aligner = aligner

  def build_index(self, data_source:str, specie:str, version:str):
    '''
    build index for genome alignment
    '''
    genome = Genome.objects.get(data_source=data_source, \
      specie=specie, version=version)

    for annot in genome.annots.values():
      fa_path = annot['file_path']
      if 'genomic.fna' in fa_path and '_from_' not in fa_path:
        index_path = os.path.join(os.path.dirname(fa_path), 'index')
        Dir(index_path).init_dir()
        index_path += f"/{self.aligner}"
        # build index
        cmd = [self.index_builder(), fa_path, index_path,]
        print(cmd)
        res = subprocess.run(cmd, capture_output=True, text=True)
        print(res.stdout)
        print(res.stderr)
        # update annot.Reference
        res = Reference.objects.load_reference(genome,
          self.aligner, fa_path, index_path)
        return res
    return None

  def index_builder(self) -> str:
    '''
    update db.Annotation based on db.Genome
    '''
    if self.aligner == 'bowtie2':
      return os.path.join(EXTERNALS_DIR, 'bowtie2-build')
    if self.aligner == 'bowtie':
      return os.path.join(EXTERNALS_DIR, 'bowtie-build')
    if self.aligner == 'hisat2':
      return os.path.join(EXTERNALS_DIR, 'hisat2-build')

  def genome_alignment(self, params:dict=None):
    '''
    genome alignment
    '''
    cmd = self.aligner_cmd()
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout)
    print(res.stderr)
    return res

  def aligner_cmd(self, params:dict, input_data:dict):
    '''

    '''
    cmd = []
    if self.aligner == 'bowtie2':
      cmd += [
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

