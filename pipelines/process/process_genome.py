'''
process genome annotations before analysis
'''
from biofile import Wrap
import os

from pipelines.connector.connect_ncbi import ConnectNCBI
from pipelines.utils.handle_json import HandleJson
from pipelines.utils.dir import Dir

from rna_seq.models import Genome, Annotation, MolecularAnnotation

class ProcessGenome:

  def __init__(self, data_source=None, specie=None, version=None):
    self.data_source = data_source
    self.specie = specie
    self.version = version

  def retrieve_assembly_summary(self):
    res = {}
    if self.data_source == "NCBI":
      # download assembly_summary.txt
      client = ConnectNCBI()
      text_files = client.download_assembly_summary()
      for text_file in text_files.values():
        HandleJson(text_file).from_text()
      #load samples into db.Specie
      res['specie'] = client.load_species()
      #load samples into db.Genome
      res['genome'] = client.load_genomes()
    return res

  def download_genome(self, overwrite:bool=None):
    '''
    genome DNA sequences
    '''
    local_files = None
    if self.data_source == "NCBI":
      # download data
      client = ConnectNCBI()
      local_path, local_files = client.download_genome(
        self.specie, self.version, overwrite
      )

    res = {}      
    # update database
    if local_files:
      res['local_files'] = local_files
      # update db.Genome 
      obj = Genome.objects.filter(specie=self.specie, version=self.version)
      obj.update(is_ready=True, local_path=local_path)
      # update db.Annotation
      Annotation.objects.load_annotations(obj[0], local_files)
    return res

  def molecular_annotation(self, overwrite:bool=None) -> list:
    '''
    retrieve annoations accordding to molecular type from fa and gtf/gff
    Save metadata to MolecularAnnotation
    '''
    genomes = Genome.objects.filter(local_path__isnull=False)
    for genome in genomes:
      local_files = [os.path.join(genome.local_path, name) for \
        name in os.listdir(genome.local_path)]
      outdir = os.path.join(genome.local_path, 'features')
      Dir(outdir).init_dir()

      # process genome annotations
      wrapper = Wrap(local_files, outdir)
      meta = wrapper.load_output()
      if overwrite or (not meta):
        meta = wrapper.ncbi_fa_gff()
        wrapper.save_output(meta, True)
      # load meta into MolecularAnnotation
      meta = wrapper.load_output()
      res = MolecularAnnotation.objects.load(meta)
      return res
      
      







