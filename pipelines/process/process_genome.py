'''
process genome annotations before analysis
'''
from biofile import Wrap
from bioomics import NCBI, ANATOMY_GROUPS
from biosequtils import Dir, HandleJson
from django.conf import settings
import os
from typing import Iterable

from rna_seq.models import Genome, Specie, Annotation, MolecularAnnotation


class ProcessGenome:

  def __init__(self, data_source=None, specie=None, version=None, overwrite:bool=None):
    self.data_source = data_source
    self.specie = specie
    self.version = version
    self.overwrite = True if overwrite else False
    self.ref_dir = getattr(settings, 'REFERENCES_DIR')

  def ncbi_assembly_summary(self, groups:list=None):
    # download assembly_summary.txt
    conn = NCBI(self.ref_dir, self.overwrite)
    output_dir, text_files = conn.download_assembly_summary(groups)
    for text_file in text_files.values():
      HandleJson(text_file).from_text()

    #load samples into Specie and Genome
    # Note: the loarding order matters
    res = {}
    for anatomy, obj_iter in self.scan_meta_json(output_dir):
      meta = self.load_specie_genome(obj_iter)
      res[anatomy] = meta
    return res

  def scan_meta_json(self, output_dir:str) -> Iterable:
    for antonomy in ANATOMY_GROUPS:
      json_file = os.path.join(output_dir, 'assembly_summary', antonomy, 'assembly_summary.json')
      if os.path.isfile(json_file):
        obj_iter = HandleJson(json_file).read_json()
        yield (antonomy, obj_iter)

  def load_specie_genome(self, obj_iter:Iterable):
    '''
    work on model Specie
    1. load anatomy from assembly_summary
    2. delete all records before insertion
    3. load data into Specie and Genome
    '''
    # truncate all data
    Specie.objects.all().delete()
    Genome.objects.all().delete()

    meta_names = ['genome_size', 'genome_size_ungapped', 'gc_percent',\
      'total_gene_count', 'protein_coding_gene_count', 'non_coding_gene_count']
    species, genomes = [], []
    for _, summary in obj_iter:
      # load Specie
      specie = Specie.objects.load_specie(summary)
      if specie:
        species.append(specie)
      # load genome
      specie_name = summary['organism_name'].replace(' ', '_')
      data = {
        'version': summary['assembly_accession'],
        'ftp_path': summary['ftp_path'],
        'specie': specie_name,
        'data_source': 'NCBI',
        'metadata': dict([(n, summary[n]) for n in meta_names]),
      }
      genome = Genome.objects.load_genome(data)
      if genome:
        genomes.append(genome)
    
    meta = {
      'specie_records': len(species),
      'genome_records': len(genomes),
    }
    return meta
  
  def download_genome(self):
    '''
    1- Download genome DNA sequences and annotations
    2. update Genome
    Note: It is supposed that genome metadata has been in Genome.
    '''
    # get genome object
    genome = Genome.objects.get_genome(self.data_source, self.specie, self.version)
    if not genome:
      return None

    # download data
    local_files = None
    if self.data_source == "NCBI":
      conn = NCBI(self.ref_dir, self.overwrite)
      local_path, local_files = conn.download_genome(
        genome.ftp_path, self.specie, self.version
      )

    # update database
    if local_files:
      # update db.Genome 
      genome.is_ready = True
      genome.local_path = local_path
      genome.save()
      # update db.Annotation
      Annotation.objects.load_annotations(genome, local_files)
      return {
        'local_files': local_files,
      }
    return None

  def molecular_annotation(self) -> list:
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
      if self.overwrite or (not meta):
        meta = wrapper.ncbi_fa_gff()
        wrapper.save_output(meta, True)
      # load meta into MolecularAnnotation
      # TODO: debuggging wrapper.load_output(): None is included sometimes
      meta = [i for i in wrapper.load_output() if i]
      res = MolecularAnnotation.objects.load(meta)
      return res
      
      







