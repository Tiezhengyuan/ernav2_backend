"""
download data from NCIB FTP
"""
import os
from ftplib import FTP
from django.conf import settings

from rna_seq.models import Specie, Genome
from utils import Dir, HandleJson
from .connect_ftp import ConnectFTP
from .connect_ftp2 import ConnectFTP2

# lock the scope of groups in local version
# ANATOMY_GROUPS = ['archaea', 'bacteria', 'fungi', 'invertebrate', 'plant',
#     'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral', ]
ANATOMY_GROUPS = ['vertebrate_mammalian',]

class ConnectNCBI(ConnectFTP):
    url = 'ftp.ncbi.nlm.nih.gov'

    def __init__(self):
        super(ConnectNCBI, self).__init__(self.url)
        ref_dir = getattr(settings, 'REFERENCES_DIR')
        self.dir_local = os.path.join(ref_dir, "NCBI")

    def download_assembly_summary(self):
        '''
        download assembly_summary.txt
        '''
        res = {}
        for antonomy in ANATOMY_GROUPS:
            outdir = os.path.join(self.dir_local, 'assembly_summary', antonomy)
            local_file = self.download_file(
                endpoint = f'genomes/refseq/{antonomy}/',
                file_name = 'assembly_summary.txt',
                local_path = outdir
            )
            res[antonomy] = local_file
        return res

    def download_genome(self, specie:str, version:str, overwrite:bool=None):
        '''
        download genome including subdirectories and files
        '''
        ftp_path = Genome.objects.get_ftp_path(specie, version)

        # download annotations
        local_path = os.path.join(self.dir_local, 'genome', specie, version)
        local_files = self.download_files(
            endpoint = ftp_path.replace('https://ftp.ncbi.nlm.nih.gov/', ''),
            match = '.gz',
            local_path = local_path,
            overwrite=overwrite,
        )
        # download index for alignment
        return local_path, local_files


    def load_genomes(self):
        '''
        work on model Genome
        1. delete all
        2. load genome from assembly_summary
        '''
        # truncate all data
        Genome.objects.all().delete()

        # retrieve data from json
        n = {}
        meta_names = ['genome_size', 'genome_size_ungapped', 'gc_percent',\
            'total_gene_count', 'protein_coding_gene_count', 'non_coding_gene_count']
        for antonomy in ANATOMY_GROUPS:
            n[antonomy] = []
            json_file = os.path.join(self.dir_local, 'assembly_summary',\
                antonomy, 'assembly_summary.json')
            obj = HandleJson(json_file).read_json()
            for _, summary in obj:
                specie_name = summary['organism_name'].replace(' ', '_')
                data = {
                    'version': summary['assembly_accession'],
                    'ftp_path': summary['ftp_path'],
                    'specie': specie_name,
                    'data_source': 'NCBI',
                    'metadata': dict([(n, summary[n]) for n in meta_names]),
                }
                res = Genome.objects.load_genome(data)
                if res:
                    n[antonomy].append(data['specie'])
            n[antonomy] = len(n[antonomy])
        return n
    
    def load_species(self):
        '''
        work on model Specie
        1. delete all
        2. load anatomy from assembly_summary
        '''
        # truncate all data
        Specie.objects.all().delete()

        # retrieve data from json
        n = {}
        names = ['taxid', 'organism_name', 'group']
        for antonomy in ANATOMY_GROUPS:
            n[antonomy] = []
            json_file = os.path.join(self.dir_local, 'assembly_summary',\
                antonomy, 'assembly_summary.json')
            obj = HandleJson(json_file).read_json()
            for _, summary in obj:
                data = dict([(n, summary[n]) for n in names])
                data['organism_name'] = data['organism_name'].replace(' ', '_')
                if data['organism_name'] not in n[antonomy]:
                    Specie.objects.create(**data)
                    n[antonomy].append(data['organism_name'])
            n[antonomy] = len(n[antonomy])
        return n
    
    def download_gene_data(self):
        '''
        download /gene/DATA including subdirectories and files
        '''
        local_files = self.download_tree(
            local_name = os.path.join(self.dir_local, 'gene', 'DATA'),
            endpoint = 'gene/DATA',
            match = '.gz$'
        )
        return local_files

    def download_pubmed(self):
        '''
        download /PubMed including subdirectories and files
        '''
        ConnectFTP2.download_tree(
            ftp_url = self.url,
            endpoint = '/pubmed',
            match = '.gz',
            local_path = self.dir_local
        )