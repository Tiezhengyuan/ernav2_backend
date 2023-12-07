'''
'''
from copy import deepcopy
from Bio import SeqIO
import os

from rna_seq.models import RNA
from pipelines.connector.connect_mirbase import ConnectMirbase
from pipelines.connector.connect_rnacentral import ConnectRNACentral
from pipelines.utils.dir import Dir


class ProcessNCRNA:

    def __init__(self, overwrite:bool=None):
        self.overwrite = True if overwrite else False
        # records indexed by specie
        self.data = {}

    def load_mirna(self):
        '''
        download data from miRBase
        load data into eRNA
        '''
        # download data
        local_path, local_files = ConnectMirbase().download_data(self.overwrite)
        # split data by specie and Update
        for rna_type, fa_path in local_files:
            self.split_fasta(fa_path, self.desc_mirbase)
            # export to fasta
            output_dir = os.path.join(local_path, rna_type)
            Dir(output_dir).init_dir()
            meta = {
                'rna_type': rna_type,
                'database': 'miRBase'
            }
            metadata = self.save_split(output_dir, meta)
            # update db.RNA
            RNA.objects.filter(database='miRBase', rna_type=rna_type).delete()
            RNA.objects.load_data(metadata)

    def load_lncrna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        conn = ConnectRNACentral()
        output_dir = os.path.join(conn.dir_local, 'lncrna')
        Dir(output_dir).init_dir()
        # download data: lncbook
        local_file = conn.download_data('lncbook.fasta', self.overwrite)
        self.split_fasta(local_file, self.desc_rna_central)
        
        # split local fasta
        meta = {
            'rna_type': 'lncRNA',
            'database': 'RNACentral'
        }
        metadata = self.save_split(output_dir, meta)
        # update db.RNA
        RNA.objects.filter(database='RNACentral', rna_type='lncRNA').delete()
        RNA.objects.load_data(metadata)

    def load_piwirna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        conn = ConnectRNACentral()
        output_dir = os.path.join(conn.dir_local, 'pirbase')
        Dir(output_dir).init_dir()

        # download data: pirbase
        local_file = conn.download_data('pirbase.fasta', self.overwrite)
        self.split_fasta(local_file, self.desc_pirbase)

        # split local fasta
        meta = {
            'rna_type': 'piwiRNA',
            'database': 'RNACentral'
        }
        metadata = self.save_split(output_dir, meta)
        # update db.RNA
        RNA.objects.filter(**meta).delete()
        RNA.objects.load_data(metadata)

    def load_rrna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        conn = ConnectRNACentral()
        local_file = conn.download_data('5srrnadb.fasta', self.overwrite)
        meta = {
            'rna_type': '5SrRNA',
            'database': 'RNACentral',
        }
        RNA.objects.filter(**meta).delete()
        meta['fa_path'] = local_file
        RNA.objects.load_data([meta,])

    def load_trna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        conn = ConnectRNACentral()
        local_file = conn.download_data('gtrnadb.fasta', self.overwrite)
        meta = {
            'rna_type': 'tRNA',
            'database': 'RNACentral',
        }
        RNA.objects.filter(**meta).delete()
        meta['fa_path'] = local_file
        RNA.objects.load_data([meta,])


    def split_fasta(self, local_file, func):
        '''
        split by specie
        update self.data
        '''
        for record in SeqIO.parse(local_file, 'fasta'):
            record.seq = record.seq.back_transcribe()
            # update record if necessary
            key = func(record)
            if key not in self.data:
                self.data[key] = []
            self.data[key].append(record)

    def save_split(self, output_dir, meta):
        '''
        '''
        metadata = []
        for key, records in self.data.items():
            _meta = deepcopy(meta)
            specie_name, organism_name, other_names, abb = key
            outfile = os.path.join(output_dir, f"{specie_name}.fa")
            # export fasta
            with open(outfile, 'w') as f:
                writer = SeqIO.FastaIO.FastaWriter(f, wrap=0)
                writer.write_file(records)
            # update metadata
            _meta.update({
                'fa_path': outfile,
                'specie_name': specie_name,
                'organism_name': organism_name,
                'other_names': other_names,
                'abb': abb,
            })
            metadata.append(_meta)
        return metadata

    def desc_pirbase(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[1]} {desc[2]}"
        specie_name = f"{desc[1]}_{desc[2]}"
        other_names = desc[3].replace(')', '').replace('(', '')
        abb = f"{desc[1][0].lower()}{desc[2][:2].lower()}"
        names = (specie_name, organism_name, other_names, abb)
        # update record
        desc = record.description.split(' ')
        record.id += f"|{desc[-1]}"
        record.description = ''
        return names

    def desc_rna_central(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[1]} {desc[2]}"
        specie_name = f"{desc[1]}_{desc[2]}"
        other_names = desc[3].replace(')', '').replace('(', '')
        abb = f"{desc[1][0].lower()}{desc[2][:2].lower()}"
        names = (specie_name, organism_name, other_names, abb)
        # update record
        record.id = record.description.replace(' ', '_')
        record.description = ''
        return names
    
    def desc_mirbase(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[2]} {desc[3]}"
        specie_name = f"{desc[2]}_{desc[3]}"
        abb = f"{desc[2][0].lower()}{desc[3][:2].lower()}"
        names = (specie_name, organism_name, None, abb)
        # update record
        desc = record.description.replace(' ', '_')
        record.id += f"_{desc}"
        record.description = ''
        return names
    