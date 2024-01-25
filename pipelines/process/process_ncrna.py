'''
'''
import json
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
        self.annot = {}
        self.data_source = None
        self.annot_type = None

    def load_mirna(self, name:str):
        '''
        download data from miRBase
        load data into eRNA
        '''
        self.data_source = 'miRBase'
        # download data
        rna_type, local_path, fa_path = ConnectMirbase().download_data(name, self.overwrite)
        # split data by specie and Update
        if fa_path:
            self.annot_type = rna_type
            self.split_fasta(fa_path, self.desc_mirbase)
            # export to fasta
            output_dir = os.path.join(local_path, rna_type)
            metadata = self.save_split(output_dir)
            # update db.RNA
            RNA.objects.filter(data_source=self.data_source, annot_type=rna_type).delete()
            RNA.objects.load_data(metadata)

    def load_lncrna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.annot_type = 'lncRNA'
        self.data_source = 'RNACentral'

        conn = ConnectRNACentral()
        output_dir = os.path.join(conn.dir_local, 'lncrna')
        Dir(output_dir).init_dir()
        # download data: lncbook
        local_file = conn.download_data('lncbook.fasta', self.overwrite)
        self.split_fasta(local_file, self.desc_lncrna)
        
        # split local fasta
        metadata = self.save_split(output_dir)
        # update db.RNA
        RNA.objects.filter(data_source=self.data_source, annot_type=self.annot_type).delete()
        RNA.objects.load_data(metadata)

    def load_piwirna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.annot_type = 'piRNA'
        self.data_source = 'RNACentral'

        conn = ConnectRNACentral()

        # download data: pirbase
        output_dir = os.path.join(conn.dir_local, 'pirna')
        Dir(output_dir).init_dir()
        local_file = conn.download_data('pirbase.fasta', self.overwrite)
        self.split_fasta(local_file, self.desc_pirbase)

        # split local fasta
        metadata = self.save_split(output_dir)
        # update db.RNA
        RNA.objects.filter(data_source=self.data_source, annot_type=self.annot_type).delete()
        RNA.objects.load_data(metadata)


    def split_fasta(self, local_file, func):
        '''
        split by specie
        update self.data, self.annot
        '''
        for record in SeqIO.parse(local_file, 'fasta'):
            record.seq = record.seq.back_transcribe()
            # update record if necessary
            names, annot_rec = func(record)
            if names not in self.data:
                self.data[names] = []
                self.annot[names] = {}
            self.data[names].append(record)
            self.annot[names][annot_rec['ID']] = annot_rec
        
    def save_split(self, output_dir):
        '''
        '''
        Dir(output_dir).init_dir()
        metadata = []
        for key, records in self.data.items():
            specie_name, organism_name, other_names, abb = key
            # export sequences in fasta
            fa_path = os.path.join(output_dir, f"{specie_name}.fa")
            with open(fa_path, 'w') as f:
                writer = SeqIO.FastaIO.FastaWriter(f, wrap=0)
                writer.write_file(records)
            # export annotations in json
            annot_json = os.path.join(output_dir, f"{specie_name}.json")
            with open(annot_json, 'w') as f:
                json.dump(self.annot[key], f, indent=4)
            # update metadata
            meta = {
                'data_source': self.data_source,
                'annot_type': self.annot_type,
                'fa_path': fa_path,
                'annot_json': annot_json,
                'specie_name': specie_name,
                'organism_name': organism_name,
                'other_names': other_names,
                'abb': abb,
            }
            metadata.append(meta)
        return metadata


    def desc_mirbase(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[2]} {desc[3]}"
        specie_name = f"{desc[2]}_{desc[3]}"
        abb = f"{desc[2][0].lower()}{desc[3][:2].lower()}"
        names = (specie_name, organism_name, None, abb)
        # update record
        record.description = ''
        annot_rec = {
            'ID': record.id,
            'data_source': self.data_source,
            'annot_type': self.annot_type,
            'accession': desc[1],
            'organism_name': organism_name,
            'specie_name': specie_name,
        }
        return names, annot_rec
    
    def desc_pirbase(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[1]} {desc[2]}"
        specie_name = f"{desc[1]}_{desc[2]}"
        other_names = desc[3].replace(')', '').replace('(', '')
        abb = f"{desc[1][0].lower()}{desc[2][:2].lower()}"
        names = (specie_name, organism_name, other_names, abb)
        # update record
        record.description = ''
        annot_rec = {
            'ID': record.id,
            'data_source': self.data_source,
            'annot_type': self.annot_type,
            'accession': desc[-1],
            'organism_name': organism_name,
            'specie_name': specie_name,
            'other_names': other_names,
        }
        return names, annot_rec

    def desc_lncrna(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[1]} {desc[2]}"
        specie_name = f"{desc[1]}_{desc[2]}"
        other_names = desc[3].replace(')', '').replace('(', '')
        abb = f"{desc[1][0].lower()}{desc[2][:2].lower()}"
        names = (specie_name, organism_name, other_names, abb)
        # update record
        record.description = ''
        annot_rec = {
            'ID': record.id,
            'data_source': self.data_source,
            'annot_type': self.annot_type,
            'name': ' '.join(desc[4:]),
            'organism_name': organism_name,
            'specie_name': specie_name,
        }
        return names, annot_rec

    def load_rrna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.data_source = 'RNACentral'
        self.annot_type = '5SrRNA'

        conn = ConnectRNACentral()
        local_file = conn.download_data('5srrnadb.fasta', self.overwrite)
        output_dir = os.path.join(conn.dir_local, 'rrna')
        meta = self.format_fasta(local_file, output_dir)

        # updata db
        RNA.objects.filter(data_source=self.data_source, annot_type=self.annot_type).delete()
        RNA.objects.load_data([meta,])

    def format_fasta(self, local_fa:str, output_dir:str):
        res = {}
        Dir(output_dir).init_dir()
        file_name = os.path.basename(local_fa)
        # to fasta
        fa_file = os.path.join(output_dir, file_name)
        with open(fa_file, 'w') as f:
            for record in SeqIO.parse(local_fa, 'fasta'):
                annot_rec = self.desc_rna(record)                
                SeqIO.write(record, f, 'fasta-2line')
                res[annot_rec['ID']] = annot_rec
        # to json
        annot_json = os.path.join(output_dir, f"{os.path.splitext(file_name)[0]}.json")
        with open(annot_json, 'w') as f:
            json.dump(res, f, indent=4)
        # meta data
        meta = {
            'annot_type': self.annot_type,
            'data_source': self.data_source,
            'fa_path': fa_file,
            'annot_json': annot_json,
        }
        return meta

    def load_trna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.data_source = 'RNACentral'
        self.annot_type = 'tRNA'

        conn = ConnectRNACentral()
        local_file = conn.download_data('gtrnadb.fasta', self.overwrite)
        output_dir = os.path.join(conn.dir_local, 'trna')
        meta = self.format_fasta(local_file, output_dir)

        # updata db
        RNA.objects.filter(data_source=self.data_source, annot_type=self.annot_type).delete()
        RNA.objects.load_data([meta,])


    def desc_rna(self, record):
        desc = record.description.split(' ')
        organism_name = f"{desc[1]} {desc[2]}"
        specie_name = f"{desc[1]}_{desc[2]}"
        # update record
        record.description = ''
        annot_rec = {
            'ID': record.id,
            'data_source': self.data_source,
            'annot_type': self.annot_type,
            'name': ' '.join(desc[3:]),
            'organism_name': organism_name,
            'specie_name': specie_name,
        }
        return annot_rec
    
    