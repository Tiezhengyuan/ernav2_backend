'''
'''
from Bio import SeqIO
from bioomics import Mirbase, RNACentral
from biosequtils import Dir
from django.conf import settings
import json
import os

from rna_seq.models import RNA


class ProcessNCRNA:

    def __init__(self, overwrite:bool=None):
        self.ref_dir = getattr(settings, 'REFERENCES_DIR')
        self.overwrite = overwrite
        # records indexed by specie
        self.data = {}
        self.annot = {}
        self.data_source = None
        self.annot_type = None

    def load_mirna(self, rna_type:str):
        '''
        download data from miRBase
        load data into eRNA
        '''
        self.data_source = 'miRBase'
        self.annot_type = rna_type

        # download sequence
        output_dir, fa_file = None, None
        conn = Mirbase(self.ref_dir, self.overwrite)
        if self.annot_type == 'miRNA_hairpin':
            output_dir, fa_file = conn.download_hairpin()
        elif self.annot_type == 'miRNA_mature':
            output_dir, fa_file = conn.download_mature()
        # split data by specie and Update models
        if fa_file:
            self.split_fasta(fa_file, self.desc_mirbase)
            metadata = self.save_split(output_dir)
            self.load_db(metadata)


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
        
    def load_lncrna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.data_source = 'RNACentral'
        self.annot_type = 'lncRNA'

        # download data: lncbook
        conn = RNACentral(self.ref_dir, self.overwrite)
        output_dir, local_file = conn.download_sequence('lncbook.fasta')
        # split local fasta
        self.split_fasta(local_file, self.desc_lncrna)
        metadata = self.save_split(output_dir)
        self.load_db(metadata)
    
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
    
    def load_piwirna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.data_source = 'RNACentral'
        self.annot_type = 'piRNA'

        # download data: pirbase
        conn = RNACentral(self.ref_dir, self.overwrite)
        output_dir, local_file = conn.download_sequence('pirbase.fasta')
        # split local fasta
        self.split_fasta(local_file, self.desc_pirbase)
        metadata = self.save_split(output_dir)
        # update model RNA
        self.load_db(metadata)

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

    def load_rrna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.data_source = 'RNACentral'
        self.annot_type = '5SrRNA'

        conn = RNACentral(self.ref_dir, self.overwrite)
        output_dir, local_file = conn.download_sequence('5srrnadb.fasta')
        metadata = self.format_fasta(local_file, output_dir)
        # update model RNA
        self.load_db([metadata,])

    def load_trna(self):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        self.data_source = 'RNACentral'
        self.annot_type = 'tRNA'

        # download
        conn = RNACentral(self.ref_dir, self.overwrite)
        output_dir, local_file = conn.download_sequence('gtrnadb.fasta')
        metadata = self.format_fasta(local_file, output_dir)
        # update model RNA
        self.load_db([metadata,])

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
        output_dir = os.path.join(output_dir, self.annot_type)
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

    def format_fasta(self, local_fa:str, output_dir:str):
        '''
        format fasta without split data
        '''
        output_dir = os.path.join(output_dir, self.annot_type)
        Dir(output_dir).init_dir()

        res = {}
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

    def load_db(self, metadata:dict):
        '''
        refresh model RNA
        '''
        RNA.objects.filter(data_source=self.data_source, \
            annot_type=self.annot_type).delete()
        RNA.objects.load_data(metadata)
    
    