'''
'''
from copy import deepcopy
from Bio import SeqIO
import os

from rna_seq.models import NonCodingRNA
from pipelines.connector.connect_mirbase import ConnectMirbase
from pipelines.connector.connect_rnacentral import ConnectRNACentral
from pipelines.utils.dir import Dir


class ProcessNCRNA:

    def __init__(self):
        pass

    def load_mirbase(self, overwrite:bool=None):
        '''
        download data from miRBase
        load data into eRNA
        '''
        # download data
        local_path, local_files = ConnectMirbase().download_data(overwrite)
        print(local_files)
        # split data by specie and Update
        for rna_type, fa_path in local_files:
            data = self.split_fasta(fa_path, self.desc_mirbase)
            # export to fasta
            output_dir = os.path.join(local_path, rna_type)
            Dir(output_dir).init_dir()
            meta = {
                'rna_type': rna_type,
                'database': 'miRBase'
            }
            metadata = self.save_split(data, output_dir, meta)
            # update db.NonCodingRNA
            NonCodingRNA.objects.filter(database='miRBase', rna_type=rna_type).delete()
            NonCodingRNA.objects.load_data(metadata)

    def load_lncrnadb(self, overwrite:bool=None):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        # download data
        c = ConnectRNACentral()
        local_path, local_file = c.download_data('lncrnadb.fasta', overwrite)
        data = self.split_fasta(local_file, self.desc_lncrnadb)
        output_dir = os.path.join(local_path, 'lncrna')
        Dir(output_dir).init_dir()
        meta = {
            'rna_type': 'lncRNA',
            'database': 'RNACentral'
        }
        metadata = self.save_split(data, output_dir, meta)
        # update db.NonCodingRNA
        NonCodingRNA.objects.filter(database='RNACentral', rna_type='lncRNA').delete()
        NonCodingRNA.objects.load_data(metadata)

    def load_pirbase(self, overwrite:bool=None):
        '''
        download data from https://rnacentral.org/
        load data into eRNA
        '''
        # download data
        c = ConnectRNACentral()
        local_path, local_file = c.download_data('pirbase.fasta', overwrite)
        data = self.split_fasta(local_file, self.desc_pirbase)
        output_dir = os.path.join(local_path, 'pirbase')
        Dir(output_dir).init_dir()
        meta = {
            'rna_type': 'piwiRNA',
            'database': 'RNACentral'
        }
        metadata = self.save_split(data, output_dir, meta)
        # update db.NonCodingRNA
        NonCodingRNA.objects.filter(database='RNACentral', rna_type='piwiRNA').delete()
        NonCodingRNA.objects.load_data(metadata)
    
    def split_fasta(self, local_file, func):
        '''
        split by specie
        '''
        res = {}
        for record in SeqIO.parse(local_file, 'fasta'):
            record.seq = record.seq.back_transcribe()
            # update record if necessary
            key = func(record)
            if key not in res:
                res[key] = []
            res[key].append(record)
        return res


    def save_split(self, data, output_dir, meta):
        '''
        '''
        metadata = []
        for key, records in data.items():
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

    def desc_lncrnadb(self, record):
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
    