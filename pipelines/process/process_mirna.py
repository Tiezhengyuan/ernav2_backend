'''
'''
from Bio import SeqIO
import os

from rna_seq.models import NonCodingRNA
from connector.connect_mirbase import ConnectMirbase
from utils.dir import Dir


class ProcessMiRNA:

    def __init__(self):
        pass

    def load_mirbase(self, overwrite:bool=None):
        '''
        download data from miRBase
        load data into eRNA
        '''
        # download data
        local_path, local_files = ConnectMirbase().download_mirbase(overwrite)
        print(local_path, local_files)

        # split data by specie and Update
        for rna_type, fa_path in local_files:
            data = {}
            for rec in SeqIO.parse(fa_path, 'fasta'):
                rec.seq = rec.seq.back_transcribe()
                abb = rec.id.split('-', 1)[0]
                if abb not in data:
                    data[abb] = []
                data[abb].append(rec)
            # export to fasta
            metadata = []
            output_dir = os.path.join(local_path, rna_type)
            Dir(output_dir).init_dir()
            for abb, records in data.items():
                outfile = os.path.join(output_dir, f"{rna_type}_{abb}.fa")
                SeqIO.write(records, outfile, 'fasta')
                metadata.append({
                    'abb': abb,
                    'fa_path': outfile,
                    'rna_type': rna_type,
                    'db': 'miRBase',
                })
            # update db.NonCodingRNA
            NonCodingRNA.objects.load_data(metadata)

