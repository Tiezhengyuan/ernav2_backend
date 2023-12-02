'''
https://rnacentral.org
'''
import os
from Bio import SeqIO
from django.conf import settings
from utils.dir import Dir
from .connect_ftp import ConnectFTP


class ConnectRNACentral:
    url = "ftp.ebi.ac.uk"
    endpoint = '/pub/databases/RNAcentral/current_release/sequences/by-database/'

    def __init__(self):
        ref_dir = getattr(settings, 'REFERENCES_DIR')
        self.dir_local = os.path.join(ref_dir, "RNACentral")
        Dir(self.dir_local).init_dir()
    
    def download_data(self, file_name, overwrite:bool):
        local_file = ConnectFTP(self.url).download_file(
            endpoint=self.endpoint,
            file_name = file_name,
            local_path=self.dir_local,
            run_gunzip=False,
            overwrite=overwrite
        )
        return self.dir_local, local_file