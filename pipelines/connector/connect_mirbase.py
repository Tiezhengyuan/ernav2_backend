import os
from Bio import SeqIO
from io import StringIO
import lxml.html as html
from urllib.request import urlopen

from django.conf import settings
from utils.dir import Dir


class ConnectMirbase:
    url = 'https://mirbase.org/download/CURRENT/'

    def __init__(self):
        ref_dir = getattr(settings, 'REFERENCES_DIR')
        self.dir_local = os.path.join(ref_dir, "miRBase")
        Dir(self.dir_local).init_dir()
    
    def download_mirbase(self, overwrite:bool=None):
        local_files = []
        for file_name in ['hairpin.fa', 'mature.fa']:
            local_file = os.path.join(self.dir_local, file_name)
            self.parse_fasta(f"{self.url}/{file_name}",local_file)
            local_files.append(local_file)
        return self.dir_local, local_files

    def parse_fasta(self, url_path, local_file):
        # parse xml
        parsed = html.parse(urlopen(url_path))
        doc = parsed.getroot()
        for br in doc.xpath("*//br"):
            br.tail = "\n" + br.tail if br.tail else "\n"
        # parse fasta
        fa_io = StringIO(doc.text_content())
        records = SeqIO.parse(fa_io, 'fasta')
        with open(local_file, 'w') as f:
            for rec in records:
                rec.seq = rec.seq.back_transcribe()
                # print(rec, str(rec.seq))
                SeqIO.write(rec, f, 'fasta')



