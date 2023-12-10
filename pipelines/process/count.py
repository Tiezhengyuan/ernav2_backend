'''
retrieve data for further analysis
'''
from copy import deepcopy
import os
import pandas as pd
import pysam
from typing import Iterable
from .process import Process

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_seq.models import SampleProject
from pipelines.utils.utils import Utils

class Count:
    def __init__(self, params:dict):
        self.params = params

    def merge_read_counts(self):
        '''
        method: merge_read_counts
        '''
        rc_files = self.scan_rc_files()
        # merge RC.txt if they exist
        if rc_files.get('rc'):
            self.merge_rc_files(rc_files['rc'])
        if rc_files.get('stringtie'):
            self.stringtie_merge(rc_files['stringtie'], 'TPM')
            self.stringtie_merge(rc_files['stringtie'], 'FPKM')
        return None
    
    def scan_rc_files(self) -> Iterable:
        rc_files = {'rc': [], 'stringtie': []}
        for parent in self.params['parents']:
            outputs = parent.task_execution.get_output()
            for output in outputs:
                if 'abundance_file' in output:
                    item = (output['sample_name'], output['abundance_file'])
                    rc_files['stringtie'].append(item)
                elif 'RC' in output:
                    rc_files['rc'].append(output['RC'])
        return rc_files
                    
   
    def merge_rc_files(self, rc_files:list):
        '''
        merge multiple RC files into RC.txt
        '''
        df = pd.read_csv(rc_files[0], sep='\t', header=0)
        if len(rc_files) > 1:
            for rc_file in rc_files[1:]:
                tmp = pd.read_csv(rc_file, sep='\t', header=0)
                df = pd.merge(df, tmp, how='outer').fillna(0)
        df = df.convert_dtypes()
        outfile = os.path.join(self.params['output_dir'], 'RC.txt')
        df.to_csv(outfile, sep='\t', index=False)
        meta = {
            'count': 'RC',
            'outfile': outfile,
        }
        self.params['output'].append(meta)
        return df


    def stringtie_merge(self, rc_files:list, rc_type:str):
        '''
        stringtie
        rc_type: 'TPM' or 'FPKM'
        '''
        outfile = os.path.join(self.params['output_dir'], f'{rc_type}.txt')
        meta = {
            'count': rc_type,
            'outfile': outfile,
            'samples': [i[0] for i in rc_files],
        }

        sample_name, abund_file = rc_files.pop(0)
        df = self.read_abund(sample_name, abund_file, rc_type)
        for sample_name, abund_file in rc_files:
            tmp = self.read_abund(sample_name, abund_file, rc_type)
            df = pd.merge(df, tmp, how='outer').fillna(0)
        # 
        df = df.convert_dtypes()
        df.to_csv(outfile, index=False, sep='\t')
        meta['total'] = df[meta['samples']].sum().to_dict()
        self.params['output'].append(meta)
        return df

    def read_abund(self, sample_name, abund_file, rc_type):
        '''
        stringtie
        '''
        col_names = ['Gene ID', 'Gene Name', 'Reference', 'Strand',	'Start', 'End', 'Coverage']
        col_names.append(rc_type)
        df = pd.read_csv(abund_file, sep='\t')
        df = df[col_names]
        df.columns = df.columns.str.replace(rc_type, sample_name)
        return df


    def count_reads(self):
        '''
        reads counting from *.sam
        '''
        for parent in self.params['parents']:
            output = parent.task_execution.get_output()
            for item in [i for i in output if 'sam_file' in i]:
                sample_name = item['sample_name']
                output_prefix = os.path.join(self.params['output_dir'], sample_name)
                unaligned_file = output_prefix + ".unaligned.fa"
                rc = self.analyze_samfile(item['sam_file'], unaligned_file)
                # to txt
                df = pd.DataFrame.from_dict(rc, orient='index', columns=[sample_name,])
                outfile = output_prefix + ".RC.txt"
                df.to_csv(outfile, sep='\t', index_label='reference')
                self.params['output'].append({
                    'sample_name': sample_name,
                    'RC': outfile,
                    'unaligned': unaligned_file,
                })
  
    def analyze_samfile(self, sam_file, unaligned_file):
        '''
        count reads
        collect unaligned files
        '''
        rc = {}
        with open(unaligned_file, 'w') as f:
            infile = pysam.AlignmentFile(sam_file, 'r')
            for rec in infile.fetch():
                if rec.reference_name:
                    if rec.reference_name not in rc:
                        rc[rec.reference_name] = 1
                    else:
                        rc[rec.reference_name] += 1
                else:
                    if rec.seq:
                        record = SeqRecord(
                            Seq(rec.seq),
                            id=rec.qname,
                            description='',
                        )
                        SeqIO.write(record, f, 'fasta')
        return rc



  