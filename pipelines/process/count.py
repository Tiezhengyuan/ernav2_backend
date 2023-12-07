'''
retrieve data for further analysis
'''
from copy import deepcopy
import os
import pandas as pd
import pysam
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
        run method: merge read counts
        '''
        if self.params['tool'].tool_name == 'stringtie':
            self.stringtie_merge_read_counts()
        return None

    def stringtie_merge_read_counts(self):
        for parent in self.params['parents']:
            parent_output = parent.task_execution.get_output()
            first = parent_output[0]
            df_tpm, df_fpkm = self.read_abund(first['abundance_file'], first['sample_name'])
        for item in parent_output[1:]:
            df1, df2 = self.read_abund(item['abundance_file'], item['sample_name'])
            df_tpm = pd.merge(df_tpm, df1, how='outer').fillna(0)
            df_fpkm = pd.merge(df_fpkm, df2, how='outer').fillna(0)
        #
        tpm_file = os.path.join(self.params['output_dir'], 'TPM.txt')
        df_tpm.to_csv(tpm_file, index=False, sep='\t')
        fpkm_file = os.path.join(self.params['output_dir'], 'FPKM.txt')
        df_fpkm.to_csv(fpkm_file, index=False, sep='\t')
        self.params['output'].append({
            'TPM': tpm_file,
            'FPKM': fpkm_file,
        })


    def read_abund(self, infile, sample_name):
        df=pd.read_csv(infile, sep='\t')
        # FPKM
        df1=df.loc[:, df.columns != 'TPM']
        df1.columns = df1.columns.str.replace('FPKM', sample_name)
        # TPM
        df2=df.loc[:, df.columns != 'FPKM']
        df2.columns = df2.columns.str.replace('TPM', sample_name)
        return df1, df2


    def count_reads(self):
        '''
        reads counting
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
        with open(unaligned_file, 'w') as outfile:
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
                        SeqIO.write(record, outfile, 'fasta')
        return rc


    def merge_read_counts(self):
        '''
        reads counting
        '''
        rc_files, meta = [], {}
        for parent in self.params['parents']:
            output = parent.task_execution.get_output()
            rc_files += [i['RC'] for i in output if 'RC' in i]
        # merge RC files
        if rc_files:
            df = self.merge_rc_files(rc_files)
            outfile = os.path.join(self.params['output_dir'], 'RC.txt')
            df.to_csv(outfile, sep='\t', index=False)
            meta['RC_file'] = outfile
            # TODO: normalization

        self.params['output'].append(meta)
    
    def merge_rc_files(self, rc_files):
        '''
        merge multiple RC files into RC.txt
        '''
        df = pd.read_csv(rc_files[0], sep='\t', header=0)
        if len(rc_files) > 1:
            for rc_file in rc_files[1:]:
                tmp = pd.read_csv(rc_file, sep='\t', header=0)
                df = pd.merge(df, tmp, how='outer').fillna(0)
        df = df.convert_dtypes()
        return df

  