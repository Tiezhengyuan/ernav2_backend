'''
retrieve data for further analysis
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biofile import GFF
import os
import pandas as pd
import pysam
from rnaseqdata import NodeData, dump_seqdata
from typing import Iterable

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
            gene_json = self.params['genome_annot']['gff_json'].get('gene')
            if gene_json:
                gene = GFF(gene_json).lift_attribute('ID')
                self.params['seqdata'].put_variables(gene)
            self.stringtie_merge(rc_files['stringtie'], 'TPM')
            self.stringtie_merge(rc_files['stringtie'], 'FPKM')
        
        # update seqdata.obj
        dump_seqdata(self.params['seqdata'], self.params['seqdata_path'])
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
                    path = (output['sample_name'], output['RC'])
                    rc_files['rc'].append(path)
        return rc_files
                    
   
    def merge_rc_files(self, rc_files:list):
        '''
        merge multiple RC files into RC.txt
        '''
        # print(rc_files[-4:])
        # update SeqData
        rc_node = NodeData(self.params['seqdata'].root, 'RC')
        for sample_name, rc_file in rc_files[-4:]:
            rc = pd.read_csv(rc_file, sep='\t', index_col=0, header=0)
            rc.name = sample_name
            rc_node.put_data(rc.iloc[:,0])
        self.params['seqdata'].nodes['RC'] = rc_node

        # export
        df = self.params['seqdata'].to_df('RC', 1)
        outfile = os.path.join(self.params['output_dir'], 'RC.txt')
        df.to_csv(outfile, sep='\t', index=True, header=True)
        meta = {
            'count': 'RC',
            'RC': outfile,
            'shape': df.shape,
        }
        df = df.T
        self.params['output'].append(meta)
        outfile = os.path.join(self.params['output_dir'], 'RC_T.txt')
        df.to_csv(outfile, sep='\t', index=True, header=True)
        meta = {
            'count': 'RC',
            'RC': outfile,
            'shape': df.shape,
        }
        self.params['output'].append(meta)
        return rc_node


    def stringtie_merge(self, rc_files:list, rc_type:str):
        '''
        stringtie
        rc_type: 'TPM' or 'FPKM'
        '''
        outfile = os.path.join(self.params['output_dir'], f'{rc_type}.txt')
        meta = {
            'count': rc_type,
            rc_type: outfile,
            'samples': [i[0] for i in rc_files],
        }
        # update seqdata
        node = NodeData(self.params['seqdata'].root, rc_type)
        for sample_name, abund_file in rc_files:
            df = pd.read_csv(abund_file, sep='\t', index_col=0, header=0)
            node.put_data(pd.Series(df[rc_type], name=sample_name))
        # remove zeros
        node.X = node.X.loc[:,node.X.sum(axis=0)>0]
        self.params['seqdata'].nodes[rc_type] = node

        # sample in columns
        df = self.params['seqdata'].to_df(rc_type, 1).T
        df = df.convert_dtypes()
        df.to_csv(outfile, index=True, index_label='ID', header=True, sep='\t')

        meta['total'] = df[meta['samples']].sum().to_dict()
        self.params['output'].append(meta)
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

