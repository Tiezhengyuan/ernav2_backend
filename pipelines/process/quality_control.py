'''
quality control
'''
import os
import scanpy

from .process import Process


class QualityControl:
    def __init__(self, params:dict) -> None:
        self.params = params
    
    def __call__(self):
        if self.params['tool']:
            if self.params['tool'].tool_name == 'fastqc':
                for item in self.params['parent_outputs']:
                    fastq_files = item.get('R1', []) + item.get('R2', [])
                    for fq_file in fastq_files:
                        self.fastqc(item['sample_name'], fq_file)
                        Process.run_subprocess(self.params)

    def fastqc(self, sample_name, fq_file):
        cmd = [
            self.params['tool'].exe_path,
            fq_file,
            '-o', self.params['output_dir'],
            '-d', self.params['output_dir'],
        ]
        self.params['cmd'] = cmd
        fq_prefix = os.path.splitext(os.path.basename(fq_file))[0]
        self.params['output'].append({
            'sample_name': sample_name,
            'html': os.path.join(self.params['output_dir'], f"{fq_prefix}.html"),
            'zip': os.path.join(self.params['output_dir'], f"{fq_prefix}.html"),
        })
        return cmd