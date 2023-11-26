'''
quality control
'''
import os
from .process import Process


class QualityControl:
    def __init__(self, params:dict) -> None:
        self.params = params
    
    def __call__(self):
        if self.params['tool']:
            if self.params['tool'].tool_name == 'fastqc':
                for sample_file in self.params['sample_files']:
                    self.fastqc(sample_file)
                    # run process
                    Process.run_subprocess(self.params)

    def fastqc(self, sample_file):
        raw = sample_file.raw_data
        cmd = [
            self.params['tool'].exe_path,
            os.path.join(raw.file_path, raw.file_name),
            '-o', self.params['output_dir'],
            '-d', self.params['output_dir'],
        ]
        print(cmd)
        self.params['cmd'] = cmd
        fq_prefix = os.path.splitext(raw.file_name)[0]
        self.params['output'].append({
            'sample_name': sample_file.sample.sample_name,
            'html': os.path.join(self.params['output_dir'], f"{fq_prefix}.html"),
            'zip': os.path.join(self.params['output_dir'], f"{fq_prefix}.html"),
        })
        return cmd