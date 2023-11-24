import os
from .process import Process

class ConvertFormat:

    def __init__(self, params:dict) -> None:
        self.params = params
    
    def sam_to_bam(self):
        '''
        convert SAM to BAM
        '''
        for parent in self.params['parents']:
            parent_output = parent.task_execution.get_output()
            for item in parent_output:
                sample_name = item.get('sample_name', '_')
                output_prefix = os.path.join(self.params['output_dir'], sample_name)
                bamfile = output_prefix + '.sam'
                self.params['cmd'] = [
                    self.params['tool'].exe_path,
                    'view',
                    '-b',
                    item['sam_file'],
                    '>',
                    bamfile
                ]
                Process.run_subprocess(self.params)
                self.params['output'].append({
                    'cmd': ' '.join(self.params['cmd']),
                    'sample_name': sample_name,
                    'bam_file': bamfile,
                })