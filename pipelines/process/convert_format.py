import os
from .process import Process

class ConvertFormat:

    def __init__(self, params:dict) -> None:
        self.params = params
    
    def sam_to_bam(self):
        '''
        convert SAM to BAM
        '''
        for parent_output in self.params['parent_outputs']:
            sample_name = parent_output.get('sample_name', '_')
            output_prefix = os.path.join(self.params['output_dir'], sample_name)

            # to BAM
            bamfile = output_prefix + '.bam'
            self.params['cmd'] = [
                self.params['tool'].exe_path,
                'view', '-b', parent_output['sam_file'],
                '>', bamfile
            ]
            self.params['force_run'] = False if os.path.isfile(bamfile) else True
            Process.run_subprocess(self.params)

            # sort BAM
            sorted_bamfile = output_prefix + '.bam.sorted'
            self.params['cmd'] = [
                self.params['tool'].exe_path,
                'sort',
                '-o', sorted_bamfile,
                bamfile,
            ]
            self.params['force_run'] = False if os.path.isfile(sorted_bamfile) else True
            Process.run_subprocess(self.params)
            
            self.params['output'].append({
                'cmd': ' '.join(self.params['cmd']),
                'sample_name': sample_name,
                'bam_file': bamfile,
                'sorted_bam_file': sorted_bamfile,
            })