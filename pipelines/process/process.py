import os
import subprocess


class Process:

    @staticmethod
    def run_subprocess(params:dict):
        '''
        Note: don't update object "params"
        '''
        if params.get('cmd') and params['force_run']:
            cmd = ' '.join(params['cmd'])
            res = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if res.stdout:
                log = f"{params.get('output_prefix', '_')}.out.log"
                with open(log, 'w') as f:
                    f.writelines(res.stdout)
                print(res.stdout)
            if res.stderr:
                err_log = f"{params.get('output_prefix', '_')}.err.log"
                with open(err_log, 'w') as f:
                    f.writelines(res.stderr)
                print(res.stderr)
            return res
        return None

    @staticmethod
    def uncompress_gz(gz_file:str):
        if not os.path.isfile(gz_file):
            return None
        # run gzip
        print(f"Try to uncompress {gz_file}")
        cmd = ' '.join(['gzip', '-d', gz_file])
        res = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if res.stdout:
            print(res.stdout)
        if res.stderr:
            print(res.stderr)
        return res