import subprocess


class Process:

    @staticmethod
    def run_subprocess(params:dict):
        '''
        Note: don't update object "params"
        '''
        if params.get('cmd') and params['force_run']:
            cmd = ' '.join(params['cmd'])
            # print(cmd)
            res = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if res.stdout:
                with open(f"{params.get('output_prefix', '_')}.out.log", 'w') as f:
                    f.writelines(res.stdout)
                print(res.stdout)
            if res.stderr:
                with open(f"{params.get('output_prefix', '_')}.err.log", 'w') as f:
                    f.writelines(res.stderr)
                print(res.stderr)
            return res
        return None