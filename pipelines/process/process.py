import subprocess


class Process:

    @staticmethod
    def run_subprocess(params:dict):
        '''
        Note: don't update object "params"
        '''
        if params.get('cmd'):
            res = subprocess.run(params['cmd'], capture_output=True, text=True)
            if res.stdout:
                print(res.stdout)
                with open(f"{params.get('output_prefix', '_')}.out.log", 'w') as f:
                    f.writelines(res.stdout)
            if res.stderr:
                print(res.stderr)
                with open(f"{params.get('output_prefix', '_')}.err.log", 'w') as f:
                    f.writelines(res.stderr)
            return res
        return None