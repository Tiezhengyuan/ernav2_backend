'''
prepare commands used by external tools
'''

class ProcessCMD:
    
    @staticmethod
    def test():
        pass

    @staticmethod
    def aligner_build_index(tool, input_data:dict):
        '''
        aligner: hisat2, bowtie2 etc.
        '''
        cmd = [
            tool.exe_path,
            input_data['fa_path'],
            input_data['index_path'],
        ]
        output = {
            'cmd': ' '.join(cmd),
            'index_path': input_data['index_path'],
        }
        return cmd, output

   
    @staticmethod
    def hisat2_align(tool, input_data:dict):
        cmd = [
            tool.exe_path,
            '-x', input_data['index_path'],
        ]
        if input_data.get('bam'):
            cmd += ['-b', input_data['bam']]
        else:
            if input_data.get('R1'):
                cmd += ['-1', ','.join(input_data['R1'])]
            if input_data.get('R2'):
                cmd += ['-2', ','.join(input_data['R2'])]

        sam_file = input_data['output_prefix'] + '.sam'
        cmd += ['-S', sam_file]
        
        output_data = {
            'sample_name': input_data['sample_name'],
            'cmd': ' '.join(cmd),
            'sam_file': sam_file,
        }
        return cmd, output_data

    @staticmethod
    def bowtie2_align(tool, input_data:dict):
        cmd = [
            tool.exe_path,
            '-x', input_data['index_path'],
        ]
        if input_data.get('R1') and input_data.get('R2'):
            cmd += [
                '-1', ','.join(input_data['R1']),
                '-2', ','.join(input_data['R2']),
            ]
        elif input_data.get('bam'):
            cmd += ['-b', input_data['bam']]
        elif input_data.get('unaligned'):
            cmd += ['-f', input_data['unaligned']]
        else:
            raw_data = input_data.get('R1', []) + input_data.get('R2', [])
            cmd += ['-U', ','.join(raw_data)]

        sam_file = input_data['output_prefix'] + '.sam'
        # suppress @SQ header
        cmd += ['-S', sam_file]
        
        output_data = {
            'sample_name': input_data['sample_name'],
            'cmd': ' '.join(cmd),
            'sam_file': sam_file,
        }
        return cmd, output_data

    @staticmethod
    def stringtie_assemble(tool, input_data:dict):
        gtf_file = input_data['output_prefix'] + '.gtf'
        output = {
            'sample_name': input_data['sample_name'],
            'gtf_file': gtf_file,
        }
        cmd = [
            tool.exe_path,
            input_data['sorted_bam_file'],
            '-o', input_data['gtf_file'],
        ]
        if 'annot_file' in input_data:
            ballgown_file = input_data['output_prefix'] + '.ctab'
            abundance_file = input_data['output_prefix'] + '.abund'
            covered_file = input_data['output_prefix'] + '.cov'
            cmd += [
                '-G', input_data['annot_file'],
                '-e', '-b', ballgown_file,
                '-A', abundance_file,
                '-C', covered_file,
            ]            
            output.update({
                'annotation_file': input_data['annotation_file'],
                'ballgown_file': ballgown_file,
                'abundance_file': abundance_file,
                'covered_file': covered_file,
            })
        output['cmd'] = ' '.join(cmd)
        return cmd, output

    @staticmethod
    def star_build_index(tool, input_data:dict):
        cmd = [
          tool.exe_path,
          '--runThreadN', '6',
          '--runMode', 'genomeGenerate',
          '--genomeDir', input_data['index_path'],
          '--genomeFastaFiles', input_data['fa_path'],
          '--sjdbGTFfile', input_data['gtf_path'],
          '--sjdbOverhang', '99',
        ]
        output = {
            'cmd': ' '.join(cmd),
            'index_path': input_data['index_path'],
        }
        return cmd, output
    
    @staticmethod
    def star_align(tool, input_data:dict):
        cmd = [
            tool.exe_path,
            '--runThreadN', '6',
            '--genomeDir', input_data['index_path'],
            '--outFileNamePrefix', input_data['output_prefix'],
        ]
        if input_data.get('R1') or input_data.get('R2'):
            fq_files = input_data.get('R1', []) + input_data.get('R2', [])
            cmd += ['--readFilesIn', ','.join(fq_files),]
        output_data = {
            'sample_name': input_data['sample_name'],
            'cmd': ' '.join(cmd),
            'output_prefix': input_data['output_prefix'],
        }
        return cmd, output_data