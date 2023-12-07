TOOL_EXE = {
    'bowtie': [
        {'name': 'bowtie2',},
        {'name': 'bowtie2-build',},
        {'name': 'bowtie2-inspect',},
        {'name': 'bowtie',},
        {'name': 'bowtie-build',},
        {'name': 'bowtie-inspect',},
    ],
    'fastqc': [
        {'name': 'fastqc',},
    ],
    'hisat2': [
        {'name': 'hisat2',},
        {'name': 'hisat2-build',},
        {'name': 'hisat2-inspect'},
    ],
    'minimap2': [
        {'name': 'minimap2',},
    ],
    'samtools': [
        {'name': 'samtools',},
    ],
    'stringtie': [
        {'name': 'stringtie',},
    ],
    'tophat': [
        {'name': 'tophat',},
    ],
    'star': [
        {
            'name': 'STAR',
            'params': {
                '--sjdbOverhang': '99',
                '--runTreadN': '6',
            },
        },
    ],
}