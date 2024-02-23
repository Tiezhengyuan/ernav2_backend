METHODS = [
  {
    "method_name": "import_data",
    "description": "bound with task 00",
  },
  {
    "method_name": "align_transcriptome",
    "description": "RNA-seq",
    "child_method": ["convert_format",],
    "exe_name": [
        "hisat2",
        # "tophat2",
        # "minimap2",
        # 'STAR',
    ],
  },
  # {
  #   "method_name": "align_genome",
  #   "description": "",
  #   "child_method": [],
  #   "exe_name": ["bowtie2", "minimap2",],
  # },
  {
    "method_name": "align_short_reads",
    "description": "",
    "child_method": ["count_reads", "convert_format",],
    "exe_name": ["bowtie2",],
  },
  # {
  #   "method_name": "align_long_reads",
  #   "description": "",
  #   "child_method": ["count_reads"],
  #   "exe_name": ["minimap2",],
  # },
  {
    "method_name": "build_genome_index",
    "description": "build index against genome DNA for sequencing alignment",
    "child_method": [
        "align_transcriptome",
        # "align_genome",
    ],
    "exe_name": [
        "hisat2-build",
        # "bowtie2-build",
        # 'STAR',
    ],
  },
  {
    "method_name": "build_index",
    "description": "build index for sequencing alignment",
    "child_method": [
        # "align_long_reads",
        "align_short_reads",
    ],
    "exe_name": [
        # "bowtie1-build",
        "bowtie2-build",
    ],
  },
  {
    "method_name": "assemble_transcripts",
    "description": "assemble transcripts after transcripts alignment",
    "child_method": ["merge_transcripts", "merge_read_counts"],
    "exe_name": ["stringtie"],
  },
  {
    "method_name": "count_reads",
    "description": "count reads and collect unaligned reads",
    "child_method": ["merge_read_counts",],
    "exe_name": [],
  },
  {
    "method_name": "merge_read_counts",
    "description": "count reads for differential expression",
    "child_method": [],
    "exe_name": [],
  },
  {
    "method_name": "merge_transcripts",
    "description": "merge transcripts into a non-redundant set of transcripts",
    "child_method": [],
    "exe_name": ["stringtie"],
  },
  {
    "method_name": "trim_sequences",
    "description": "miRNA-seq",
    "child_method": ["align_short_reads",],
    "exe_name": [],
    "default_params": {'adapter_3end': 'TGGAATTCTCGGGTGCCAAGG',},
  },
  {
    "method_name": "quality_control",
    "description": "sequencing quality control",
    "child_method": [],
    "exe_name": ["fastqc"],
  },
  {
    "method_name": "convert_format",
    "description": "sam-bam",
    "child_method": ['assemble_transcripts'],
    "exe_name": ["samtools"],
  },
]

