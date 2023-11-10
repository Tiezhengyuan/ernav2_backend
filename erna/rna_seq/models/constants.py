METHODS = [
  {
    "method_name": "align_transcriptome",
    "description": "RNA-seq",
    "child_method": ["assemble_transcripts",],
    "tool_name": ["tophat2", "hisat2", "minimap2",],
  },
  {
    "method_name": "align_genome",
    "description": "",
    "child_method": [],
    "tool_name": ["bowtie2", "minimap2",],
  },
  {
    "method_name": "align_short_reads",
    "description": "",
    "child_method": ["count_reads",],
    "tool_name": ["bowtie2",],
  },
  {
    "method_name": "align_long_reads",
    "description": "",
    "child_method": ["count_reads",],
    "tool_name": ["minimap2",],
  },
  {
    "method_name": "build_index",
    "description": "build index for sequencing alignment",
    "child_method": ["align_transcriptome", "align_genome", "align_short_reads",],
    "tool_name": ["bowtie2-build", "hisat2-build",],
  },
  {
    "method_name": "assemble_transcripts",
    "description": "rna-seq",
    "child_method": ["count_reads",],
    "tool_name": ["stringtie"],
  },
  {
    "method_name": "count_reads",
    "description": "count reads for differential expression",
    "child_method": [],
    "tool_name": [],
  },
  {
    "method_name": "trim_sequences",
    "description": "miRNA-seq",
    "child_method": ["align_short_reads",],
    "tool_name": [],
  },
  {
    "method_name": "quality_control",
    "description": "sequencing quality control",
    "child_method": [],
    "tool_name": ["fastqc"],
  },
  {
    "method_name": "convert_format",
    "description": "sam-bam",
    "child_method": [],
    "tool_name": ["samtools"],
  },
]