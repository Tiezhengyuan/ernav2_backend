'''
initialize models
example:
    python3 erna/manage.py shell < erna/scripts/init_erna.py
'''
from rna_seq.models import *
from annot.models import *
from process.process_raw_data import ProcessRawData
from process.process_genome import ProcessGenome



# annotation-related tables

# refresh Specie and Genome
species = ProcessGenome('NCBI').retrieve_assembly_summary()

#refresh downloaded genome
genomes = Genome.objects.refresh()

# refresh reference
references = Reference.objects.refresh()


'''
# refresh Method
methods = Method.objects.refresh()

# refresh Tool
tools = Tool.objects.refresh()


# refresh MethodTool
method_tools = MethodTool.objects.refresh()

# refresh MethodRelation
method_relations = MethodRelation.objects.refresh()


# import default pipelines
Pipeline.objects.all().delete()
mrna_seq = [
    ("align_transcriptome", 'hisat2', None),
    ("assemble_transcripts", 'stringtie', None), 
    ("count_reads", None, None),
]
steps = Pipeline.objects.load_pipeline('mRNA-Seq', mrna_seq)


mirna_seq = [
    ("trim_sequences", None, None),
    ("align_short_reads", "bowtie2", None),
    ("count_reads", None, None),
]
steps = Pipeline.objects.load_pipeline('miRNA-Seq', mirna_seq)


# Delete all data in RawData, Sample, SampleFile, SampleProject
# refresh RawData
res = ProcessRawData().reset_sample()
'''