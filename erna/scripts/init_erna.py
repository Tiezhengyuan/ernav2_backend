'''
initialize models
example:
    python3 erna/manage.py shell < erna/scripts/init_erna.py
'''
import sys
from rna_seq.models import *
from process.process_raw_data import ProcessRawData
from process.process_genome import ProcessGenome
from process.process_mirna import ProcessMiRNA

enter = 1
run = False
print('\n\n###Begin to refresh/update database###\n\n')


if enter == 1 or run:
    print('Refresh Method...')
    methods = Method.objects.refresh()
    run=True

if enter == 2 or run:
    print('Refresh Tool...')
    tools = Tool.objects.refresh()
    run=True

if enter == 3 or run:
    print('refresh MethodTool...')
    method_tools = MethodTool.objects.refresh()
    run=True

if enter == 4 or run:
    print('refresh MethodRelation...')
    method_relations = MethodRelation.objects.refresh()
    run=True

if enter == 5 or run:
    print('import default pipelines of mRNA-Seq...')
    mrna_seq = [
        ("build_genome_index", 'hisat2-build', None),
        ("align_transcriptome", 'hisat2', None),
        ("convert_format", 'samtools', None),
        ("assemble_transcripts", 'stringtie', None), 
        ("merge_transcripts", 'stringtie', None), 
        ("count_reads", None, None),
    ]
    Pipeline.objects.filter(pipeline_name='mRNA-Seq').delete()
    steps = Pipeline.objects.load_pipeline('mRNA-Seq', mrna_seq)
    run=True

if enter == 6 or run:
    print('import default pipelines of miRNA-Seq...')
    mirna_seq = [
        ("build_index", 'bowtie2-build', None),
        ("trim_sequences", None, None),
        ("align_short_reads", "bowtie2", None),
        ("count_reads", None, None),
    ]
    Pipeline.objects.filter(pipeline_name='miRNA-Seq').delete()
    steps = Pipeline.objects.load_pipeline('miRNA-Seq', mirna_seq)
    run=True

if enter == 7 or run:
    print('refresh RawData...')
    # Delete all data in RawData, Sample, SampleFile, SampleProject
    res = ProcessRawData().reset_sample()
    run=True

if enter == 8 or run:
    print('refresh Specie and Genome...')
    species = ProcessGenome('NCBI').retrieve_assembly_summary()
    run=True

if enter == 9 or run:
    print('refresh Genome...')
    genomes = Genome.objects.refresh()
    run=True

if enter == 10 or run:
    print("Download default genome...")
    ProcessGenome('NCBI', 'Homo_sapiens', 'GCF_000001405.40').download_genome()
    ProcessGenome('NCBI', 'Homo_sapiens', 'GCF_009914755.1').download_genome()
    run=True

if enter == 11 or run:
    print('refresh Reference...')
    Reference.objects.refresh()
    run=True

if enter == 12 or run:
    print('Refresh NonCodingRNA...')
    ProcessMiRNA().load_mirbase(False)
    run=True

