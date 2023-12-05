'''
initialize models
example:
    python3 erna/manage.py shell < erna/scripts/init_reference.py
'''
# import sys
# print(sys.path)
from rna_seq.models import Genome, Reference

from pipelines.process.process_genome import ProcessGenome
from pipelines.process.process_ncrna import ProcessNCRNA

enter = 11
run = False
print('\n\n###Begin to refresh/update database###\n\n')


if enter == 1 or run:
    print('refresh Specie and Genome...')
    species = ProcessGenome('NCBI').retrieve_assembly_summary()
    run=True

if enter == 2 or run:
    print('refresh Genome...')
    genomes = Genome.objects.refresh()
    run=True

if enter == 3 or run:
    print("Download default genome...")
    ProcessGenome('NCBI', 'Homo_sapiens', 'GCF_000001405.40').download_genome()
    ProcessGenome('NCBI', 'Homo_sapiens', 'GCF_009914755.1').download_genome()
    run=True

if enter == 4 or run:
    print('refresh Reference...')
    Reference.objects.refresh()
    run=True

if enter == 5 or run:
    # print('load miRNA...')
    ProcessNCRNA().load_mirbase(False)
    print('Process long non-coding RNA...')
    ProcessNCRNA().load_lncrnadb(False)
    print('Process piwiRNA...')
    ProcessNCRNA().load_pirbase(False)
    run=True
