'''
initialize models
example:
    python3 erna/manage.py shell < erna/scripts/init_reference.py
'''
# import sys
# print(sys.path)
from rna_seq.models import Genome, AlignerIndex

from pipelines.process.process_genome import ProcessGenome
from pipelines.process.process_ncrna import ProcessNCRNA


print('\n\n###Begin to refresh/update database###\n\n')

start_end = '5'.split('-')
start, end = int(start_end[0]), int(start_end[-1])
pool = range(start, end + 1)
for enter in pool:
    match enter:
        case 1:
            print('refresh Specie and Genome...')
            species = ProcessGenome('NCBI').retrieve_assembly_summary()
        case 2:
            print("Download human genome from NCBI...")
            ProcessGenome('NCBI', 'Homo_sapiens', 'GCF_000001405.40').download_genome()
        case 3:
            print("Download human genome from NCBI...")
            ProcessGenome('NCBI', 'Homo_sapiens', 'GCF_009914755.1').download_genome()
        case 4:
            print('refresh Genome...')
            genomes = Genome.objects.refresh()
        case 5:
            print('refresh index...')
            AlignerIndex.objects.refresh()
        case 6:
            print('load miRNA...')
            ProcessNCRNA().load_mirna('hairpin')
            ProcessNCRNA().load_mirna('mature')
        case 7:
            print('Process piwiRNA...')
            ProcessNCRNA().load_piwirna()
        case 8:
            print('Process long non-coding RNA...')
            ProcessNCRNA().load_lncrna()
        case 9:
            print('Process long rRNA...')
            ProcessNCRNA().load_rrna()
        case 10:
            print('Process long tRNA...')
            ProcessNCRNA().load_trna()


