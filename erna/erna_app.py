'''
entrance into pipelines
'''
import os
import sys
from dotenv import load_dotenv
load_dotenv()
DIR_PROJECT = os.path.dirname(os.path.dirname(__file__))
sys.path.append(DIR_PROJECT)

import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'erna.settings')
django.setup()

def main(args):
  if len(args) < 1:
    print("method should be defined.")
    sys.exit(1)

  match args[0]:
    case 'scan_raw_data':
      '''
      Get batch_name and path from raw data 
      example:
      python3 erna_app.py scan_raw_data
      '''
      from pipelines.process.process_raw_data import ProcessRawData
      return ProcessRawData().scan_raw_data()

    case 'refresh_raw_data':
      '''
      scan raw data and load to db.RawData
      example:
      python3 erna_app.py scan_raw_data
      '''
      from pipelines.process.process_raw_data import ProcessRawData
      return ProcessRawData().refresh_raw_data()

    case 'assembly_summary':
      '''
      download assembly_summary.txt
      example:
      python3 erna_app.py assembly_summary NCBI
      '''
      if len(args) >= 2:
        from pipelines.process.process_genome import ProcessGenome
        return ProcessGenome(args[1]).retrieve_assembly_summary()

    case 'download_genome':
      '''
      download refseq of genome fro FTP
      example:
      python3 erna_app.py download_genome NCBI "Homo sapiens" GCF_000001405.40
      '''
      if len(args)>=4:
        from pipelines.process.process_genome import ProcessGenome
        data_source, specie, version = args[1:]
        p = ProcessGenome(data_source, specie, version)
        return p.download_genome()

    case 'build_index':
      '''
      build index for aligner namely bowtie
      example: 
      python3 erna_app.py build_index NCBI "Homo sapiens" GCF_000001405.40 bowtie2
      '''
      if len(args)>=5:
        from pipelines.process.align import Align
        data_source, specie, version, aligner = args[1:]
        p = Align(aligner)
        return p.build_index(data_source, specie, version)

    case 'genome_alignment':
      if len(args)>=5:
        from pipelines.process.align import Align
        data_source, specie, version, aligner = args[1:]
        p = Align(data_source, specie, version, aligner)
        return p.build_index()

    case 'trim_adapter':
      '''
      trim adapter sequences from reads for miRNA-seq
      example: 
        python3 erna_app.py trim_adapter ATGGCG \
          /home/yuan/bio/erna_v2/pipelines/tests/data/reads_1.fq \
          /home/yuan/bio/erna_v2/pipelines/tests/temp/reads_1_trimmed.fq
      '''
      if len(args) >=4:
        from pipelines.process.trim_adapter import TrimAdapter
        params =  dict([(k,v) for k, v in zip(['adapter', 'input',\
          'output',], args[1:])])
        return TrimAdapter(params)()






    
    case 'genome_assembly':
      if len(args)>=2:
        from pipelines.process import Assembly
        p = Assembly(args[1])
        return p.genome_aseembly()
    case 'count_reads':
      if len(args)>=3:
        from pipelines.process import Collect
        p = Collect(args[1])
        return p.count_reads()
  print('wrong arguments')

if __name__ == '__main__':
  args = sys.argv[1:]
  main(args)
