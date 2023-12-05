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

from pipelines.process.process_raw_data import ProcessRawData
from pipelines.process.process_genome import ProcessGenome
from pipelines.process.process_ncrna import ProcessNCRNA
from pipelines.process.execute_tasks import ExecuteTasks

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
      return ProcessRawData().scan_raw_data()

    case 'refresh_raw_data':
      '''
      scan raw data and load to db.RawData
      example:
      python3 erna_app.py scan_raw_data
      '''
      return ProcessRawData().refresh_raw_data()

    case 'assembly_summary':
      '''
      download assembly_summary.txt
      example: python3 erna/erna_app.py assembly_summary NCBI
      '''
      if len(args) >= 2:
        return ProcessGenome(args[1]).retrieve_assembly_summary()

    case 'download_genome':
      '''
      download refseq of genome fro FTP
      example:
      python3 erna_app.py download_genome NCBI "Homo sapiens" GCF_000001405.40
      '''
      if len(args)>=4:
        data_source, specie, version = args[1:]
        p = ProcessGenome(data_source, specie, version)
        return p.download_genome()

    case 'load_mirbase':
      '''
      download data from miRBase and load them into eRNA
      example: python3 erna_app.py load_mirbase
      '''
      overwrite = eval(args[1]) if len(args)>1 and args[1] else False
      return ProcessNCRNA().load_mirbase(overwrite)
      
    case 'execute_tasks':
      '''
      example:
      python3 erna/erna_app.py execute_tasks P00001 T01
      '''
      if len(args)>=2:
        project_id = args[1]
        task_id = args[2] if len(args)==3 else None
        chain = True if len(args)==4 else False
        return ExecuteTasks(project_id, task_id, chain)()


    case 'build_genome_index':
      '''
      build index for aligner namely bowtie
      example: 
      python3 erna/erna_app.py build_index "Homo sapiens" GCF_000001405.40 bowtie 2.5.2
      python3 erna/erna_app.py build_index "Homo sapiens" GCF_000001405.40 hisat2 2.2.1
      python3 erna/erna_app.py build_index "Homo sapiens" GCF_009914755.1 bowtie 2.5.2
      python3 erna/erna_app.py build_index "Homo sapiens" GCF_009914755.1 hisat2 2.2.1
      '''
      if len(args)>=4:
        from pipelines.process.align import Align
        specie, genome_version, aligner, aligner_version = args[1:]
        c = Align()
        c.index_builder(specie, genome_version, aligner, aligner_version)
        return c.build_genome_index()

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

  print('wrong arguments')

if __name__ == '__main__':
  args = sys.argv[1:]
  main(args)
