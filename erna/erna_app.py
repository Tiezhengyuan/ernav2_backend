'''
entrance into pipelines
'''
import argparse
import os
import sys
from dotenv import load_dotenv
load_dotenv()
DIR_PROJECT = os.path.dirname(os.path.dirname(__file__))
sys.path.append(DIR_PROJECT)

import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'erna.settings')
django.setup()

from pipelines.process import ProcessRawData, ProcessGenome, \
  ProcessNCRNA, ExecuteTasks, Align, TrimAdapter


def main(args):
  match args.get('method', ''):
    case 'scan_raw_data':
      '''
      Get batch_name and path from raw data 
      example:
      python3 erna_app.py -m scan_raw_data
      '''
      return ProcessRawData().scan_raw_data()

    case 'refresh_raw_data':
      '''
      scan raw data and load to db.RawData
      example:
      python3 erna_app.py -m scan_raw_data
      '''
      return ProcessRawData().refresh_raw_data()

    case 'assembly_summary':
      '''
      download assembly_summary.txt
      example: python3 erna/erna_app.py -m assembly_summary -d NCBI
      '''
      data_source = args.get('data_source')
      if data_source:
        return ProcessGenome(data_source).retrieve_assembly_summary()

    case 'download_genome':
      '''
      download refseq of genome fro FTP
      example:
      python3 erna_app.py download_genome -d NCBI -d Homo_sapiens -g GCF_000001405.40
      '''
      data_source = args.get('data_source')
      specie = args.get('specie')
      genome_version = args.get('genome_version')
      if data_source and specie and genome_version:
        return ProcessGenome(data_source, specie, genome_version).download_genome()

    case 'load_mirbase':
      '''
      download data from miRBase and load them into eRNA
      example: python3 erna_app.py -m load_mirbase
      '''
      return ProcessNCRNA().load_mirbase()
      
    case 'execute_tasks':
      '''
      example:
      python3 erna/erna_app.py -m execute_tasks -p P00001 -t T01
      '''
      project_id = args.get('project_id')
      task_id = args.get('task_id')
      chain = args.get('chain')
      if project_id:
        return ExecuteTasks(project_id, task_id, chain)()

    case 'build_genome_index':
      '''
      build index for aligner namely bowtie
      example: python3 erna/erna_app.py -m build_index -s Homo_sapiens \
        -g GCF_000001405.40 --tool_name hisat2 --tool_version 2.2.1
      '''
      specie = args.get('specie')
      genome_version = args.get('genome_version')
      tool_name = args.get('tool_name')
      tool_version = args.get('tool_version')
      if specie and genome_version and tool_name and tool_version:
        c = Align()
        c.index_builder(specie, genome_version, tool_name, tool_version)
        return c.build_genome_index()

    case 'trim_adapter':
      '''
      trim adapter sequences from reads for miRNA-seq
      example: 
        python3 erna_app.py -m trim_adapter -a ATGGCG \
          -i /home/yuan/bio/erna_v2/pipelines/tests/data/reads_1.fq \
          -o /home/yuan/bio/erna_v2/pipelines/tests/temp/reads_1_trimmed.fq
      '''
      adapter_3end = args.get('adapter_3end')
      input =args.get('input')
      output =args.get('output')
      if adapter_3end and input and output:
        params =  {
          'adapter_3end': adapter_3end,
          'input': input,
          'output': output,
        }
        return TrimAdapter(params)()

  print('wrong arguments')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='eRNAv2: Analyze RNA-Seq data.')
  
  # required
  parser.add_argument('-m', '--method',
    help='Method or certain analytic step')

  # optional relatd to method
  parser.add_argument('-d', '--data_source',
    help="Data source of annotations namely NCBI")
  parser.add_argument('-s', '--specie',
    help="Specie name defined in data source")
  parser.add_argument('-g', '--genome_version',
    help="Genome version or accession given a genome defined in data source")
  parser.add_argument('--tool_name',
    help="third-party tool or built-in tool")
  parser.add_argument('--tool_version',
    help="tool version given a tool")
  parser.add_argument('-p', '--project_id', help="Project ID")
  parser.add_argument('-t', '--task_id', default='T00',
    help="Task ID given project_id")
  parser.add_argument('-c', '--chain', default=False,
    help="Execute current task and all remaining tasks followed by this task.")
  parser.add_argument('-a', '--adapter_3end',
    help='adapter sequences')
  parser.add_argument('-i', '--input', help='Input path')
  parser.add_argument('-o', '--output', help='Output path')
  
  # general optional
  parser.add_argument('-f', '--force_run', default=False, action='store_false', 
    help="Force the app to execute task and overwrite current results.")
    
  # pass arguments
  args = parser.parse_args()
  main(vars(args))
