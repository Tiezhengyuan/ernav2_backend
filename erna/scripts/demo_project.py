'''
example:
    python3 erna/manage.py shell < erna/scripts/demo_project.py
'''
import json
from rna_seq.models import *
from annot.models import *
from commons.models import CustomUser


project_id = "P00001"

# cretae project
new_project = {
    "project_id": project_id,
    "project_name": "test_mrna_seq",
    "description": "for testing mRNA-seq",
    "status": "A",
    "sequencing": "M",
    'owner': CustomUser.objects.get(pk=1),
    'specie': Specie.objects.get(organism_name='Homo sapiens'),
}
Project.objects.filter(project_id=project_id).delete()
project = Project.objects.create(**new_project)
print(project)

# add tasks
tasks_data = [
    {
        'task_id': 'T01',
        'task_name': 'build index',
        'method_name': 'build_index',
        'tool_name': 'hisat2-build',
        'params': {
            'data_source': 'NCBI',
            'specie': 'Homo sapiens',
            'version': 'GCF_000001405.40',
        },
        'child': ['T02'],
    },
    {
        'task_id': 'T02',
        'task_name': '',
        'method_name': 'align_transcriptome',
        'tool_name': 'hisat2',
        'params': {
            'index_path': 'aaa',
        },
        'child': ['T03'],
    },
    {
        'task_id': 'T03',
        'task_name': '',
        'method_name': 'assemble_transcripts',
        'tool_name': 'stringtie',
        'params': {},
        'child': ['T04'],
    },
    {
        'task_id': 'T04',
        'task_name': '',
        'method_name': 'count_reads',
        'parent': ['T03'],
    },
    {
        'task_name': '',
        'method_name': 'quality_control',
        "tool_name": "fastqc",
    },
    # empty task without method_name
    {
        'task_id': 'T06',
    },
    # update task
    {
        'task_id': 'T06',
        'task_name': 'test update',
    },
]
tasks = Task.objects.load_tasks(project_id, tasks_data)
print(tasks)

tasks_tree = TaskTree.objects.load_tasks_tree(project_id, tasks_data)
print(tasks_tree)


samples = {
    'sample_1': {
        'R1': 'sample_1_R1.fastq',
        'R2': 'sample_1_R2.fastq',
    },
    'sample_2': {
        'R1': 'sample_2_R1.fastq',
        'R2': 'sample_2_R2.fastq',
    }
}
