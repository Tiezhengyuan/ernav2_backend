'''
example:
    python3 erna/manage.py shell < erna/scripts/demo_mirnaseq.py
'''
from rna_seq.models import *
from commons.models import CustomUser

user = CustomUser.objects.get(pk=1)

print("Cretae project...") 
project_id = "P00002"
project_data = {
    "project_name": "test_mirna_seq",
    "description": "test miRNA-seq pipeline",
    "status": "active",
    "sequencing": "mirna-seq",
    'owner': user,
}
project = Project.objects.update_or_create(
    project_id = project_id,
    defaults = project_data
)
print(project)

print('Load samples...')
study_name = 'test_mirnaseq'
sample_names = ['LW_AB_1', 'LW_AI_1', 'LW_AN_1']
sample_data = [{'study_name':study_name, 'sample_name':s, 'metadata':{},} \
    for s in sample_names]
samples = Sample.objects.load_samples(user, sample_data)
print(samples)

print('Load RawData...')
batch_names = ['demo_mirnaseq',]
sample_files = SampleFile.objects.parse_sample_rawdata([study_name,], batch_names)

print('Update SampleProject...')
res = SampleProject.objects.load_project_sample_files(
    project_id, [s.id for s,_ in sample_files]
)
print(res)

print('Add tasks...')
tasks_data = [
    {
        'task_id': 'T01',
        'method_name': 'trim_sequences',
        'params': {
            'adapter_3end': 'TATGC',
        },
    },
    {
        'task_id': 'T02',
        'method_name': 'build_index',
        'tool': {
            'tool_name': 'bowtie',
            'exe_name': 'bowtie2-build',
            'version': '2.5.2',
        },
        'params': {
            'model': 'NonCodingRNA',
            'query': {
                'specie': "Homo_sapiens",
                'rna_type': 'mature',
            }
        },
    },
    {
        'task_id': 'T03',
        'method_name': 'align_short_reads',
        'tool': {
            'tool_name': 'bowtie',
            'exe_name': 'bowtie2',
            'version': '2.5.2',
        },
    },
    {
        'task_id': 'T04',
        'method_name': 'convert_format',
        'tool': {
            'tool_name': 'samtools',
            'exe_name': 'samtools',
            'version': '1.18',
        },
    },
    {
        'task_id': 'T05',
        'method_name': 'count_reads',
    },
    {
        'task_id': 'T06',
        'method_name': 'quality_control',
        'tool': {
            "tool_name": "fastqc",
            "exe_name": "fastqc",
            'version': '0.12.1'
        }
    },
]
Task.objects.filter(project_id=project_id).delete()
tasks = Task.objects.load_tasks(project_id, tasks_data)
print(tasks)

print('Add Task Tree...')
task_pair = [
    ('T00', 'T01'),('T00', 'T02'),('T00', 'T06'),
    ('T01', 'T03'),('T02', 'T03'),
    ('T03', 'T04'),('T03', 'T05'),
]
tasks_tree = TaskTree.objects.load_tasks_tree(project_id, task_pair)
print(tasks_tree)

