'''
example:
    python3 erna/manage.py shell < erna/scripts/demo_mirnaseq_LW.py
'''
from rna_seq.models import *
from commons.models import CustomUser

user = CustomUser.objects.get(pk=1)

print("Cretae project...") 
project_id = "P00005"
project_data = {
    "project_name": "test_iterative_mirna_seq",
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
sample_names = ['C1S', 'P25']
sample_data = [{'study_name':study_name, 'sample_name':s,} for s in sample_names]
samples = Sample.objects.load_samples(user, sample_data)
print(samples)

print('Load RawData...')
batch_names = ['miRNA_2023',]
sample_files = SampleFile.objects.parse_sample_rawdata([study_name,], batch_names)

print('Update SampleProject...')
res = SampleProject.objects.load_project_sample_files(
    project_id, [s.id for s,_ in sample_files]
)
print(res)

print('Add tasks...')
specie_name = "Homo_sapiens"
builder = {
    'tool_name': 'bowtie',
    'exe_name': 'bowtie2-build',
    'version': '2.5.2',
}
aligner = {
    'tool_name': 'bowtie',
    'exe_name': 'bowtie2',
    'version': '2.5.2',
}
tasks_data = [
    {
        'task_id': 'T01',
        'task_name': 'trim 3-end sequences',
        'method_name': 'trim_sequences',
        'params': {
            'adapter_3end': 'AACTGTAGGCACCATCAAT',
        },
    },
    # mature miRNA
    {
        'task_id': 'T02',
        'task_name': 'miRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'RNA',
            'query': {
                'specie': specie_name,
                'annot_type': 'miRNA_mature',
                'database': 'miRBase',
            }
        },
    },
    {
        'task_id': 'T03',
        'task_name': 'align miRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T04',
        'task_name': 'count miRNA',
        'method_name': 'count_reads',
    },
    # hairpin miRNA
    {
        'task_id': 'T05',
        'task_name': 'hairpin miRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'RNA',
            'query': {
                'specie': specie_name,
                'annot_type': 'miRNA_hairpin',
            }
        },
    },
    {
        'task_id': 'T06',
        'task_name': 'align hairpin miRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T07',
        'task_name': 'count hairpin miRNA',
        'method_name': 'count_reads',
    },
    # piwiRNA
    {
        'task_id': 'T08',
        'task_name': 'piwiRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'RNA',
            'query': {
                'specie': specie_name,
                'annot_type': 'piwiRNA',
            }
        },
    },
    {
        'task_id': 'T09',
        'task_name': 'align piwiRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T10',
        'task_name': 'count piwiRNA',
        'method_name': 'count_reads',
    },
    # lncRNA
    {
        'task_id': 'T11',
        'task_name': 'lncRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'RNA',
            'query': {
                'specie': specie_name,
                'annot_type': 'lncRNA',
            }
        },
    },
    {
        'task_id': 'T12',
        'task_name': 'align lncRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T13',
        'task_name': 'count lncRNA',
        'method_name': 'count_reads',
    },
    # rRNA
    {
        'task_id': 'T14',
        'task_name': '5SrRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'RNA',
            'query': {'annot_type': '5SrRNA',}
        },
    },
    {
        'task_id': 'T15',
        'task_name': 'align 5SrRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T16',
        'task_name': 'count 5SrRNA',
        'method_name': 'count_reads',
    },
    # tRNA
    {
        'task_id': 'T17',
        'task_name': 'tRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'RNA',
            'query': {'annot_type': 'tRNA',}
        },
    },
    {
        'task_id': 'T18',
        'task_name': 'align tRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T19',
        'task_name': 'count tRNA',
        'method_name': 'count_reads',
    },
    # mRNA
    {
        'task_id': 'T20',
        'task_name': 'mRNA',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'Annotation',
            'query': {
                "file_path": "/home/yuan/bio/references/NCBI/genome/Homo_sapiens/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_rna.fna",
            }
        },
    },
    {
        'task_id': 'T21',
        'task_name': 'align mRNA',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T22',
        'task_name': 'count mRNA',
        'method_name': 'count_reads',
    },
    # pseudo
    {
        'task_id': 'T23',
        'task_name': 'pseudo genes',
        'method_name': 'build_index',
        'tool': builder,
        'params': {
            'model': 'Annotation',
            'query': {
                "file_path": "/home/yuan/bio/references/NCBI/genome/Homo_sapiens/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_pseudo_without_product.fna",
            }
        },
    },
    {
        'task_id': 'T24',
        'task_name': 'align pseudo genes',
        'method_name': 'align_short_reads',
        'tool': aligner,
    },
    {
        'task_id': 'T25',
        'task_name': 'count pseudo genes',
        'method_name': 'count_reads',
    },

    # merge all read counts
    {
        'task_id': 'T30',
        'task_name': 'merge RC',
        'method_name': 'merge_read_counts',
    },
]
# Task.objects.filter(project_id=project_id).delete()
tasks = Task.objects.load_tasks(project_id, tasks_data)
print(tasks)


'''
               T00
   /   /    /   |    \   \   \  \   \
T01  T02  T05  T08  T11 T14 T17 T20 T23
  \  /     |    |    |   |   |   |   |
   T03     |    |    |   |   |   |   |
    |     /     |    |   |   |   |   |
   T04   /      /    |   |   |   |   |
    |\  /      /     |   |   |   |   |
    | T06     /      /   |   |   |   |
    |  |     /      /    |   |   |   |
    | T07   /      /     |   |   |   |
    |  |\  /      /      /   |   |   |
    |  | T09     /      /    |   |   |
    |  |  |     /      /     |   |   |
    |  | T10   /      /      /   |   |
    |  |  | \ /      /      /    |   |
    |  |  | T12     /      /     |   |
    |  |  |  |     /      /      /   |
    |  |  | T13   /      /      /    |
    |  |  |  | \ /      /      /     |
    |  |  |  | T15     /      /      |
    |  |  |  |  |     /      /       /
    |  |  |  | T16   /      /       /
    |  |  |  |  | \ /      /       /
    |  |  |  |  | T18     /       /
    |  |  |  |  |  |     /       /
    |  |  |  |  | T19   /       /
    |  |  |  |  |  | \ /       /
    |  |  |  |  |  | T21      /
    |  |  |  |  |  |  |      /
    |  |  |  |  |  | T22    /
    |  |  |  |  |  |  | \  /
    \  |  |  |  |  |  | T24
     \ \  |  |  |  /  /  |
      \ \ \  |  / /  /  T25
       \ \ \ | / /  /  /
        \ \ \|/ /  /  /
            T30
'''
print('Add Task Tree...')
task_pair = [
    ('T00', 'T01'),
    ('T00', 'T02'), ('T00', 'T05'), ('T00', 'T08'), ('T00', 'T11'),
    ('T00', 'T14'), ('T00', 'T17'), ('T00', 'T20'), ('T00', 'T23'),
    # iterative
    ('T01', 'T03'), ('T02', 'T03'), ('T03', 'T04'),
    ('T04', 'T06'), ('T05', 'T06'), ('T06', 'T07'),
    ('T07', 'T09'), ('T08', 'T09'), ('T09', 'T10'),
    ('T10', 'T12'), ('T11', 'T12'), ('T12', 'T13'),
    ('T13', 'T15'), ('T14', 'T15'), ('T15', 'T16'),
    ('T16', 'T18'), ('T17', 'T18'), ('T18', 'T19'),
    ('T19', 'T21'), ('T20', 'T21'), ('T21', 'T22'),
    ('T22', 'T24'), ('T23', 'T24'), ('T24', 'T25'),
    # merge
    ('T04', 'T30'), ('T07', 'T30'), ('T10', 'T30'), ('T13', 'T30'),
    ('T16', 'T30'), ('T19', 'T30'), ('T22', 'T30'), ('T25', 'T30'),
]
tasks_tree = TaskTree.objects.load_tasks_tree(project_id, task_pair)
print(tasks_tree)

