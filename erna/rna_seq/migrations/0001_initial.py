# Generated by Django 4.2.6 on 2023-12-08 02:22

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('contenttypes', '0002_remove_content_type_name'),
    ]

    operations = [
        migrations.CreateModel(
            name='Annotation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file_path', models.CharField(max_length=256)),
                ('file_format', models.CharField(blank=True, max_length=8, null=True)),
                ('annot_type', models.CharField(blank=True, max_length=48, null=True)),
            ],
            options={
                'ordering': ('genome', 'file_path'),
            },
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('data_source', models.CharField(choices=[('NCBI', 'NCBI'), ('ENSEMBL', 'ENSEMBL'), ('other', 'other')], default='NCBI', max_length=10)),
                ('version', models.CharField(max_length=56)),
                ('ftp_path', models.CharField(blank=True, max_length=512, null=True)),
                ('local_path', models.CharField(blank=True, max_length=1028, null=True)),
                ('metadata', models.CharField(blank=True, max_length=1256, null=True)),
                ('is_ready', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['specie', 'version'],
            },
        ),
        migrations.CreateModel(
            name='Method',
            fields=[
                ('method_name', models.CharField(max_length=96, primary_key=True, serialize=False)),
                ('description', models.CharField(blank=True, max_length=1028, null=True)),
            ],
            options={
                'ordering': ['method_name'],
            },
        ),
        migrations.CreateModel(
            name='MethodTool',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('method', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='tools', to='rna_seq.method')),
            ],
            options={
                'ordering': ['method', 'tool'],
            },
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('project_id', models.CharField(max_length=10, primary_key=True, serialize=False, verbose_name='project ID')),
                ('status', models.CharField(choices=[('active', 'active'), ('ready', 'ready'), ('locked', 'locked'), ('deleted', 'deleted')], default='active', max_length=10)),
                ('sequencing', models.CharField(choices=[('mrna-seq', 'mRNA-Seq'), ('mirna-seq', 'miRNA-Seq'), ('scrna-seq', 'scRNA-Seq'), ('other', 'other')], default='mrna-seq', max_length=10, verbose_name='Sequencing technique')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('project_name', models.CharField(blank=True, max_length=50, null=True, verbose_name='Project name')),
                ('description', models.CharField(blank=True, max_length=526, null=True)),
                ('genome', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='rna_seq.genome')),
                ('owner', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL, verbose_name='owner identified by user_id')),
            ],
            options={
                'ordering': ['project_id', 'owner'],
            },
        ),
        migrations.CreateModel(
            name='RawData',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file_path', models.CharField(max_length=512)),
                ('file_name', models.CharField(max_length=128)),
                ('file_format', models.CharField(blank=True, max_length=10, null=True)),
                ('file_type', models.CharField(blank=True, max_length=10, null=True)),
                ('batch_name', models.CharField(blank=True, max_length=20, null=True)),
                ('parsed', models.BooleanField(default=False)),
            ],
            options={
                'ordering': ['batch_name', 'file_path', 'file_name'],
                'unique_together': {('batch_name', 'file_path', 'file_name')},
            },
        ),
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('study_name', models.CharField(max_length=50)),
                ('sample_name', models.CharField(max_length=100)),
                ('metadata', models.CharField(blank=True, max_length=3000, null=True)),
                ('creator', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ('study_name', 'sample_name'),
            },
        ),
        migrations.CreateModel(
            name='SampleFile',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('raw_data', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.rawdata')),
                ('sample', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.sample')),
            ],
            options={
                'ordering': ('sample', 'raw_data'),
            },
        ),
        migrations.CreateModel(
            name='Specie',
            fields=[
                ('specie_name', models.CharField(max_length=128, primary_key=True, serialize=False)),
                ('organism_name', models.CharField(max_length=128)),
                ('other_names', models.CharField(blank=True, max_length=128, null=True)),
                ('group', models.CharField(blank=True, max_length=56, null=True)),
                ('taxid', models.IntegerField(blank=True, null=True)),
                ('abbreviation', models.CharField(blank=True, max_length=10, null=True)),
            ],
            options={
                'ordering': ['group', 'organism_name'],
            },
        ),
        migrations.CreateModel(
            name='Task',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('task_id', models.CharField(max_length=10, verbose_name='Task ID')),
                ('is_ready', models.BooleanField(default=False)),
                ('task_name', models.CharField(blank=True, max_length=100, null=True, verbose_name='Task Name')),
                ('params', models.CharField(blank=True, max_length=1256, null=True, verbose_name='Parameters (in json string format)')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('modified_time', models.DateTimeField(auto_now=True)),
                ('annotation', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='rna_seq.annotation')),
                ('method_tool', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='rna_seq.methodtool', verbose_name='Method and Tool')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='tasks', to='rna_seq.project', verbose_name='Project ID')),
            ],
            options={
                'ordering': ('project', 'task_id'),
            },
        ),
        migrations.CreateModel(
            name='Tool',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tool_name', models.CharField(max_length=32, verbose_name='Tool name')),
                ('version', models.CharField(max_length=32)),
                ('exe_name', models.CharField(max_length=32, verbose_name='Executable name')),
                ('tool_path', models.CharField(max_length=256, verbose_name='Path of tool')),
                ('exe_path', models.CharField(max_length=256, verbose_name='Path of the executable')),
                ('default_params', models.CharField(blank=True, max_length=1028, null=True, verbose_name='Default parameters')),
            ],
            options={
                'ordering': ['tool_name', 'version', 'exe_name'],
                'unique_together': {('exe_name', 'version')},
            },
        ),
        migrations.CreateModel(
            name='TaskTree',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('child', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='child_tasks', to='rna_seq.task', verbose_name='Child task')),
                ('task', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.task', verbose_name='Task')),
            ],
            options={
                'ordering': ('task', 'child'),
            },
        ),
        migrations.CreateModel(
            name='TaskExecution',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('status', models.CharField(choices=[('suspend', 'suspend'), ('ready', 'ready'), ('pause', 'pause'), ('run', 'run'), ('finish', 'finish'), ('skip', 'skip'), ('fail', 'fail')], default='suspend', max_length=10)),
                ('start_time', models.DateTimeField(blank=True, null=True)),
                ('end_time', models.DateTimeField(blank=True, null=True)),
                ('output', models.CharField(blank=True, max_length=5000, null=True, verbose_name='Path of output files')),
                ('task', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='task_executions', to='rna_seq.task')),
            ],
        ),
        migrations.AddField(
            model_name='task',
            name='task_execution',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='task_execution', to='rna_seq.taskexecution'),
        ),
        migrations.CreateModel(
            name='SampleProject',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sample_projects', to='rna_seq.project')),
                ('sample_file', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='sample_files', to='rna_seq.samplefile')),
            ],
            options={
                'ordering': ['project', 'sample_file'],
            },
        ),
        migrations.CreateModel(
            name='RNA',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file_path', models.CharField(max_length=512)),
                ('file_format', models.CharField(blank=True, max_length=8, null=True)),
                ('annot_type', models.CharField(blank=True, max_length=48, null=True)),
                ('database', models.CharField(blank=True, max_length=10, null=True)),
                ('specie', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='rna_seq.specie')),
            ],
        ),
        migrations.CreateModel(
            name='ProjectUser',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.project', verbose_name='Executed project identified by project_id')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL, verbose_name='executor identified by user_id')),
            ],
            options={
                'ordering': ['project', 'user'],
            },
        ),
        migrations.CreateModel(
            name='Pipeline',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pipeline_name', models.CharField(max_length=128)),
                ('next_step', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='next_steps', to='rna_seq.methodtool')),
                ('step', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='steps', to='rna_seq.methodtool')),
            ],
            options={
                'ordering': ['pipeline_name', 'step'],
            },
        ),
        migrations.AddField(
            model_name='methodtool',
            name='tool',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='methods', to='rna_seq.tool'),
        ),
        migrations.CreateModel(
            name='MethodRelation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('child', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='children', to='rna_seq.method')),
                ('method', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.method')),
            ],
            options={
                'ordering': ['method', 'child'],
            },
        ),
        migrations.AddField(
            model_name='genome',
            name='specie',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.specie'),
        ),
        migrations.CreateModel(
            name='ExecutionTree',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('child_execution', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='child_executions', to='rna_seq.taskexecution', verbose_name='Child of executed task')),
                ('execution', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.taskexecution', verbose_name='executed task')),
            ],
        ),
        migrations.AddField(
            model_name='annotation',
            name='genome',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='annots', to='rna_seq.genome'),
        ),
        migrations.CreateModel(
            name='AlignerIndex',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('index_path', models.CharField(max_length=256, verbose_name='index path used by aligner')),
                ('object_id', models.PositiveIntegerField()),
                ('content_type', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='contenttypes.contenttype')),
                ('tool', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.tool')),
            ],
            options={
                'ordering': ['tool', 'index_path'],
            },
        ),
        migrations.CreateModel(
            name='TaskAnnotation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('annotation', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.annotation')),
                ('task', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='rna_seq.task')),
            ],
            options={
                'ordering': ['annotation', 'task'],
                'unique_together': {('annotation', 'task')},
            },
        ),
        migrations.AlterUniqueTogether(
            name='task',
            unique_together={('project', 'task_id')},
        ),
        migrations.AlterUniqueTogether(
            name='methodtool',
            unique_together={('method', 'tool')},
        ),
        migrations.AlterUniqueTogether(
            name='genome',
            unique_together={('specie', 'version')},
        ),
    ]
