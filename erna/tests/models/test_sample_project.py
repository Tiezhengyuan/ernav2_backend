from django.test import TestCase, override_settings
from ddt import ddt, data, unpack
from unittest.mock import patch
import os
from rna_seq.models import Project
from commons.models import User
from sample.models import Sample, SampleFile, SampleProject

@ddt
class TestSampleFile(TestCase):
    def setUp(self):
        user = User.objects.create(
            user_name="tester",
            email="a@abc.com",
            password="pass"
        )
        self.project = Project.objects.create(
            project_id='P00001',
            project_name='test1',
            owner=user
        )
        sample = Sample.objects.create(
            batch_name='batch1',
            sample_name='sample1',
            creator=user
        )
        self.R1 = SampleFile.objects.create(
            sample = sample,
            file_name = 'a_R1.fq',
            file_format = 'fastq',
            file_type = 'R1',
            batch_no = 'PS1'
        )
        self.R2 = SampleFile.objects.create(
            sample = sample,
            file_name = 'a_R2.fq',
            file_format = 'fastq',
            file_type = 'R2',
            batch_no = 'PS1',
            exist = True
        )
        SampleProject.objects.create(project=self.project, sample_file=self.R1)
        SampleProject.objects.create(project=self.project, sample_file=self.R2)

    def test_CRUD(self):
        res = SampleProject.objects.get(project=self.project, sample_file=self.R1)
        assert str(res) == 'P00001'

        res = [ i.sample_file.file_name for i in SampleProject.objects.all()]
        assert res == ['a_R1.fq', 'a_R2.fq']
        
        res = [ i.sample_file.file_name for i in \
            SampleProject.objects.filter(project=self.project)]
        assert res == ['a_R1.fq', 'a_R2.fq']

        res = SampleProject.objects.filter(sample_file=self.R1)\
            .update(sample_file=self.R2)
        assert res == 1

        res = SampleProject.objects.filter(sample_file=self.R2).delete()
        assert res[0] == 2
        res = SampleProject.objects.filter(sample_file=self.R2).delete()
        assert res[0] == 0
    
    @data(
        ['P00001', 'batch1', 2],
    )
    @unpack
    def test_get_sample_files(self, project_id, batch_name, expect):
        res = SampleProject.objects.get_sample_files(project_id)
        assert batch_name in res
        assert len(res[batch_name]) == expect

    @data(
        ['P00001', 'batch1', 2],
    )
    @unpack
    def test_get_samples(self, project_id, batch_name, expect):
        res = SampleProject.objects.get_samples(project_id)
        assert len(res[batch_name]) == expect

    @data(
        ['P00001', {'batch1': ['sample1', 'sample1']}],
    )
    @unpack
    def test_get_sample_names(self, project_id, expect):
        res = SampleProject.objects.get_sample_names(project_id)
        assert res == expect
    
    def test_load_project_sample_file(self):
        data=[
            {'project_id': 'P00001', 'batch_name': 'batch1', 'sample_name': 'sample1'},
        ]
        res = SampleProject.objects.load_project_sample_file(data)
        assert len(res) == 2
