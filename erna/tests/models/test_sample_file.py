from django.test import TestCase, override_settings
from ddt import ddt, data, unpack
from unittest.mock import patch
import os
from commons.models import User
from sample.models import Sample, SampleFile

@ddt
class TestSampleFile(TestCase):
    def setUp(self):
        user = User.objects.create(
            user_name="tester",
            email="a@abc.com",
            password="pass"
        )
        self.sample = Sample.objects.create(
            batch_name='batch1',
            sample_name='sample1',
            creator=user
        )
        self.R1 = SampleFile.objects.create(
            sample = self.sample,
            file_name = 'a_R1.fq',
            file_format = 'fastq',
            file_type = 'R1',
            batch_no = 'PS1'
        )
        self.R2 = SampleFile.objects.create(
            sample = self.sample,
            file_name = 'a_R2.fq',
            file_format = 'fastq',
            file_type = 'R2',
            batch_no = 'PS1',
            exist = True
        )

    @data(
        ['batch1', 'sample1', 2],
        ['batch1', 'sample2', 0],
        ['batch2', 'sample1', 0],
    )
    @unpack
    def test_get_sample_files(self, batch_name, sample_name, expect):
        res = SampleFile.objects.get_sample_files(batch_name, sample_name)
        assert len(res) == expect


    def test_CRUD(self):
        res = SampleFile.objects.get(sample=self.sample, file_type='R1')
        assert res.sample.batch_name == 'batch1'
        assert res.sample.sample_name == 'sample1'
        assert res.file_type == 'R1'
        assert res.exist == False

        res = [ i.file_name for i in SampleFile.objects.all()]
        assert res == ['a_R1.fq', 'a_R2.fq']
        
        res = SampleFile.objects.filter(sample=self.sample).update(exist=True)
        assert res == 2

        res = SampleFile.objects.filter(sample=self.sample).delete()
        assert res[0] == 2
        res = SampleFile.objects.filter(sample=self.sample).delete()
        assert res[0] == 0

    @override_settings(DATA_DIR='/mnt/data')
    def test_storage_path(self):
        assert SampleFile().storage_path == '/mnt/data'

    @override_settings(DATA_DIR='data')
    def test_file_path(self):
        expect = os.path.join('data', 'fastq', 'R1', 'a_R1.fq')
        res = self.R1.file_path
        assert res == expect

        res = SampleFile.objects.get(sample=self.sample, \
            file_type='R1').file_path
        assert res == expect
    
    @patch('os.path.isfile')
    @data(
        [True, True],
        [False, False],
    )
    @unpack
    def test_file_exists(self, input, expect, mock_isfile):
        mock_isfile.return_value=input
        res = SampleFile.objects.get(sample=self.sample, \
            file_type='R1').file_exists()
        assert res == expect
        
    
    def test_add_file(self):
        R1 = {
            'sample': self.sample,
            'file_name': 'b_R1.fq',
            'file_format': 'fastq',
            'file_type': 'R1',
            'batch_no': 'PS1'
        }
        res = SampleFile.objects.add_file(R1)
        assert res.file_name == 'b_R1.fq'

        R2 = {
            'sample': self.sample,
            'file_name': 'b_R2.fq',
            'file_format': 'fastq',
            'file_type': 'R2',
            'batch_no': 'PS1'
        }
        res = SampleFile.objects.add_file(R2)
        assert res.file_name == 'b_R2.fq'

        res = SampleFile.objects.filter(file_type='R1')
        assert len(res) == 2

        
    @data(
        ['batch1', None, 2],
        ['batch1', 'R1', 1],
        ['batch1', 'R2', 1],
    )
    @unpack
    def test_get_batch_files(self, batch_name, file_type, expect):
        res = SampleFile.objects.get_batch_files(batch_name, file_type)
        assert 'sample1' in res
        assert len(res['sample1']) == expect

    @data(
        ['batch1', 2],
    )
    @unpack
    def test_check_not_found_files(self, batch_name, expect):
        res = SampleFile.objects.check_not_found_files(batch_name)
        assert 'sample1' in res
        assert len(res['sample1']) == expect
