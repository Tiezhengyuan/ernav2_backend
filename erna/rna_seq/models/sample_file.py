'''
one file one record.
'''
import re
from django.db import models
from django.conf import settings
import os
from typing import Iterable

from .sample import Sample
from .raw_data import RawData

class SampleFileManager(models.Manager):
    def parse_sample_rawdata(self, study_names:list, batch_names:list):
        '''
        parse samples ~ raw_data
        one sample may have multiple raw_data
        args: Sample.study_name
        args: RawData.batch_name
        '''
        # collect samples and raw data
        samples = []
        for study_name in study_names:
            samples += Sample.objects.filter(study_name=study_name)
        raw = []
        for batch_name in batch_names:
            raw += RawData.objects.filter(batch_name=batch_name)
        
        # parse Sample~RawData
        sample_files = []
        for sample in samples:
            for raw_data in raw:
                if raw_data.file_name.startswith(sample.sample_name):
                    obj = self.update_or_create(sample=sample, raw_data=raw_data)
                    sample_files.append(obj)
        return sample_files

    def load_sample_files(self, sample_files:list):
        '''
        '''
        res = []
        for data in sample_files:
            raw_data = RawData.objects.get(
                batch_name=data['batch_name'],
                file_path=data['file_path'],
                file_name=data['file_name'],
            )
            sample = Sample.objects.get(
                study_name=data['study_name'],
                sample_name=data['sample_name']
            )
            obj = self.update_or_create(sample=sample, raw_data=raw_data)
            res.append(obj)
        return res
    
    def iterate_raw_data(self, study_name:str)->Iterable:
        '''
        Iteration of all files defined in RawData
        given a study name defined in Sample
        '''
        samples = Sample.objects.filter(study_name=study_name)
        for sample in samples:
            sample_files = self.filter(sample=sample)
            for sample_file in sample_files:
                raw_data = sample_file.raw_data
                yield sample, raw_data, sample_file
                 
    def get_study_files(self, study_name:str):
        '''
        get relations of sample~sample file given a study_name
        '''
        res = {}
        iter = self.iterate_raw_data(study_name)
        for sample, raw_data, sample_file in iter:
            sample_name = sample.sample_name
            if sample_name not in res:
                res[sample_name] = {
                    'sample_file_id': sample_file.id,
                    'sample_id': sample.id,
                    'sample_name': sample_name,
                    'raw_data': [],
                }
            sample_file = raw_data.to_dict()
            res[sample_name]['raw_data'].append(sample_file)
        return res


    def add_sample_file(self, sample_file:dict):
        '''
        post a new record
        '''
        existing = self.model.objects.filter(
            sample_id=sample_file['sample'],
            raw_data_id=sample_file['raw_data']
        )
        if not existing:
            # SampleFile
            file_obj = self.model.objects.create(
                sample_id=sample_file['sample'],
                raw_data_id=sample_file['raw_data']
            )
            file_obj.save()
            return file_obj
        return None
    
    def get_sample_files(self, sample_name:str):
        '''
        one sample may include 0-many files
        '''
        try:
            sample = Sample.objects.get(sample_name=sample_name)
            return self.model.objects.filter(sample=sample.id)
        except Exception as e:
            pass
        return []

    def detect_unparsed_data(self, study_name, reg):
        '''
        one raw data match one sample
        one sample may match to 1-many raw data
        '''
        # get unparsed raw data
        raw_data = RawData.objects.all()
        parsed_raw_data = {i[1] for i in self.iterate_raw_data(study_name)}
        unparsed_raw_data = set(raw_data).difference(parsed_raw_data)
        # print(len(unparsed_raw_data), len(parsed_raw_data), len(raw_data))

        # try to parsing raw data and sample name
        sample_data = []
        for raw_data in unparsed_raw_data:
            file_name = raw_data.file_name
            for sample in Sample.objects.filter(study_name=study_name):
                target = reg.replace('<S>', sample.sample_name)
                if re.search(target, file_name):
                    pair = {
                        'sample_id': sample.id,
                        'sample_name': sample.sample_name,
                        'raw_data_id': raw_data.id,
                        'full_path': os.path.join(raw_data.file_path, raw_data.file_name),
                        'file_name': raw_data.file_name,
                        'batch_name': raw_data.batch_name,
                        # used for UI, key name is fixed
                        '_showDetails': False,
                    }
                    sample_data.append(pair)
                    break
        sample_data = sorted(sample_data, key=lambda x: x['sample_name'])
        return sample_data

class SampleFile(models.Model):
    # one sample may include multiple files
    sample = models.ForeignKey(
        Sample,
        on_delete = models.CASCADE
    )
    # path of raw data
    raw_data = models.ForeignKey(
        RawData,
        on_delete = models.CASCADE
    )

    objects = SampleFileManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('sample', 'raw_data')

    def get_sample(self):
        res = Sample.objects.get(id=self.sample.id)
        return res.to_dict()
    
    def get_raw_data(self):
        res = RawData.objects.get(id=self.raw_data.id)
        return res.to_dict()

