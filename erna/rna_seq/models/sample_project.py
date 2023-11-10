import os
from django.db import models
from typing import Iterable
from .project import Project
from .sample_file import SampleFile


class SampleProjectManager(models.Manager):

    def iterate_sample_files(self, project_id)->Iterable:
        project = Project.objects.get_project_by_project_id(project_id)
        for sample_project_obj in self.filter(project=project):
            file_obj = sample_project_obj.sample_file
            batch_name = file_obj.sample.batch_name
            yield batch_name, file_obj

    def get_sample_files(self, project_id)->dict:
        '''
        project_id ~ batch_name -> [file_path,...]
        '''
        res = {}
        sample_files = self.iterate_sample_files(project_id)
        for batch_name, file_obj in sample_files:
            if batch_name not in res:
                res[batch_name] = []
            res[batch_name].append(file_obj.file_path)
        return res
    
    def get_samples(self, project_id)->dict:
        '''
        project_id ~ batch_name -> [sample,...]
        '''
        res = {}
        sample_files = self.iterate_sample_files(project_id)
        for batch_name, sample_file in sample_files:
            if batch_name not in res:
                res[batch_name] = []
            sample_name = sample_file.sample.sample_name
            if sample_name not in res[batch_name]:
                res[batch_name].append(sample_file.sample)
        return res


    def get_sample_names(self, project_id)->dict:
        '''
        project_id ~ batch_name -> [sample_name,...]
        '''
        res = {}
        project_samples = self.get_samples(project_id)
        for batch_name, samples in project_samples.items():
            res[batch_name] = [s.sample_name for s in samples]
        return res

    def get_project_sample_files(self, project_id):
        project_files = SampleProject.objects.filter(project_id=project_id)
        res = []
        for pf in project_files:
            study_name = pf.sample_file.sample.study_name
            sample_name = pf.sample_file.sample.sample_name
            full_path = os.path.join(pf.sample_file.raw_data.file_path,
                pf.sample_file.raw_data.file_name)
            tag = 0
            for el in res:
                if study_name == el['study_name'] and sample_name == el['sample_name']:
                    el['raw_data'].append(full_path)
                    tag = 1
                    break
            if tag == 0:
                item = {
                    'study_name': study_name,
                    'sample_name': sample_name,
                    'raw_data': [full_path,],
                }
                res.append(item)
        return res

    def load_project_sample_file(self, data_iters:list):
        '''
        args: data_iters[{},{}]
        '''
        res = []
        for data in data_iters:
            if 'project_id' in data:
                project = Project.objects.get_project_by_project_id(data['project_id'])
                if 'batch_name' in data and 'sample_name' in data:
                    sample_files = SampleFile.objects.get_sample_files(
                        data['batch_name'], data['sample_name'])
                    for sample_file in sample_files:
                        sp_obj = self.create(project=project, sample_file=sample_file)
                        res.append(sp_obj)
        return res

    def add_sample_file(self, project_sample):
        existing = self.filter(**project_sample)
        print(existing[0].project, existing[0].sample_file)
        if not existing:
            obj = self.create(**project_sample)
            obj.save()
            return obj
        return None


class SampleProject(models.Model):
    '''
    one project may include many samples
    Samples may come from same or different batch
    '''
    project = models.ForeignKey(
        Project,
        related_name='sample_projects',
        on_delete=models.CASCADE
    )
    sample_file = models.ForeignKey(
        SampleFile,
        related_name='sample_files',
        on_delete=models.CASCADE
    )
    
    objects = SampleProjectManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ['project', 'sample_file']
    
    def __str__(self):
        return self.project.project_id