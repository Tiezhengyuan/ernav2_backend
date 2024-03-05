import os
from django.db import models
from typing import Iterable

from .project import Project
from .sample_file import SampleFile
from .sample import Sample


class SampleProjectManager(models.Manager):

    def load_data(self, project_id:str, samples:list):
        res = []
        project = Project.objects.get(project_id=project_id)
        for item in samples:
            sample = Sample.objects.get(
                study_name=item['study_name'],
                sample_name=item['sample_name'],
            )
            obj = self.update_or_create(
                project = project,
                sample = sample,
            )
            res.append(obj)
        return res

    def iterate_samples(self, project_id)->Iterable:
        project = Project.objects.get(project_id=project_id)
        for obj in self.filter(project=project):
            yield obj.sample

    def get_project_samples(self, project_id) -> list:
        '''
        '''
        res = []
        project = Project.objects.get(project_id=project_id)
        for obj in self.filter(project=project):
            item = obj.sample.to_dict()
            item['sample_project_id'] = obj.id
            sample_files = SampleFile.objects.filter(sample=obj.sample)
            item['raw_data'] = [i.get_raw_data() for i in sample_files]
            res.append(item)
        return res
    
    def get_unassigned_samples(self, project_id, study_name):
        project = Project.objects.get(project_id=project_id)
        _assigned = {i.sample for i in self.filter(project=project)}
        study_samples = {i for i in Sample.objects.filter(study_name = study_name)}

        _unassigned = study_samples.difference(_assigned)
        res = []
        for sample in _unassigned:
            item = sample.to_dict()
            sample_files = SampleFile.objects.filter(sample=sample)
            item['raw_data'] = [i.get_raw_data() for i in sample_files]
            res.append(item)
        return res



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
    sample = models.ForeignKey(
        Sample,
        related_name='project_samples',
        on_delete=models.CASCADE
    )
    
    objects = SampleProjectManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ['project', 'sample']
    