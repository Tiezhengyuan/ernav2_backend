import os
from django.db import models
from django.conf import settings

from .genome import Genome

FILE_TYPES = [
    ('DNA', 'DNA'),
    ('transcript', 'transcript'),
    ('protein', 'protein'),
    ('CDS', 'CDS'),
    ('RNA', 'RNA'),
    ('other', 'other'),
]

class AnnotationManager(models.Manager):
    def load_annotations(self, genome, local_files):
        res = []
        for file_path in local_files:
            data = {
                'file_format': os.path.splitext(file_path)[1],
                'annot_type': self.annot_type(file_path)
            }
            obj = self.update_or_create(
                genome=genome,
                file_path=file_path,
                defaults=data
            )
            res.append(obj)
        return res

    def annot_type(self, file_path:str):
        '''
        orders of if-statements do matter
        '''
        if 'rna' in file_path:
            return 'RNA'
        if 'cds' in file_path:
            return 'CDS'
        if file_path.endswith('fna'):
            return 'DNA'
        if file_path.endswith('gtf') or file_path.endswith('gff'):
            return 'transcript'
        return 'other'

class Annotation(models.Model):
    genome = models.ForeignKey(
        Genome,
        related_name = 'annots',
        on_delete = models.CASCADE
    )
    file_path = models.CharField(max_length = 256)
    file_format = models.CharField(
        max_length = 8,
        blank = True,
        null = True
    )
    annot_type = models.CharField(
        max_length = 24,
        default = 'other',
        choices = FILE_TYPES,
    )

    objects = AnnotationManager()

    class Meta:
        app_label = 'rna_seq'
        unique_together = ('genome', 'file_path')
        ordering = ('genome', 'file_path')
    

