import os
from django.db import models
from django.conf import settings


class AnnotationManager(models.Manager):
    def load_annotations(self, genome, local_files):
        res = []
        for file_path in local_files:
            extension = os.path.splitext(file_path)[1]
            data = {
                'file_format': extension.replace('.', ''),
                'annot_type': self.annot_type(file_path, genome)
            }
            obj = self.update_or_create(
                genome=genome,
                file_path=file_path,
                defaults=data
            )
            res.append(obj)
        return res

    def annot_type(self, file_path:str, genome:Genome):
        '''
        not file format but type of annotation 
        '''
        file_name = os.path.basename(file_path)
        file_name = os.path.splitext(file_name)[0]
        file_name = file_name.replace(f"{genome.version}_", "")
        return file_name.split('_', 1)[-1]

class Annotation(models.Model):
    genome = models.ForeignKey(
        'rna_seq.Genome',
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
        max_length = 48,
        null=True,
        blank=True
    )

    objects = AnnotationManager()

    class Meta:
        app_label = 'rna_seq'
        unique_together = ('genome', 'file_path')
        ordering = ('genome', 'file_path')
    

