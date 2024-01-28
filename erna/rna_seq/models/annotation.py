'''
annotation files determined by NCBI or Ensemble
file format could be fata, gtf/gff, gbk etc.
'''
import os
from django.db import models
from django.conf import settings
from django.contrib.contenttypes.fields import GenericRelation

from .genome import Genome
from .aligner_index import AlignerIndex
from .tool import Tool

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

    def genome_annot(self, genome:Genome, file_format:str):
        '''
        retrieve annotation of genome DNA in FASTA
        '''
        annotations = self.filter(genome=genome)
        for annot in annotations:
            if annot.annot_type == 'genomic' and annot.file_format == file_format:
                return annot
        return None


class Annotation(models.Model):
    genome = models.ForeignKey(
        'rna_seq.Genome',
        related_name = 'annots',
        on_delete = models.CASCADE
    )
    file_path = models.CharField(max_length = 512)
    file_format = models.CharField(max_length = 24)
    annot_type = models.CharField(
        max_length = 48,
        null=True,
        blank=True
    )

    indexes = GenericRelation(AlignerIndex)
    objects = AnnotationManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('genome', 'file_path')
    
    def get_index_path(self, tool=None):
        '''
        args: tool_query: {'exec_name':<>, 'version':<>}
        '''
        if tool:
            for aligner_index in self.indexes.all():
                if aligner_index.tool == tool and aligner_index.index_path:
                    return aligner_index.index_path
        return None
