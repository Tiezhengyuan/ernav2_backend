'''
molecular annotations derive from downloaded annotations
parsed into .fa and .json
'''
from django.db import models
from django.contrib.contenttypes.fields import GenericRelation

from .annotation import Annotation
from .genome import Genome
from .aligner_index import AlignerIndex

class MolecularAnnotationManager(models.Manager):
    def load(self, meta:list) -> list:
        res = []
        for item in meta:
            annotations = Annotation.objects.filter(file_path=item['infile'])
            data = {}
            if item['file_format'] == 'json':
                data['annot_json'] = item['outfile']
            else:
                data['file_path'] = item['outfile']
            if len(annotations) >= 1 and data:
                obj = self.update_or_create(
                    genome = annotations[0].genome,
                    molecular_type = item['molecular_type'],
                    defaults = data,
                )
                res.append(obj)
        return res

class MolecularAnnotation(models.Model):
    genome = models.ForeignKey(
        Genome,
        related_name = 'mol_annots',
        on_delete = models.CASCADE
    )
    molecular_type = models.CharField(max_length=24)
    # NOTE: the field names showed as the below should
    # be compatible with that used in other models
    file_path = models.CharField(
        max_length=512,
        blank = True,
        null = True
    )
    annot_json = models.CharField(
        max_length=512,
        blank = True,
        null = True
    )

    indexes = GenericRelation(AlignerIndex)
    objects = MolecularAnnotationManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('genome', 'molecular_type')
        unique_together = ('genome', 'molecular_type')
    
    def get_index_path(self, tool=None):
        '''
        args: tool_query: {'exec_name':<>, 'version':<>}
        '''
        if tool:
            for aligner_index in self.indexes.all():
                if aligner_index.tool == tool and aligner_index.index_path:
                    return aligner_index.index_path
        return None
