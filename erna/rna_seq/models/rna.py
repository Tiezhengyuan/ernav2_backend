'''
non-coding RNA
'''
import os
from typing import Iterable
from django.db import models
from django.contrib.contenttypes.fields import GenericRelation


from .specie import Specie
from .aligner_index import AlignerIndex
from .tool import Tool

class RNAManager(models.Manager):
    def load_data(self, data:Iterable):
        for item in data:
            specie = None
            if item.get('specie_name'):
                species = Specie.objects.filter(specie_name=item['specie_name'])
                if species:
                    if len(species) == 1:
                        specie = species[0]
                    else:
                        print(f"Duplicated specie detectd. {item}")
                else:
                    new_data = {
                        'specie_name': item['specie_name'],
                        'organism_name': item['organism_name'],
                        'other_names': item.get('other_names'),
                        'abbreviation': item['abb'],
                    }
                    specie = Specie.objects.create(**new_data)
            
            #Note: include complete data if specie is None
            self.update_or_create(
                file_path = item['fa_path'],
                defaults = {
                    'specie': specie,
                    'data_source': item['data_source'],
                    'annot_type': item['annot_type'],
                    'annot_json': item.get('annot_json'),
                }
            )


class RNA(models.Model):
    file_path = models.CharField(max_length=512)
    specie = models.ForeignKey(
        'rna_seq.Specie',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    # path of json file
    annot_json = models.CharField(
        max_length=512,
        null=True,
        blank=True
    )
    # miRNA, tRNA, rRNA, mRNA etc
    annot_type = models.CharField(
        max_length = 48,
        null=True,
        blank=True
    )
    data_source = models.CharField(
        max_length=10,
        null=True,
        blank=True,
    )

    indexes = GenericRelation(AlignerIndex)
    objects = RNAManager()

    class Meta:
        app_label = 'rna_seq'
    
    def get_index_path(self, tool=None):
        '''
        args: tool_query: {'exec_name':<>, 'version':<>}
        '''
        if tool:
            for aligner_index in self.indexes.all():
                if aligner_index.tool == tool:
                    return aligner_index.index_path
        return None
