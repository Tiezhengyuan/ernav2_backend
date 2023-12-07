'''
non-coding RNA
'''
import os
from django.db import models
from typing import Iterable

from .specie import Specie

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
                specie = specie,
                rna_type = item['rna_type'],
                defaults = {
                    'fa_path': item['fa_path'],
                    'database': item['database'],
                }
            )


class RNA(models.Model):
    specie = models.ForeignKey(
        'rna_seq.Specie',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    fa_path = models.CharField(max_length=512)
    rna_type = models.CharField(max_length=16)
    database = models.CharField(
        max_length=10,
        null=True,
        blank=True,
    )

    objects = RNAManager()

    class Meta:
        app_label = 'rna_seq'
    
    def index_path(self, tool_name:str, tool_version:str):
        index_dir_path = os.path.join(os.path.dirname(self.fa_path), 'index')
        file_prefix = os.path.splitext(os.path.basename(self.fa_path))[0]        
        index_path = os.path.join(
            index_dir_path,
            f"{tool_name}_{tool_version}_{file_prefix}"
        )
        return index_dir_path, index_path
