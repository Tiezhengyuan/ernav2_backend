from typing import Iterable
from django.db import models

class SpecieManager(models.Manager):
    def load_species(self, obj_iter:Iterable):
        res = []
        names = ['taxid', 'organism_name', 'group']
        for _, summary in obj_iter:
            data = dict([(n, summary[n]) for n in names])
            name1, name2 = data['organism_name'].split(' ')[:2]
            data['abbreviation'] = name1[0].lower() + name2[:2].lower()
            data['organism_name'] = data['organism_name'].replace(' ', '_')
            obj = self.update_or_create(
                organism_name = data['organism_name'],
                defaults = data
            )
            res.append(obj)
        return res


class Specie(models.Model):
    organism_name = models.CharField(
        primary_key=True,
        max_length=256
    )
    group = models.CharField(
        max_length=56,
        blank=True,
        null=True
    )
    taxid = models.IntegerField(
        blank=True,
        null=True
    )
    abbreviation = models.CharField(
        max_length=10,
        blank=True,
        null=True
    )

    objects = SpecieManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ['group', 'organism_name',]

    def __str__(self):
         return self.organism_name