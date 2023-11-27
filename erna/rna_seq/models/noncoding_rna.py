'''
non-coding RNA
'''
from django.db import models
from typing import Iterable
from .specie import Specie

class NonCodingRNAManager(models.Manager):
    def load_data(self, data:Iterable):
        for item in data:
            specie = Specie.objects.filter(abbreviation=item['abb'])
            if specie:
                if len(specie) == 1:
                    self.update_or_create(
                        specie = specie[0],
                        rna_type = item['rna_type'],
                        defaults = {
                            'fa_path': item['fa_path'],
                            'database': item['db'],
                        }
                    )
                else:
                    print(f"Duplicated specie detectd. {item}")
            else:
                pass
                # print(f"No specie detectd. {item}")

class NonCodingRNA(models.Model):
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

    objects = NonCodingRNAManager()

    class Meta:
        app_label = 'rna_seq'
