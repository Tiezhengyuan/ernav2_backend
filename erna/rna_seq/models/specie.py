from django.db import models

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

    class Meta:
        app_label = 'rna_seq'
        ordering = ['group', 'organism_name',]

    def __str__(self):
         return self.organism_name