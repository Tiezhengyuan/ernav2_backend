import os
import json
from django.db import models
from django.conf import settings

from .specie import Specie


class GenomeManager(models.Manager):
    def get_genome(self, data_source:str, specie_name:str, version:str=None):
        try:
            specie = Specie.objects.get(specie_name=specie_name)
            if version is None:
                return self.filter(data_source=data_source, specie=specie).last()
            return self.get(data_source=data_source, specie=specie, version=version)
        except Exception as e:
            pass
        return None
    
    def get_versions(self, data_source:str, specie_name:str):
        specie = Specie.objects.get(specie_name=specie_name)
        query = self.filter(data_source=data_source, specie=specie)
        return [i.version for i in query]

    def version_genomes(self):
        '''
        group genomes by data_source and specie_name
        '''
        res = {}
        for genome in self.all():
            data_source = genome.data_source if genome.data_source else 'other'
            if data_source not in res:
                res[data_source] = {}
            specie_name = genome.specie.specie_name
            if specie_name not in res[data_source]:
                res[data_source][specie_name] = []
            item = {
                'value': genome.version,
                'text': genome.version,
            }
            res[data_source][specie_name].append(item)
        return res

    def get_ftp_path(self, data_source:str, specie:str, version:str=None):
        obj = self.filter(data_source=data_source, specie=specie, \
            version=version).first()
        return obj.ftp_path

    def load_genome(self, data:dict):
        '''
        post new genome or update metadata of the genome
        '''
        if 'specie' in data and 'version' in data:
            if 'metadata' in data:
                data['metadata'] = json.dumps(data['metadata'])
            existing = Genome.objects.filter(specie=data['specie'],\
                version=data['version'])
            if not existing:
                data['specie'] = Specie.objects.get(specie_name=data['specie'])
                return Genome.objects.create(**data)
            else:
                if 'metadata' in data:
                    return existing.update(metadata=data['metadata'])
        return None

    def refresh(self):
        '''
        update Genome if genome was downloaded and is available
        '''
        res = []
        ref_dir = settings.REFERENCES_DIR
        for data_source in os.listdir(ref_dir):
            genome_dir = os.path.join(ref_dir, data_source, 'genome')
            if os.path.isdir(genome_dir):
                for specie in os.listdir(genome_dir):
                    specie_dir = os.path.join(genome_dir, specie)
                    for version in os.listdir(specie_dir):
                        local_path = os.path.join(specie_dir, version)
                        obj = self.filter(specie=specie, version=version)\
                            .update(is_ready=True, local_path=local_path)
                        res.append(obj)
        return res
         
class Genome(models.Model):
    '''
    Genome DNA, one chromosome one fasta file
    '''
    specie = models.ForeignKey(
        Specie,
        on_delete=models.CASCADE
    )
    data_source = models.CharField(
        max_length=10,
        default = "NCBI",
        choices=[
            ('NCBI', 'NCBI'),
            ('ENSEMBL', 'ENSEMBL'),
            ('other', 'other'),
        ] 
    )
    version = models.CharField(max_length=56)
    ftp_path = models.CharField(
        max_length=512,
        blank=True,
        null=True
    )
    local_path = models.CharField(
        max_length=1028,
        blank=True,
        null=True
    )
    # str type from json format
    metadata = models.CharField(
        max_length=1256,
        blank=True,
        null=True
    )
    is_ready = models.BooleanField(default=False)

    objects = GenomeManager()

    class Meta:
        app_label = 'rna_seq'
        unique_together = ('specie', 'version')
        ordering = ['specie', 'version']

    def __str__(self):
        return self.specie.specie_name