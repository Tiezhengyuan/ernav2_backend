'''
store external tools
'''
import os
import sys
from django.db import models
from django.conf import settings
from pipelines.utils.dir import Dir

TOOL_EXE = {
    'bowtie': ['bowtie2', 'bowtie2-build', 'bowtie2-inspect',
        'bowtie', 'bowtie-build', 'bowtie-inspect'],
    'fastqc': ['fastqc',],
    'hisat2': ['hisat2', 'hisat2-build', 'hisat2-inspect'],
    'minimap2': ['minimap2',],
    'samtools': ['samtools',],
    'stringtie': ['stringtie'],
    'tophat': ['tophat'],
}

class ToolManager(models.Manager):
    def refresh(self):
        '''
        refresh table Tool
        '''
        # delete all data
        self.all().delete()
        # add tools
        res = []
        externals_dir = settings.EXTERNALS_DIR
        for tool_name in os.listdir(externals_dir):
            versions_path = os.path.join(externals_dir, tool_name)
            for version in os.listdir(versions_path):
                tool_path = os.path.join(versions_path, version)
                for exe_name in TOOL_EXE.get(tool_name, []):
                    exe_path = os.path.join(tool_path, exe_name)
                    if os.path.isfile(exe_path):
                        tool = self.update_or_create(
                            tool_name=tool_name,
                            version=version,
                            exe_name=exe_name,
                            defaults = {
                                'tool_path': tool_path,
                                'exe_path': exe_path,
                            }
                        )
                        res.append(tool)
        return res

    def get_tool(self, exe_name:str, version:str=None):
        if version:
            return self.get(exe_name=exe_name, version=version)
        return self.filter(exe_name=exe_name).last()

class Tool(models.Model):
    # required
    tool_name = models.CharField(
        max_length=32,
        verbose_name='Tool name',
    )
    version = models.CharField(max_length=32)
    exe_name = models.CharField(
        max_length=32,
        verbose_name='Executable name',
    )
    tool_path = models.CharField(
        max_length=256,
        verbose_name='Path of tool',
    )
    exe_path = models.CharField(
        max_length=256,
        verbose_name='Path of the executable',
    )
    # optional
    default_params = models.CharField(
        max_length=1028,
        blank=True,
        null=True,
        verbose_name='Default parameters'
    )

    objects = ToolManager()

    class Meta:
        app_label = 'rna_seq'
        unique_together = ['exe_name', 'version']
        ordering = ['tool_name', 'version', 'exe_name']
    
    def __str__(self):
        return self.tool_name
    
    def add_tool(self, tool_info:dict):
        '''
        new tool would be added
        '''
        if 'tool_name' not in tool_info:
            raise ValueError('Name of the tool should be defined')
        if 'tool_path' not in tool_info:
            raise ValueError('path of the tool should be defined')
        dup = Tool.objects.filter(
            tool_name=tool_info['tool_name'],
            version=tool_info['version']
        )
        if len(dup) > 0:
            raise ValueError('duplicate on the name and version')
        # post a new tool
        tool_obj = Tool.objects.create(**tool_info)
        return tool_obj