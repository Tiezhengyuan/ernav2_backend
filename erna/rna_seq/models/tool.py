'''
store external tools
'''
import os
import sys
from django.db import models
from django.conf import settings
from pipelines.process.utils.dir import Dir


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
        for exe_path in Dir(externals_dir).recrusive_files():
            tool_path = os.readlink(exe_path)
            tool_name = os.path.basename(tool_path)
            version = os.path.basename(os.path.dirname(tool_path))
            tool = self.update_or_create(
                tool_name=tool_name,
                version=version,
                defaults= {'tool_path':tool_path, 'exe_path':exe_path}
            )
            res.append(tool)
        return res

    def get_tool(self, tool_name:str, version:str=None):
        if version:
            return self.get(tool_name=tool_name, version=version)
        return self.filter(tool_name=tool_name).last()

class Tool(models.Model):
    # required
    tool_name = models.CharField(max_length=32)
    version = models.CharField(max_length=32)
    tool_path = models.CharField(max_length=256)
    exe_path = models.CharField(max_length=256)
    # optional
    default_params = models.CharField(
        max_length=1028,
        blank=True,
        null=True,
        verbose_name='default parameters of the tool'
    )

    objects = ToolManager()

    class Meta:
        app_label = 'rna_seq'
        unique_together = ['tool_name', 'version']
        ordering = ['tool_name', 'version']
    
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