'''
Note: order of the classes does matter
'''
from rest_framework import serializers
from rna_seq.models import *

class MethodToolSerializer(serializers.ModelSerializer):
    class Meta:
        model = MethodTool
        fields = '__all__'
        depth = 1
    # def get_queryset(self, *args):
    #     pass

class ToolSerializer(serializers.ModelSerializer):
    methods = MethodToolSerializer(many=True)

    class Meta:
        model = Tool
        fields = ['tool_name', 'version', 'tool_path', \
            'exe_path', 'default_params', 'methods']

class MethodRelationSerializer(serializers.ModelSerializer):
    class Meta:
        model = MethodRelation
        fields = '__all__'

class MethodSerializer(serializers.ModelSerializer):
    children = MethodRelationSerializer(many=True)
    tools = MethodToolSerializer(many=True)

    class Meta:
        model = Method
        fields = ['method_name', 'description', \
            'children', 'tools']

class PipelineSerializer(serializers.ModelSerializer):
    class Meta:
        model = Pipeline
        fields = '__all__'
        depth = 2

# class MethodNameSerializer(serializers.ListSerializer):
#     class Meta:
#         model = Method
#         fields = ['method_name',]