from rest_framework import serializers

from rna_seq.models import *


class TaskTreeSerializer(serializers.ModelSerializer):
    class Meta:
        model = TaskTree
        fields = '__all__'
        depth = 1

class TaskSerializer(serializers.ModelSerializer):
    child_tasks = TaskTreeSerializer(many=True)
    
    class Meta:
        model = Task
        fields = '__all__'

class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = '__all__'
        depth = 1

class ProjectUserSerializer(serializers.ModelSerializer):
    class Meata:
        model = ProjectUser
        fields = '__all__'


class TaskExecutionSerializer(serializers.ModelSerializer):
    class Meta:
        model = TaskExecution
        fields = '__all__'

class ExecutionTreeSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExecutionTree
        fields = '__all__'
        depth = 1

class TaskAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = TaskAnnotation
        fields = '__all__'
