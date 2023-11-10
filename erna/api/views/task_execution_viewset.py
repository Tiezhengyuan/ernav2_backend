'''
scheduled task execution
'''

from django.shortcuts import render
from rest_framework import viewsets, permissions

from rna_seq.models import *
from api.serializers import *


class TaskExecutionViewSet(viewsets.ModelViewSet):
    queryset = TaskExecution.objects.all()
    serializer_class = TaskExecutionSerializer
    permission_classes = [permissions.IsAuthenticated]

