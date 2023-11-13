'''
scheduled task execution
'''

from django.shortcuts import render
from rest_framework import viewsets, permissions

from rna_seq.models import TaskAnnotation
from api.serializers import TaskAnnotationSerializer


class TaskAnnotationViewSet(viewsets.ModelViewSet):
    queryset = TaskAnnotation.objects.all()
    serializer_class = TaskAnnotationSerializer
    permission_classes = [permissions.IsAuthenticated]

