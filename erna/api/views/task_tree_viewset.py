from django.shortcuts import render
from rest_framework import viewsets, permissions

from rna_seq.models import *
from api.serializers import *

class TaskTreeViewSet(viewsets.ModelViewSet):
    queryset = TaskTree.objects.all()
    serializer_class = TaskTreeSerializer
    permission_classes = [permissions.IsAuthenticated]


