'''
scheduled task execution
'''

from django.shortcuts import render
from rest_framework import viewsets, permissions

from rna_seq.models import ExecutionTree
from api.serializers import ExecutionTreeSerializer


class ExecutionTreeViewSet(viewsets.ModelViewSet):
    queryset = ExecutionTree.objects.all()
    serializer_class = ExecutionTreeSerializer
    permission_classes = [permissions.IsAuthenticated]

