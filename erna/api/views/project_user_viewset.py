from django.shortcuts import render
from rest_framework import viewsets, permissions

from rna_seq.models import *
from api.serializers import *


class ProjectUserViewSet(viewsets.ModelViewSet):
    queryset = ProjectUser.objects.all()
    serializer_class = ProjectUserSerializer
    permission_classes = [permissions.IsAuthenticated,]


