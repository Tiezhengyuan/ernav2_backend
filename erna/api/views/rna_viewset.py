import os
from rest_framework.decorators import action
from rest_framework import viewsets, permissions, status, serializers
from rest_framework.response import Response

from rna_seq.models import RNA
from api.serializers import RNASerializer


class RNAViewSet(viewsets.ModelViewSet):
    queryset = RNA.objects.all()
    serializer_class = RNASerializer
    permission_classes = [permissions.IsAuthenticated,]

    @action(detail=False, methods=['get'])
    def count(self, request):
        count = RNA.objects.count()
        return Response({'count': count})

    @action(detail=False, methods=['get'])
    def type_rnas(self, request):
        res = {}
        for rna in RNA.objects.all():
            if rna.annot_type not in res:
                res[rna.annot_type] = []
            item = {
                'value': rna.id,
                'text': os.path.basename(rna.file_path),
            }
            res[rna.annot_type].append(item)
        return Response(res)
