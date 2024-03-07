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
    def front_rnas(self, request):
        type_rnas = {}
        for rna in RNA.objects.all():
            if rna.annot_type not in type_rnas:
                type_rnas[rna.annot_type] = []
            _specie_name = rna.specie.specie_name if rna.specie else 'unknown'
            item = {
                'value': rna.id,
                'text': f"{rna.data_source}_{_specie_name}",
            }
            type_rnas[rna.annot_type].append(item)

        res = {
            'type_rnas': type_rnas,
            'rna_types': [{'value':i, 'text':i} for i in list(type_rnas)],
        }
        return Response(res)
