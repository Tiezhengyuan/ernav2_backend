from rest_framework.decorators import action
from rest_framework import viewsets, permissions, status, serializers
from rest_framework.response import Response

from rna_seq.models import NonCodingRNA
from api.serializers import NonCodingRNASerializer


class NonCodingRNAViewSet(viewsets.ModelViewSet):
    queryset = NonCodingRNA.objects.all()
    serializer_class = NonCodingRNASerializer
    permission_classes = [permissions.IsAuthenticated,]

    @action(detail=False, methods=['get'])
    def count(self, request):
        count = NonCodingRNA.objects.count()
        return Response({'count': count})
