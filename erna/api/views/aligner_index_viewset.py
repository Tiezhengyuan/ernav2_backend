
from rest_framework import viewsets, response, permissions, decorators

from rna_seq.models import Specie, Genome, Annotation, AlignerIndex
from api.serializers import AlignerIndexSerializer


class AlignerIndexViewSet(viewsets.ModelViewSet):
    queryset = AlignerIndex.objects.all()
    serializer_class = AlignerIndexSerializer
    permission_classes = [permissions.IsAuthenticated]
