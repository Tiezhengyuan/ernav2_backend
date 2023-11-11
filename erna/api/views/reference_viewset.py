
from rest_framework import viewsets, response, permissions, decorators

from rna_seq.models import Specie, Genome, Annotation, Reference
from api.serializers import ReferenceSerializer


class ReferenceViewSet(viewsets.ModelViewSet):
    queryset = Reference.objects.all()
    serializer_class = ReferenceSerializer
    permission_classes = [permissions.IsAuthenticated]
