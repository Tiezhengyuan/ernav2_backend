from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import Annotation
from api.serializers import AnnotationSerializer


class AnnotationViewSet(viewsets.ModelViewSet):
    serializer_class = AnnotationSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        genome = self.request.query_params.get('genome')
        if genome:
            return Annotation.objects.filter(genome=genome)
        return Annotation.objects.all()
    
    @action(detail=False, methods=['delete'])
    def delete_genome_annot(self, request):
        genome = request.query_params.get('genome')
        res = Annotation.objects.filter(genome=genome)
        count = res.count()
        res.delete()
        return Response({'deleted': count})

    @action(detail=False, methods=['delete'])
    def delete_all(self, request):
        res = Project.objects.all()
        count = res.count()
        res.delete()
        return Response({'deleted': count})
    
