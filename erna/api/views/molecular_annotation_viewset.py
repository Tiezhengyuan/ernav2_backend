from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import MolecularAnnotation
from api.serializers import MolecularAnnotationSerializer


class MolecularAnnotationViewSet(viewsets.ModelViewSet):
    serializer_class = MolecularAnnotationSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        annotation = self.request.query_params.get('annotation')
        if annotation:
            return MolecularAnnotation.objects.filter(annotation=annotation)
        return MolecularAnnotation.objects.all()
    
    @action(detail=False, methods=['delete'])
    def delete_genome_annot(self, request):
        annotation = self.request.query_params.get('annotation')
        res = MolecularAnnotation.objects.filter(annotation=annotation)
        count = res.count()
        res.delete()
        return Response({'deleted': count})

    @action(detail=False, methods=['delete'])
    def delete_all(self, request):
        res = MolecularAnnotation.objects.all()
        count = res.count()
        res.delete()
        return Response({'deleted': count})
    
