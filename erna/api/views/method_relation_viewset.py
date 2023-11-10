from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import MethodRelation
from api.serializers import MethodRelationSerializer

class MethodRelationViewSet(viewsets.ModelViewSet):
    serializer_class = MethodRelationSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        method_name = self.request.query_params.get('method_name', None)
        if method_name is not None:
            return MethodRelation.objects.get_method_tools(method_name)
        return MethodRelation.objects.all()
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = MethodRelation.objects.refresh()
        return Response({'created': len(res)})
