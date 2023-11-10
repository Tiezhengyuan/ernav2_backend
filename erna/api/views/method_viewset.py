from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import Method
from api.serializers import MethodSerializer

class MethodViewSet(viewsets.ModelViewSet):
    serializer_class = MethodSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        method_name = self.request.query_params.get('method_name', None)
        if method_name is not None:
            return Method.objects.get_method_tools(method_name)
        return Method.objects.all()
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = Method.objects.refresh()
        return Response({'created': len(res)})
