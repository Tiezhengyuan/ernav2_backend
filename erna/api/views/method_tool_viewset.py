from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import MethodTool
from api.serializers import MethodToolSerializer

class MethodToolViewSet(viewsets.ModelViewSet):
    serializer_class = MethodToolSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        method_name = self.request.query_params.get('method_name', None)
        if method_name is not None:
            return MethodTool.objects.get_method_tools(method_name)
        return MethodTool.objects.all()
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = MethodTool.objects.refresh()
        return Response({'created': len(res)})
