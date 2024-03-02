import json
from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import Method, MethodTool, MethodRelation
from api.serializers import MethodToolSerializer

class MethodToolViewSet(viewsets.ModelViewSet):
    serializer_class = MethodToolSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        method_name = self.request.query_params.get('method_name', None)
        if method_name is not None:
            return MethodTool.objects.filter(method=method_name)
        return MethodTool.objects.all()
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = MethodTool.objects.refresh()
        return Response({'created': len(res)})

    @action(detail=False, methods=['get'])
    def front_methods(self, request):
        '''
        responsed is used by SelectMethod.Vue
        '''
        res = {
            'method_names': Method.objects.method_names(),
            'method_tools': MethodTool.objects.get_method_tools(),
            'method_parents': MethodRelation.objects.get_parents(),
        }
        return Response(res)
            