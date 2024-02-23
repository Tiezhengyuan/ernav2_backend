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
            return Method.objects.filter(method_name=method_name)
        return Method.objects.all()
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = Method.objects.refresh()
        return Response({'created': len(res)})

    @action(detail=False, methods=['get'])
    def method_names(self, request):
        res = []
        methods = Method.objects.all()
        for method in methods:
            obj = {
                'method_name': method.method_name,
                'label': method.method_name.title().replace('_', ' '),
                'component': method.method_name.title().replace('_', ''),
            }
            res.append(obj)
        return Response(res)
    
