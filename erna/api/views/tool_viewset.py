from rest_framework import viewsets, permissions
from rest_framework.decorators import action
from rest_framework.response import Response

from rna_seq.models import Tool
from api.serializers import ToolSerializer

class ToolViewSet(viewsets.ModelViewSet):
    serializer_class = ToolSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        params = {}
        tool_name = self.request.query_params.get('tool_name')
        version = self.request.query_params.get('version')
        if tool_name:
            params['tool_name'] = tool_name
        if version:
            params['version'] = version
        if params:
            return Tool.objects.filter(**params)
        return Tool.objects.all()
    
    @action(detail=False, methods=['get'])
    def count(self, request):
        count = Tool.objects.all().count()
        return Response({'count': count})

    @action(detail=False, methods=['delete'])
    def delete_all(self, request):
        res = Tool.objects.all()
        count = res.count()
        res.delete()
        return Response({'deleted': count})
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = Tool.objects.refresh()
        return Response({'created': len(res)})