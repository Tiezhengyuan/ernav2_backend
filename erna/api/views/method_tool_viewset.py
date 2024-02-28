import json
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
            return MethodTool.objects.filter(method=method_name)
        return MethodTool.objects.all()
    
    @action(detail=False, methods=['get'])
    def refresh(self, request):
        res = MethodTool.objects.refresh()
        return Response({'created': len(res)})

    @action(detail=False, methods=['get'])
    def method_tools(self, request):
        '''
        used by Vue
        '''
        res = {}
        for obj in MethodTool.objects.all():
            method_name = obj.method.method_name
            if method_name not in res:
                res[method_name] = []
            item = {}
            value = {
                'method_tool_id': obj.id,
                'method_name': obj.method.method_name,
            }
            if obj.tool:
                default_params = json.loads(obj.tool.default_params) if \
                    obj.tool.default_params else {}
                value.update({
                    'tool_name': obj.tool.tool_name,
                    'exe_name': obj.tool.exe_name,
                    'version': obj.tool.version,
                    'params': default_params,
                })
                item = {
                    'text': f"{obj.tool.tool_name}_{obj.tool.version}",
                    'value': value,
                }
            else:
                if obj.method.default_params:
                    value['params'] = json.loads(obj.method.default_params) \
                        if obj.method.default_params else {},
                item = {
                    'text': 'built-in',
                    'value': value,
                }
            # include empty methods of which no tools are defined.
            res[method_name].append(item)
        return Response(res)
            