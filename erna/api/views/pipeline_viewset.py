from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import Pipeline
from api.serializers import PipelineSerializer

class PipelineViewSet(viewsets.ModelViewSet):
    serializer_class = PipelineSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        pipeline_name = self.request.query_params.get('pipeline_name')
        if pipeline_name:
            return Pipeline.objects.filter(pipeline_name=pipeline_name)
        return Pipeline.objects.all()
    
    @action(detail=False, methods=['delete'])
    def delete(self, request):
        pipeline_name = self.request.query_params.get('pipeline_name')
        if pipeline_name:
            return Pipeline.objects.filter(pipeline_name=pipeline_name).delete()
