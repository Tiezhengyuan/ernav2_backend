from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.renderers import JSONRenderer
from django.http import Http404, JsonResponse
from rest_framework import viewsets, permissions
from rna_seq.models import Project, SampleProject, SampleFile
from api.serializers import SampleProjectSerializer
       
class SampleProjectViewSet(viewsets.ModelViewSet):
    serializer_class = SampleProjectSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        project_id = self.request.query_params.get('project_id')
        if project_id:
            project = Project.objects.get(project_id=project_id)
            return SampleProject.objects.filter(project=project)
        return SampleProject.objects.all()

    @action(detail=False, methods=['post'])
    def load_project_samples(self, request):
        projec_id = request.data.get('project_id')
        samples = request.data.get('samples')
        res = SampleProject.objects.load_data(projec_id, samples)
        return Response({'count': len(res)})

    @action(detail=False, methods=['delete'])
    def delete_project_samples(self, request):
        project_id = self.request.query_params.get('project_id')
        res = {}
        if project_id:
            project = Project.objects.get(project_id=project_id)
            samples = SampleProject.objects.filter(project=project)
            res['deleted'] = samples.count()
            samples.delete()
        return Response(res)

    @action(detail=False, methods=['get'])
    def project_data(self, request):
        '''
        reponsee is used by UI
        '''
        project_id = self.request.query_params.get('project_id')
        study_name = self.request.query_params.get('study_name')
        res = {}
        if project_id:
            res['project_samples']= SampleProject.objects.get_project_samples(project_id)
            if study_name:
                res['unassigned_samples'] = SampleProject.objects.get_unassigned_samples(
                    project_id, study_name
                )
        return Response(res)


    


