from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.renderers import JSONRenderer
from django.http import Http404, JsonResponse
from rest_framework import viewsets, permissions
from rna_seq.models import SampleProject, SampleFile
from api.serializers import SampleProjectSerializer
       
class SampleProjectViewSet(viewsets.ModelViewSet):
    serializer_class = SampleProjectSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        project_id = self.request.query_params.get('project_id')
        project_files = SampleProject.objects.filter(project_id=project_id) \
            if project_id else SampleProject.objects.all()
        return project_files

    @action(detail=False, methods=['get'])
    def project_sample_files(self, request):
        '''
        given project_id, retrieve all sample, raw_data
        '''
        project_id = self.request.query_params.get('project_id')
        if project_id:
            res = SampleProject.objects.get_project_sample_files(project_id)
            return Response(res)
        return Response({})

    @action(detail=False, methods=['get'])
    def unassigned_sample_files(self, request):
        '''
        given project_id and study_name,
        retrieve sample/raw_data which are unassigned
        '''
        project_id = self.request.query_params.get('project_id')
        project_files = SampleProject.objects.get_project_sample_files(\
            project_id) if project_id else []
        # 
        res = []
        study_name = self.request.query_params.get('study_name')
        if study_name:
            study_files = SampleFile.objects.get_study_files(study_name).values()
            for sf in study_files:
                tag = 0
                for pf in project_files:
                    if study_name == pf['study_name'] and sf['sample_name'] == pf['sample_name']:
                        tag = 1
                        break
                if tag == 0:
                    item = {
                        'study_name': study_name,
                        'sample_file_id': sf['sample_file_id'],
                        'sample_name':sf['sample_name'],
                        'raw_data': [i['full_path'] for i in sf['raw_data']],
                    }
                    res.append(item)
        return Response(res)

    

    @action(detail=False, methods=['delete'])
    def remove_sample_files(self, request):
        '''
        give a project:
            remove all sample files
            remove one certain samplefile
        '''
        res = {'deleted': 0}
        project_id = self.request.query_params.get('project_id')
        sample_file = self.request.query_params.get('sample_file')
        if project_id:
            obj = SampleProject.objects.filter(project_id=project_id, \
                    sample_file=sample_file) if sample_file else \
                SampleProject.objects.filter(project_id=project_id)
            res['deleted'] = obj.count()
            obj.delete()
        return Response(res)

    @action(detail=False, methods=['post'])
    def load_sample_files(self, request):
        res = {'created': 0, 'skipped': 0}
        project_samples = request.data
        for project_sample in project_samples:
            obj = SampleProject.objects.add_sample_file(project_sample)
            if obj:
                res['created'] += 1
            else:
                res['skipped'] += 1
        return Response(res)
