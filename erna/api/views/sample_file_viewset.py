from rest_framework.decorators import action
from rest_framework.response import Response
from django.http import Http404
from rest_framework import viewsets, permissions
from rna_seq.models import Sample, SampleFile
from api.serializers import SampleFileSerializer
       
class SampleFileViewSet(viewsets.ModelViewSet):
    serializer_class = SampleFileSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        '''
        1. all data
        2. all files given a sample_name
        '''
        sample_name = self.request.query_params.get('sample_name')
        if sample_name is not None:
            return SampleFile.objects.get_sample_files(sample_name)
        return SampleFile.objects.all()
    
    @action(detail=False, methods=['get'])
    def count(self, request):
        sample_name = self.request.query_params.get('sample_name')
        res = SampleFile.objects.all() if sample_name is None else \
            SampleFile.objects.get_sample_files(sample_name)
        count = res.count()
        return Response({'count': count})

    @action(detail=False, methods=['get'])
    def study_files(self, request):
        '''
        given study name, get all files (raw_data)
        '''
        study_name = self.request.query_params.get('study_name')
        if study_name:
            study_files = []
            for sample in Sample.objects.filter(study_name=study_name):
                item = sample.to_dict()
                sample_files = SampleFile.objects.filter(sample=sample)
                item['raw_data'] = [i.get_raw_data() for i in sample_files]
                study_files.append(item)
            return Response(study_files)
        return []

    @action(detail=False, methods=['get'])
    def unparsed_data(self, request):
        '''
        parse samples with unparsed raw data
        '''
        study_name = self.request.query_params.get('study_name', '')
        batch_name = self.request.query_params.get('batch_name')
        reg = self.request.query_params.get('reg', '<S>')
        res = SampleFile.objects.detect_unparsed_data(study_name, reg, batch_name)
        return Response(res)
    
    @action(detail=False, methods=['post'])
    def parse_sample_files(self, request):
        '''
        parse sample with raw data by sample.id and raw_data.id
        [
            {"sample": 1, "raw_data": 3},
            {"sample": 1, "raw_data": 4}
        ]
        '''
        res = {
            'created': [],
            'failed': [],
            'skipped': [],
        }
        sample_files = request.data
        for sample_file in sample_files:
            obj = SampleFile.objects.add_sample_file(sample_file)
            if obj:
                res['created'].append(obj.id)
            else:
                res['skipped'].append(sample_file)
        return Response(res)