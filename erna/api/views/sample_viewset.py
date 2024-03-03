from rest_framework.decorators import action
from rest_framework.response import Response
from django.http import Http404
from rest_framework import viewsets, permissions
from rna_seq.models import Sample
from api.serializers import SampleSerializer
       
class SampleViewSet(viewsets.ModelViewSet):
    serializer_class = SampleSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        '''
        1. get list of all samples
        2. get list of samples given study_name/creator
        '''
        config = {}
        study_name = self.request.query_params.get('study_name')
        creator = self.request.query_params.get('creator')
        if study_name is not None:
            config['study_name'] = study_name
        if creator is not None:
            config['creator'] = creator
        if config:
            return Sample.objects.filter(**config)
        return Sample.objects.all()
    
    @action(detail=False, methods=['get'])
    def count(self, request):
        '''
        total number of samples
        '''
        res = Sample.objects.all()
        count = res.count()
        return Response({'count': count})

    @action(detail=False, methods=['get'])
    def count_study(self, request):
        '''
        number of samples per study
        '''
        res = Sample.objects.group_by_study()
        count = []
        for k in res:
            item = {
                'study_name': k,
                'count': len(res[k]),
            }
            count.append(item)
        return Response(count)

    @action(detail=False, methods=['get'])
    def study_names(self, request):
        '''
        study names
        '''
        study_names = Sample.objects.values_list('study_name', flat=True).distinct()
        res = [{'value': s, 'text': s} for s in list(set(study_names))]
        return Response(res)
    
    @action(detail=False, methods=['post'])
    def load_samples(self, request):
        '''
        insert multiple samples given a study_name. for example:
        [
            {"study_name":"prostate", "sample_name":"a", "metadata":{}},
            {"study_name":"prostate", "sample_name":"b", "metadata":{"type":"control"}}
        ]
        '''
        res = Sample.objects.load_samples(request.user, request.data)
        return Response({'loaded': len(res)})

    @action(detail=False, methods=['delete'])
    def delete_study_samples(self, request):
        '''
        deletee sample given a study name
        '''
        study_name = self.request.query_params.get('study_name')
        if study_name is not None:
            res = Sample.objects.delete_study_samples(study_name)
        return Response({'message': f"{res} are deleted."})

    @action(detail=False, methods=['delete'])
    def delete_all(self, request):
        samples = Sample.objects.all()
        res = {'deleted': len(samples)}
        samples.delete()
        return Response(res)