import os
from rest_framework.decorators import action
from rest_framework.response import Response
from django.http import Http404
from rest_framework import viewsets, permissions
from rna_seq.models import RawData
from api.serializers import RawDataSerializer

class RawDataViewSet(viewsets.ModelViewSet):
    serializer_class = RawDataSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        '''
        get all objects or objects given batch_name
        '''
        batch_name = self.request.query_params.get('batch_name', None)
        if batch_name is not None:
            return RawData.objects.get_batch_files(batch_name)
        return RawData.objects.all()

    def destroy(self, request):
        '''
        delete one by id
        '''
        try:
            instance = self.get_object()
            self.perform_destroy(instance)
        except Exception as e:
            return Response({'error': str(e)})
        return Response({'message': "The record of " + \
            f"{instance.file_path}/{instance.file_name} is deleted"})
    
    @action(detail=False, methods=['get'])
    def count(self, request):
        '''
        number of all rows
        '''
        count = RawData.objects.count()
        return Response({'count': count})
        
    @action(detail=False, methods=['get'])
    def batch_names(self, request):
      '''
      get all batch names
      '''
      res = RawData.objects.values_list('batch_name', flat=True).distinct()
      return Response({'batch_names': res})
    
    @action(detail=False, methods=['delete'])
    def delete_all(self, request):
        '''
        delete all data
        '''
        res = RawData.objects.all()
        count = res.count()
        res.delete()
        return Response({'message': f"{count} are deleted."})

    @action(detail=False, methods=['delete'])
    def delete_batch_files(self, request):
        '''
        delete all files given a batch_name
        '''
        batch_name = self.request.query_params.get('batch_name', None)
        count = 0
        if batch_name is not None:
          res = RawData.objects.filter(batch_name=batch_name)
          count = res.count()
          res.delete()
        return Response({'message': f"{count} are deleted."})

    @action(detail=False, methods=['post'])
    def load_batch_data(self, request):
        '''
        add multiple data. for example:
        {
          "batch_name": "A3",
          "full_path_list": [
            "/raw_data/A3/a.fa",
            "/raw_data/A3/b.fa"
          ]
        }
        '''
        batch_name = request.data.get('batch_name')
        full_path_list = request.data.get('full_path_list')
        count = len(full_path_list) if full_path_list else 0
        for full_file_path in full_path_list:
            file_path = os.path.dirname(full_file_path)
            file_name = os.path.basename(full_file_path)
            RawData.objects.add_data(batch_name, file_path, file_name)
        return Response({'message': f"{count} are added. batch_name={batch_name}"})
        
   
