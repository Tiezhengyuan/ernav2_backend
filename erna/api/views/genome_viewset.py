
from rest_framework.decorators import action
from rest_framework import viewsets, permissions
from rest_framework.response import Response

from annot.models import Genome
from api.serializers import GenomeSerializer


class GenomeViewSet(viewsets.ModelViewSet):
    serializer_class = GenomeSerializer
    permission_classes = [permissions.IsAuthenticated]

    def get_queryset(self):
        names = ['specie', 'data_source', 'version', 'is_ready']
        params = dict([(k, self.request.query_params.get(k)) \
            for k in names if self.request.query_params.get(k)])
        if params:
            return Genome.objects.filter(**params)
        return Genome.objects.all()


    @action(detail=False, methods=['get'])
    def data_sources(self, request):
        res = Genome.objects.values_list('data_source', flat=True)
        return Response({'names': list(set(res))})
    
    @action(detail=False, methods=['post', 'update'])
    def load_genomes(self, request):
        res = {'created': 0, 'skipped': 0}
        genomes = request.data
        if genomes:
            for genome in genomes:
                print(genome)
                obj = Genome.objects.load_genome(genome)
                if obj:
                    res['created'] += 1
                else:
                    res['skipped'] += 1
        return Response(res)

    @action(detail=False, methods=['delete'])
    def delete_all(self, request):
        qs = Genome.objects.all()
        count = qs.count()
        qs.delete()
        return Response({'deleted': count})
