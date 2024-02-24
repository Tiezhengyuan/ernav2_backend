from rest_framework.decorators import action
from rest_framework import viewsets, permissions, status, serializers
from rest_framework.response import Response

from rna_seq.models import Specie
from api.serializers import SpecieSerializer


class SpecieViewSet(viewsets.ModelViewSet):
    serializer_class = SpecieSerializer
    permission_classes = [permissions.IsAuthenticated,]

    def get_queryset(self):
        group = self.request.query_params.get('group')
        organism_name = self.request.query_params.get('organism_name')
        if group:
            if organism_name:
                return Specie.objects.filter(
                    group=group, organism=organism_name)
            return Specie.objects.filter(group=group)
        else:
            if organism_name:
                return Specie.objects.filter(organism_name=organism_name)
        return Specie.objects.all()

    # TODO debugging in the future
    # @action(detail=False, methods=['get'])
    # def reg_search(self, request):
    #     organism_name = self.request.query_params.get('organism_name')
    #     if organism_name:
    #         qs = Specie.objects.filter(
    #             organism_name__contains=organism_name)
    #         data = SpecieSerializer(qs).data
    #         print(data)
    #         return Response(data, status=status.HTTP_200_OK)
    #     return Response({})

    @action(detail=False, methods=['get'])
    def count(self, request):
        count = Specie.objects.count()
        return Response({'count': count})

    @action(detail=False, methods=['get'])
    def group_names(self, request):
        groups = Specie.objects.values_list('group', flat=True)
        res = [{'text': i.replace('_', ' '), 'value': i} for i in list(set(groups)) if i]
        return Response(res)

    @action(detail=False, methods=['get'])
    def count_by_group(self, request):
        groups = {}
        res = Specie.objects.values_list('group', flat=True).distinct()
        for group in list(set(res)):
            groups[group] = Specie.objects.filter(group=group).count()
        return Response(groups)
