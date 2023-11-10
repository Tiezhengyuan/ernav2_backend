from rest_framework import serializers
from annot.models import *


class AnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = Annotation
        fields = '__all__'

class ReferenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reference
        fields = '__all__'
        depth = 2

class GenomeSerializer(serializers.ModelSerializer):
    annots = AnnotationSerializer(many=True)
    references = ReferenceSerializer(many=True)

    class Meta:
        model = Genome
        fields = ('data_source', 'version', 'specie', 'is_ready', \
            'local_path', 'metadata', 'ftp_path', \
            'annots', 'references')


class SpecieSerializer(serializers.ModelSerializer):
    class Meta:
        model = Specie
        fields = '__all__'

