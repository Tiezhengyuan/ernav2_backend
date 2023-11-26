from rest_framework import serializers
from rna_seq.models import Annotation, Reference, Genome, \
    Specie, NonCodingRNA


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
        fields = '__all__'


class SpecieSerializer(serializers.ModelSerializer):
    class Meta:
        model = Specie
        fields = '__all__'


class NonCodingRNASerializer(serializers.ModelSerializer):
    class Meta:
        model = NonCodingRNA
        fields = '__all__'
