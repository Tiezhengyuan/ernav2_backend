from rest_framework import serializers
from rna_seq.models import Annotation, MolecularAnnotation, \
    AlignerIndex, Genome, Specie, RNA


class AnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = Annotation
        fields = '__all__'

class MolecularAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = MolecularAnnotation
        fields = '__all__'

class AlignerIndexSerializer(serializers.ModelSerializer):
    class Meta:
        model = AlignerIndex
        fields = '__all__'
        depth = 2

class GenomeSerializer(serializers.ModelSerializer):
    annots = AnnotationSerializer(many=True)

    class Meta:
        model = Genome
        fields = '__all__'


class SpecieSerializer(serializers.ModelSerializer):
    class Meta:
        model = Specie
        fields = '__all__'


class RNASerializer(serializers.ModelSerializer):
    class Meta:
        model = RNA
        fields = '__all__'
