from rest_framework import serializers
from rna_seq.models import RawData, Sample, SampleFile, SampleProject


class RawDataSerializer(serializers.ModelSerializer):
    class Meta:
        model = RawData
        fields = '__all__'

class SampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sample
        fields = '__all__'

class SampleProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleProject
        fields = '__all__'


class SampleFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleFile
        fields = '__all__'


