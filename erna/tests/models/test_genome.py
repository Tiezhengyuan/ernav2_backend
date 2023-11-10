from django.test import TestCase, override_settings
from ddt import ddt, data, unpack
from annot.models import Specie, Genome

@ddt
class TestGenome(TestCase):
    def setUp(self):
        self.specie = Specie.objects.create(
            specie_name = 'human',
            abbreviation = 'H',
            scientific_name = 'Homo Sapiens'
        )
        self.g = Genome.objects.create(
            specie = self.specie,
            version = "G39",
            file_name = 'chr12.fa',
            metadata = 'abc',
        )

    def test_CRUD(self):
        res = Genome.objects.get(specie=self.specie)
        assert str(res) == 'human'

        res = [ i.version for i in Genome.objects.all()]
        assert res == ['G39']
        
        res = Genome.objects.filter(specie=self.specie).update(metadata='11')
        assert res == 1

        res = Genome.objects.filter(specie=self.specie).delete()
        assert res[0] == 1
        res = Genome.objects.filter(specie=self.specie).delete()
        assert res[0] == 0

    @override_settings(DATA_DIR='data')
    @data(
        ['data\\genome_dna', 'data\\genome_dna\\human\\G39\\chr12.fa']
    )
    @unpack
    def test_local_dir(self, expect_local_dir, expect_dir):
        res = self.g.local_dir
        assert res == expect_local_dir
        res = self.g.full_path
        assert res == expect_dir
    
    def test_get_file_path(self):
        res = Genome.objects.get_files_path('human')
        assert res[0].endswith('chr12.fa')
    
    def test_get_versions(self):
        res = Genome.objects.get_versions('human')
        assert res == ["G39",]

    def test_get_genome(self):
        res = Genome.objects.get_genome('human')
        assert res.version == "G39"

        res = Genome.objects.get_genome('human', 'G39')
        assert res.version == "G39"
