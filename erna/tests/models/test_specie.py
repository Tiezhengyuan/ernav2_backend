from django.test import TestCase
from annot.models import Specie


class TestSpecie(TestCase):
    def setUp(self):
        Specie.objects.create(
            specie_name = 'human',
            abbreviation = 'H',
            scientific_name = 'Homo Sapiens'
        )

    def test_CRUD(self):
        res = Specie.objects.get(specie_name = 'human')
        assert res.specie_name == 'human'

        res = [ i.specie_name for i in Specie.objects.all()]
        assert res == ['human',]
        
        res = Specie.objects.filter(specie_name = 'human').update(abbreviation = 'HH')
        assert res == 1

        res = Specie.objects.filter(specie_name = 'human').delete()
        assert res[0] == 1
        res = Specie.objects.filter(specie_name = 'human').delete()
        assert res[0] == 0
