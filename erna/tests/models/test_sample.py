from django.test import TestCase
from ddt import ddt, data, unpack
from commons.models import User
from sample.models import Sample

@ddt
class TestSamples(TestCase):
    def setUp(self):
        self.user = User.objects.create(
            user_name="tester",
            email="a@abc.com",
            password="pass"
        )
        Sample.objects.create(
            study_name='study1',
            sample_name='sample1',
            creator=self.user
        )

    def test_CRUD(self):

        res = Sample.objects.get(study_name='study1')
        assert str(res) == 'study1_sample1'
        assert res.study_name == 'study1'
        assert res.sample_name == 'sample1'
        assert res.metadata == None

        res = [ i.sample_name for i in Sample.objects.all()]
        assert res == ['sample1']
        
        res = Sample.objects.filter(study_name='study1')\
            .update(sample_name='sample2')
        assert res == 1

        res = Sample.objects.filter(study_name='study1').delete()
        assert res[0] == 1
        res = Sample.objects.filter(study_name='study1').delete()
        assert res[0] == 0

    def test_unique(self):
        res = Sample.objects.get(study_name='study1',
            sample_name='sample1')
        assert str(res) == 'study1_sample1'
        assert res.metadata is None

        res = Sample.objects.create(
            study_name='study1',
            sample_name='sample2',
            creator=self.user
        )
        assert str(res) == 'study1_sample2'

        res = Sample.objects.create(
            study_name='study2',
            sample_name='sample1',
            creator=self.user
        )
        assert str(res) == 'study2_sample1'

        res = Sample.objects.create(
            study_name='study1',
            sample_name='sample1',
            creator=self.user,
            metadata = 'data'
        )
        assert res.metadata == 'data'
    
    @data(
        ['study1', True],
        ['study2', False],
    )
    @unpack
    def test_study_exists(self, input, expect):
        res = Sample.objects.study_exists(input)
        assert res == expect

    @data(
        ['sample1', True],
        ['sample2', False],
    )
    @unpack
    def test_sample_exists(self, input, expect):
        res = Sample.objects.sample_exists('study1', input)
        assert res == expect


    @data(
        ['study1', ['sample1']],
        ['study2', []],
    )
    @unpack
    def test_get_sample_names_by_study(self, input, expect):
        res = Sample.objects.get_sample_names_by_study(input)
        assert res == expect

    @data(
        ['tester', {'study1': ['sample1']}],
    )
    @unpack
    def test_get_sample_names_by_user(self, input, expect):
        res = Sample.objects.get_sample_names_by_user(input)
        assert res == expect

    def test_load_samples(self):
        samples = [
            {'study_name':'study1', 'sample_name': 'sample2'},
            {'study_name':'study1', 'sample_name': 'sample3', 'medata':{}},
            {'study_name':'study1', 'sample_name': 'sample4', 'medata':{'a':4}},
        ]
        res = Sample.objects.load_samples('tester', samples)
        assert res == 3
    
    @data(
        ['sample1', 'sample2', 'sample2'],
        ['sample1', 'sample1', None],
        ['wrong', 'sample1', None],
    )
    @unpack
    def test_update_sample_name(self, old_name, new_name, expect):
        res = Sample.objects.update_sample_name('study1', old_name, new_name)
        if expect:
            assert res.sample_name == expect
        else:
            assert res == expect

    @data(
        ['sample1', 1],
        ['wrong', 0],
    )
    @unpack
    def test_delete_sample(self, sample_name, expect):
        res = Sample.objects.delete_sample('study1', sample_name)
        assert res[0] == expect

    @data(
        ['study1', 2],
        ['wrong', 0],
    )
    @unpack
    def test_delete_study_samples(self, study_name, expect):
        Sample.objects.create(
            study_name='study1',
            sample_name='sample2',
            creator=self.user
        )
        res = Sample.objects.delete_study_samples(study_name)
        assert res[0] == expect

    @data(
        ['study1', [None, 'abc']],
        ['wrong', []],
    )
    @unpack
    def test_export_study(self, study_name, expect):
        Sample.objects.create(
            study_name='study1',
            sample_name='sample2',
            creator=self.user,
            metadata = 'abc'
        )
        res = Sample.objects.export_study(study_name)
        assert [i['metadata'] for i in res.values()] == expect

    def test_export_study(self):
        res = Sample.objects.get_study_names()
        assert res == {'study1',}
        
        Sample.objects.create(
            study_name='study1',
            sample_name='sample2',
            creator=self.user,
            metadata = 'abc'
        )
        res = Sample.objects.get_study_names()
        assert res == {'study1',}

        Sample.objects.create(
            study_name='study2',
            sample_name='sample2',
            creator=self.user,
            metadata = 'abc'
        )
        res = Sample.objects.get_study_names()
        assert res == {'study1', 'study2'}
