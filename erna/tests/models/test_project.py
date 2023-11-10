from ddt import ddt, data, unpack
from django.test import TestCase
from rna_seq.models import  Project
from commons.models import User

@ddt
class TestProject(TestCase):
    
    def setUp(self):
        self.user = User.objects.create(
            user_name="tester",
            email="a@abc.com",
            password="pass"
        )
        User.objects.create(
            user_name="user",
            email="b@abc.com",
            password="ative"
        )


    @data(
        ['test1', 'test1', 'P00001', 'tester'],
    )
    @unpack
    def test_project_model(self, project_name, 
            expect_name, expect_id, expect_user):
        Project.objects.create(
            project_id='P00001',
            project_name='test1',
            owner=self.user
        )
        res = Project.objects.get(project_name=project_name)
        assert res.project_name == expect_name
        assert str(res.owner.user_name) == expect_user
        assert res.project_id == expect_id

    def test_last_project_id(self):
        Project.objects.create(project_id='P00001', \
            project_name='test1', owner=self.user)
        res = Project.objects.get_last_project_id()
        assert res == 'P00001'

    def test_get_next_project_id(self):
        res = Project.objects.get_next_project_id()
        assert res == 'P00001'

        Project.objects.create(project_id='P00001', \
            project_name='test1', owner=self.user)
        res = Project.objects.get_next_project_id()
        assert res == 'P00002'

        Project.objects.create(project_id='P00002', \
            project_name='test2', owner=self.user)
        res = Project.objects.get_next_project_id()
        assert res == 'P00003'

    def test_insert(self):
        data = {'project_name': 'test1', 'owner':self.user}
        res = Project.objects.insert(data)
        assert res == 'P00001'
        data = {'project_name': 'test2', 'owner':self.user}
        res = Project.objects.insert(data)
        assert res == 'P00002'

    @data(
        ['tester', ['P00001', 'P00002']], 
        ['wrong', None],
        ['user',[]],
    )
    @unpack
    def test_get_project_by_owner_name(self, owner_name, expect):
        Project.objects.create(project_id='P00001', \
            project_name='test1', owner=self.user)
        Project.objects.create(project_id='P00002', \
            project_name='test2', owner=self.user)
        res = Project.objects.get_project_by_owner_name(owner_name)
        if expect is not None:
            res = [i.project_id for i in res]
        assert res == expect

    @data(
        ['P00001', True], 
        ['P00002', False], 
    )
    @unpack
    def test_delete(self, project_id, expect):
        Project.objects.create(project_id='P00001', \
            project_name='test1', owner=self.user)
        res = Project.objects.delete(project_id)
        assert res == expect
    
    @data(
        ['P00001', 'P00001'], 
    )
    @unpack
    def test_get_project_by_project_id(self, project_id, expect):
        Project.objects.create(
            project_id='P00001',
            project_name='test1',
            owner=self.user
        )
        res = Project.objects.get_project_by_project_id(project_id)
        assert str(res) == expect

    @data(
        ['test1', ['test1']], 
        ['wrong', []],
    )
    @unpack
    def test_get_project_by_project_name(self, project_name, expect):
        Project.objects.create(
            project_id='P00001',
            project_name='test1',
            owner=self.user
        )
        res = Project.objects.get_project_by_project_name(project_name)
        # print(f"##{res}")
        res = [i.project_name for i in res]
        assert res == expect
