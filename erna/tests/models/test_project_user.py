from ddt import ddt, data, unpack
import json
import os
from django.test import TestCase, override_settings
from rna_seq.models import Project, ProjectUser
from commons.models import User

@ddt
class TestTool(TestCase):
    def setUp(self):
        self.user = User.objects.create(
            user_name="tester",
            email="a@abc.com",
            password="pass"
        )
        self.project = Project.objects.create(
            project_id='P00001',
            project_name='test1',
            owner=self.user
        )

    def test_CRUD(self):
        # create
        res = ProjectUser.objects.create(
            project = self.project,
            user=self.user,
        )
        assert str(res) == 'P00001'

        # get
        res = ProjectUser.objects.get(project=self.project)
        assert res.user.user_name == 'tester'
        res = ProjectUser.objects.get(user=self.user)
        assert res.project.project_name == 'test1'

        # udpate
        new_project = Project.objects.create(
            project_id='P00002',
            project_name='test1',
            owner=self.user
        )
        res = ProjectUser.objects.filter(project=self.project)\
            .update(project=new_project)
        assert res == 1

        # delete
        res = ProjectUser.objects.filter(project=self.project).delete()
        assert res[0] == 0
        res = ProjectUser.objects.filter(project=new_project).delete()
        assert res[0] == 1

    @data(
        ['tester', ['P00001',]],
        ['wrong', []],
    )
    @unpack
    def test_get_projects_by_user_name(self, user_name, expect):
        ProjectUser.objects.create(project = self.project, user=self.user)
        res = ProjectUser.objects.get_projects_by_user_name(user_name)
        res = [str(i) for i in res]
        assert res == expect

    @data(
        ['P00001', ['a@abc.com',]],
        ['wrong', []],
    )
    @unpack
    def test_get_users_by_project_id(self, project_id, expect):
        ProjectUser.objects.create(project = self.project, user=self.user)
        res = ProjectUser.objects.get_users_by_project_id(project_id)
        res = [str(i) for i in res]
        assert res == expect
