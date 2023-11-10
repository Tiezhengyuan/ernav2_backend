from django.test import TestCase
from ddt import ddt, data, unpack
from rna_seq.models import Project,Tool, Task
from commons.models import User

@ddt
class TestTask(TestCase):
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
        self.tool = Tool.objects.create(
            tool_name='Bowtie',
            tool_path = 'bowtie/3.4.5',
            version = '3.4.5',
        )

    def test_CRUD(self):
        res = Task.objects.create(
            task_id='T01',
            project = self.project,
            executor=self.user,
            tool=self.tool
        )
        assert str(res) == 'T01'

        # get
        res = Task.objects.get(project=self.project)
        assert res.executor.user_name == 'tester'
        res = Task.objects.get(executor=self.user)
        assert res.project.project_name == 'test1'

        # udpate
        new_project = Project.objects.create(
            project_id='P00002',
            project_name='test1',
            owner=self.user
        )
        res = Task.objects.filter(project=self.project)\
            .update(project=new_project)
        assert res == 1

        # delete
        res = Task.objects.filter(project=self.project).delete()
        assert res[0] == 0
        res = Task.objects.filter(project=new_project).delete()
        assert res[0] == 1
    
    def test_get_next_task_id(self):
        res = Task.objects.get_next_task_id('P00001')
        assert str(res) == 'T01'

        Task.objects.create(
            task_id='T02',
            project = self.project,
            executor=self.user,
            tool=self.tool
        )
        res = Task.objects.get_next_task_id('P00001')
        assert str(res) == 'T03'

    def test_add_task(self):
        res = Task.objects.add_task('P00001', 'tester', 'Bowtie')
        assert str(res) == 'T01'

    def test_get_project_tasks(self):
        Task.objects.add_task('P00001', 'tester', 'Bowtie')
        res = Task.objects.get_project_tasks('P00001')
        res = [str(i) for i in res]
        assert res == ['T01']

    def test_add_config(self):
        Task.objects.add_task('P00001', 'tester', 'Bowtie')
        res = Task.objects.add_config('T01', 'params', 'input', 'output')
