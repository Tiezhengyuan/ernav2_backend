from django.test import TestCase
from ddt import ddt, data, unpack
import json
from rna_seq.models import Project,Tool, Task, TaskExecution
from commons.models import User

@ddt
class TestTaskExecution(TestCase):
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
        self.task = Task.objects.add_task(
            'P00001', 'tester', 'Bowtie'
        )

    def test_CRUD(self):
        res = TaskExecution.objects.create(task=self.task, command='abc')
        assert res.task.task_id == 'T01'
        assert res.id == 1
        

        # get
        res = TaskExecution.objects.get(id=1)
        assert res.task.task_id == 'T01'
        assert res.id == 1
        assert res.status == 'R'

        # udpate
        res = TaskExecution.objects.filter(id=1).update(status='E')
        assert res == 1

        # delete
        res = TaskExecution.objects.filter(id=1).delete()
        assert res[0] == 1
    
    def test_execution_cycle(self):
        res = TaskExecution.objects.initialize(project_id='P00001', \
            task_id='T01', command='abc')
        assert res.status == 'R'

        TaskExecution.objects.run(id=1)
        res = TaskExecution.objects.get(id=1)
        assert res.status == 'E'

        TaskExecution.objects.pause(id=1)
        res = TaskExecution.objects.get(id=1)
        assert res.status == 'P'

        TaskExecution.objects.stop(id=1)
        res = TaskExecution.objects.get(id=1)
        assert res.status == 'D'

        TaskExecution.objects.fail(id=1)
        res = TaskExecution.objects.get(id=1)
        assert res.status == 'F'
    
    def test_get_lastest_task(self):
        TaskExecution.objects.initialize(project_id='P00001', \
            task_id='T01', command='abc')
        TaskExecution.objects.initialize(project_id='P00001', \
            task_id='T01', command='abc')

        res = TaskExecution.objects.get_lastest_task('P00001', 'T01')
        assert res.id == 2
    
    def test_get_project_status(self):
        TaskExecution.objects.initialize('P00001', 'T01', 'abc')

        Task.objects.add_task('P00001', 'tester', 'Bowtie')
        TaskExecution.objects.initialize('P00001', 'T02', 'abc')
        TaskExecution.objects.run(2)

        Task.objects.add_task('P00001', 'tester', 'Bowtie')
        TaskExecution.objects.initialize('P00001', 'T03', 'abc')
        TaskExecution.objects.pause(3)

        Task.objects.add_task('P00001', 'tester', 'Bowtie')
        TaskExecution.objects.initialize('P00001', 'T04', 'abc')
        TaskExecution.objects.stop(4)

        Task.objects.add_task('P00001', 'tester', 'Bowtie')
        TaskExecution.objects.initialize('P00001', 'T05', 'abc')
        TaskExecution.objects.stop(5)

        res = TaskExecution.objects.get_project_status('P00001')
        assert res == {1: 'R', 2: 'E', 3: 'P', 4: 'D', 5: 'D'}
    
    def test_get_command(self):
        params = {'a': 4}
        TaskExecution.objects.initialize(project_id='P00001', \
            task_id='T01', command=json.dumps(params))
        res = TaskExecution.objects.get_command(1)
        assert json.loads(res) == params
