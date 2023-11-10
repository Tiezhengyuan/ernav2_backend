from ddt import ddt, data, unpack
import json
import os
from django.test import TestCase, override_settings
from rna_seq.models import Project,Tool, Task, TaskTree
from commons.models import User

@ddt
class TestTaskTree(TestCase):
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
        self.task1 = Task.objects.add_task(
            'P00001', 'tester', 'Bowtie'
        )
        self.task2 = Task.objects.add_task(
            'P00001', 'tester', 'Bowtie'
        )
        self.task3 = Task.objects.add_task(
            'P00001', 'tester', 'Bowtie'
        )
        self.task4 = Task.objects.add_task(
            'P00001', 'tester', 'Bowtie'
        )

    @data(
        ['T02', []],
        ['T03', ['T02', 'T04']]
    )
    @unpack
    def test_cascade_deletion(self, task_id, expect):
        '''
                T01
                /
                T02
              /    \
            T03     T04
        '''
        TaskTree.objects.append_task('P00001', 'T01')
        TaskTree.objects.append_task('P00001', 'T02', 'T01')
        TaskTree.objects.append_task('P00001', 'T03', 'T02')
        TaskTree.objects.append_task('P00001', 'T04', 'T02')

        Task.objects.delete_task('P00001', task_id)
        res = TaskTree.objects.BFS('P00001', 'T01')
        res = [str(c) for t,c in res]
        assert res == expect


    def test_CRUD(self):
        res = TaskTree.objects.create(task = self.task1, \
            parent = self.task2, child = self.task3)
        assert res.task.task_id == 'T01'

        # get
        res = TaskTree.objects.get(task = self.task1)
        assert res.task.task_id == 'T01'
        res = TaskTree.objects.get(parent = self.task2)
        assert res.parent.task_id == 'T02'
        res = TaskTree.objects.get(child = self.task3)
        assert res.child.task_id == 'T03'

        # udpate
        res = TaskTree.objects.filter(task=self.task1, parent=self.task2, \
            child=self.task3).update(child=self.task3)
        assert res == 1

        # delete
        res = TaskTree.objects.filter(task=self.task1, parent=self.task2, \
            child=self.task3).delete()
        assert res[0] == 1

    @data(
        ['P00001', 'T01', ['T01']],
    )
    @unpack
    def test_get_tasks(self, project_id, task_id, expect):
        TaskTree.objects.create(task = self.task1, \
            parent = self.task2, child = self.task3)
        _, res = TaskTree.objects.get_tasks(project_id, task_id)
        res = [str(i) for i in res]
        assert res == expect

    @data(
        ['P00001', 'T01', ['T02']],
    )
    @unpack
    def test_get_parents(self, project_id, task_id, expect):
        TaskTree.objects.create(task = self.task1, \
            parent = self.task2, child = self.task3)
        _, res = TaskTree.objects.get_parents(project_id, task_id)
        res = [str(i) for i in res]
        assert res == expect

    @data(
        ['P00001', 'T01', ['T03']],
    )
    @unpack
    def test_get_children(self, project_id, task_id, expect):
        TaskTree.objects.create(task = self.task1, \
            parent = self.task2, child = self.task3)
        _, res = TaskTree.objects.get_children(project_id, task_id)
        res = [str(i) for i in res]
        assert res == expect

    def test_BFS(self):
        '''
                T01
                /
                T02
              /    \
            T03     T04
        '''
        TaskTree.objects.create(task=self.task1, child=self.task2)
        TaskTree.objects.create(task=self.task2, parent=self.task1, child=self.task3)
        TaskTree.objects.create(task=self.task2, parent=self.task1, child=self.task4)
        TaskTree.objects.create(task=self.task3, parent=self.task2)
        TaskTree.objects.create(task=self.task4, parent=self.task2)

        iters = TaskTree.objects.BFS('P00001', 'T01')
        res = [str(c) for t,c in iters]
        assert res == ['T02', 'T03', 'T04']

    def test_append_task(self):
        TaskTree.objects.append_task('P00001', 'T01')
        TaskTree.objects.append_task('P00001', 'T02', 'T01')
        TaskTree.objects.append_task('P00001', 'T03', 'T02')
        TaskTree.objects.append_task('P00001', 'T04', 'T02')
        iters = TaskTree.objects.BFS('P00001', 'T01')
        res = [str(c) for t,c in iters]
        assert res == ['T02', 'T03', 'T04']

