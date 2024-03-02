from typing import Iterable
from django.db import models

from .method_relation import MethodRelation
from .task import Task
from .project import Project

class TaskTreeManager(models.Manager):

    def get_tasks(self, project_id:str, task_id:str):
        task = Task.objects.get(project_id=project_id, task_id=task_id)
        return task, self.model.objects.filter(parent=task)

    def get_parents(self, project_id:str, task_id:str):
        '''
        one task may have multiple parents
        '''
        task = Task.objects.get(project_id=project_id, task_id=task_id)
        parents = [task_tree.parent for task_tree in self.filter(child=task)]
        return task, parents

    def get_children(self, project_id:str, task_id:str):
        '''
        one task may have multiple children
        '''
        task = Task.objects.get(project_id=project_id, task_id=task_id)
        children = [task_tree.child for task_tree in self.filter(parent=task)]
        return task, children
    
    def task_tree(self, project_id:str):
        res = {}
        # possible paretns give a method name
        method_parents = MethodRelation.objects.get_parents()
        # all tasks given a project
        project_tasks = Task.objects.filter(project=project_id)
        for task in project_tasks:
            parents = TaskTree.objects.filter(child=task)
            parent_ids = [i.parent.task_id for i in parents]
            if task.method_tool and task.method_tool.method:
                possible_parents = method_parents[task.method_tool.method.method_name]
                items = []
                for task2 in project_tasks:
                    if task2.method_tool:
                        method_name = task2.method_tool.method.method_name
                        if task2 != task and method_name in possible_parents:
                            item = {
                                'value': task2.task_id,
                                'text': task2.task_id,
                                'check': True if task2.task_id in parent_ids else False,
                            }
                            items.append(item)
                    else:
                        print('wrong task2: ', task2.task_id)
            else:
                print('wrong task', task.task_id)
            res[task.task_id] = items
        return res

    def BFS(self, project_id:str, task_id:str)->Iterable:
        '''
        task_id is root task. Given task_id,
        search all children using breadth-first search 
        '''
        i, pool = 0, [task_id]
        while i < len(pool):
            task_id = pool[i]
            task, children = self.model.objects.get_children(project_id, task_id)
            for child in children:
                if child and child not in pool:
                    pool.append(child)
                    yield task, child
            i += 1

    def load_tasks_tree(self, project_id:str, tasks_data:list):
        '''
        example: tasks_data = [('T01', 'T03'), ...],
        '''
        if not project_id:
            return []

        res = []
        project = Project.objects.get(project_id=project_id)
        for parent_id, child_id in tasks_data:
            parent = Task.objects.get(project=project, task_id=parent_id)
            child = Task.objects.get(project=project, task_id=child_id)
            obj = self.update_or_create(parent=parent, child=child)
            res.append(obj)
        return res


class TaskTree(models.Model):
    '''
    task, parent and child must belong to same project
    '''
    parent = models.ForeignKey(
        Task,
        on_delete=models.CASCADE,
        related_name='parent_tasks',
        verbose_name='Parent Task',
    )
    child = models.ForeignKey(
        Task,
        on_delete=models.CASCADE,
        related_name='child_tasks',
        verbose_name='Child task',
        blank=True,
        null=True,
    )

    objects = TaskTreeManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ('parent', 'child')
   
