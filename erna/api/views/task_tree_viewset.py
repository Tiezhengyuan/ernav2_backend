from django.shortcuts import render
from rest_framework import viewsets, permissions
from rest_framework.response import Response
from rest_framework.decorators import action

from rna_seq.models import Task, TaskTree, MethodRelation
from api.serializers import *

class TaskTreeViewSet(viewsets.ModelViewSet):
    queryset = TaskTree.objects.all()
    serializer_class = TaskTreeSerializer
    permission_classes = [permissions.IsAuthenticated]



    @action(detail=False, methods=['get'])
    def task_parents(self, request):
        project_id = self.request.query_params.get('project_id')
        res = {}
        if project_id is not None:
            # possible paretns give a method name
            method_parents = MethodRelation.objects.get_parents()
            # all tasks given a project
            project_tasks = Task.objects.filter(project=project_id)
            for task in project_tasks:
                parents = TaskTree.objects.filter(child=task)
                parent_ids = [i.parent.task_id for i in parents]
                possible_parents = method_parents[task.method_tool.method.method_name]
                items = []
                for task2 in project_tasks:
                    if task2 != task and task2.method_tool.method.method_name in possible_parents:
                        item = {
                            'value': task2.task_id,
                            'text': task2.task_id,
                            'check': True if task2.task_id in parent_ids else False,
                        }
                        items.append(item)
                res[task.task_id] = items
        return Response(res)

    @action(detail=False, methods=['post'])
    def update_task_pairs(self, request):
        project_id = request.data.get('project_id')
        child_task_id = request.data.get('child')
        parent_ids = request.data.get('parents')
        res = []
        if project_id and child_task_id:
            child_task, current_parents = TaskTree.objects.get_parents(project_id, child_task_id)
            for current_parent in current_parents:
                if current_parent.parent.task_id in parent_ids:
                    parent_ids.remove(current_parent.parent.task_id)
                else:
                    current_parent.delete()
            for parent_task_id in parent_ids:
                parent_task = Task.objects.get(project_id=project_id, task_id=parent_task_id)
                obj = TaskTree.objects.create(parent=parent_task, child=child_task)
                print(obj)
            _, new_parents = TaskTree.objects.get_parents(project_id, child_task_id)
            res = [i.parent.task_id for i in new_parents]
            print(res)
        return Response(res)

