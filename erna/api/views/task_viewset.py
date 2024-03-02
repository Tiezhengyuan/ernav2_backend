import json
from django.shortcuts import render
from django.core import serializers
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework import viewsets, permissions

from rna_seq.models import Task, TaskTree
from api.serializers import TaskSerializer


class TaskViewSet(viewsets.ModelViewSet):
  serializer_class = TaskSerializer
  permission_classes = [permissions.IsAuthenticated]

  def get_queryset(self):
    config = {}
    for name in ['project_id', 'task_id']:
      val = self.request.query_params.get(name)
      if val is not None:
        config[name] = val
    if config:
      return Task.objects.filter(**config)
    return Task.objects.all()

  @action(detail=False, methods=['post'])
  def load_tasks(self, request):
      project_id = request.data.get('project_id')
      tasks = request.data.get('tasks')
      if project_id and tasks:
        res = Task.objects.load_tasks(project_id, tasks)
        return Response({'count': len(res)})
      return Response({'error': 'project_id is missing.'})
  
  @action(detail=False, methods=['delete'])
  def delete_tasks(self, request):
    project_id = request.GET.get('project_id')
    task_id = request.GET.get('task_id')
    res = Task.objects.delete_tasks(project_id, task_id)
    return Response(res)
  
  @action(detail=False, methods=['get'])
  def front_project_tasks(self, request):
    '''
    response is consumeb by NewTask.vue
    '''
    project_id = self.request.query_params.get('project_id')
    if project_id is None:
      return Response({'error': f'project_id={project_id} is missing.'})

    res = {
      'task_tree': TaskTree.objects.task_tree(project_id),
      'tasks': Task.objects.project_tasks(project_id),
      'new_task_id': Task.objects.next_task_id(project_id),
    }
    return Response(res)
         
