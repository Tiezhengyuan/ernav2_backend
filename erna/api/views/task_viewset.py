import json
from django.shortcuts import render
from django.core import serializers
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework import viewsets, permissions

from rna_seq.models import Task
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
  def next_task_id(self, request):
      project_id = self.request.query_params.get('project_id')
      if project_id is not None:
        res = Task.objects.next_task_id(project_id)
        return Response(res)
      return Response({'error': f'project_id={project_id} is missing.'})

  @action(detail=False, methods=['get'])
  def project_tasks(self, request):
    project_id = self.request.query_params.get('project_id')
    res = []
    if project_id is not None:
      for task in Task.objects.filter(project=project_id):
        item = {
          'project_id': task.project.project_id,
          'task_id': task.task_id,
          'params': task.params,
          'method_tool_id': task.method_tool.id if task.method_tool else None,
          'method_name': task.method_tool.method.method_name if task.method_tool else None,
          'status': task.task_execution.status if task.task_execution else None,
        }
        res.append(item)
    return Response(res)
         
