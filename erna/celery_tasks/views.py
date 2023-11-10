import json
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from django_celery_beat.models import IntervalSchedule, PeriodicTask

from .tasks import *

def DownloadGenomeView(request):
  data_source = request.GET.get('data_source', '')
  specie = request.GET.get('specie')
  version = request.GET.get('version')
  task_id = download_genome.delay(data_source, specie, version)
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)

def ScanRawDataView(request):
  task_id = scan_raw_data.delay()
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)

def RefreshRawDataView(request):
  task_id = refresh_raw_data.delay()
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)

def ParseSampleDataView(request):
  study_name = request.GET.get('study_name', '')
  prefix = request.GET.get('prefix')
  postfix = request.GET.get('postfix')
  task_id = parse_sample_data.delay(study_name, prefix, postfix)
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)

def ResetSampleView(request):
  task_id = reset_sample.delay()
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)

def TrimAdapterView(request):
  params = {
    'adapter': request.GET.get('adapter'),
    'min_match': request.GET.get('min_match'),
    'max_err': request.GET.get('max_err'),
    'input': request.GET.get('input'),
    'output': request.GET.get('output'),
  }
  task_id = trim_adapter.delay(params)
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)

def BuildIndexView(request):
  params = {
    'data_source': request.GET.get('data_source'),
    'specie': request.GET.get('specie'),
    'version': request.GET.get('version'),
    'aligner': request.GET.get('aligner'),
  }
  task_id = build_index.delay(**params)
  res = {
    'task_id': [str(task_id),],
  }
  return JsonResponse(res, safe=False)



'''
for debugging
'''
def test_async(request):
  '''
  celery task
  '''
  task_id = minus.delay(2, 3)
  return HttpResponse(task_id)

def test_schedule(request):
  '''
  scheduled celery task
  '''
  interval, err = IntervalSchedule.objects.get_or_create(
    every = 30,
    period = IntervalSchedule.SECONDS
  )
  task_id = PeriodicTask.objects.create(
    interval = interval,
    name="test_schedule",
    task = 'celery_tasks.tasks.minus',
    args = json.dumps([1, 2])
  )
  return HttpResponse(task_id)