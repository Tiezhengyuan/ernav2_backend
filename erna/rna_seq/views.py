from django.http import HttpResponse
from django.shortcuts import render

# Create your views here.
from .tasks import add, send_email
from celery_tasks.tasks import minus

def index(request):
  res = add(5,3)
  return HttpResponse(res)

def async_task(request):
  res_id = add.delay(23,4)
  return HttpResponse(res_id)

def async_email(request):
  # instruct celery to execute in the background
  res = send_email.delay()
  return HttpResponse(res)

def async_minus(request):
  res = minus.delay(3,4)
  return HttpResponse(res)