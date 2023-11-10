from django.urls import path
from . import views

app_name = 'rna_seq'

urlpatterns = [
  path('index/', views.index, name='index'),
  path('async/', views.async_task, name='async'),
  path('async_minus/', views.async_minus, name='async_minus'),
  path('async_email/', views.async_email, name='async_email'),
]