from django.urls import path
from .views import *

urlpatterns = [
  path('execute_tasks/', ExecuteTasksView, name='execute_tasks'),
  path('download_genome/', DownloadGenomeView, name='download_genome'),
  path('scan_raw_data/', ScanRawDataView, name='scan_raw_data'),
  path('refresh_raw_data/', RefreshRawDataView, name='refresh_raw_data'),
  path('parse_sample_data/', ParseSampleDataView, name='parse_sample_data'),
  path('reset_sample/', ResetSampleView, name='reset_sample'),
  path('trim_adapter', TrimAdapterView, name='trim_adapter'),
  path('build_index', BuildIndexView, name='build_index'),

  # for debugging
  path('test_async/', test_async, name='test_async'),
  path('test_schedule/', test_schedule, name='test_schedule'),
]