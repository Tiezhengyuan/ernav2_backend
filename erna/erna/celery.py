import os
import sys
from celery import Celery

ernav2_dir = os.path.dirname(
  os.path.dirname(os.path.dirname(__file__))
)
sys.path.append(ernav2_dir)
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'erna.settings')

app = Celery('erna');
app.config_from_object('django.conf:settings', namespace='CELERY')
app.autodiscover_tasks()

@app.task(bind=True)
def debug_task(self):
  print(f'Reqyest: {self.request}')