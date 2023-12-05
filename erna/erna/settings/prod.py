"""
Django settings for erna project in PROD.
"""
import os
from pathlib import Path
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv('.env.prod'))
from celery.schedules import crontab
import json

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = os.path.dirname(BASE_DIR)

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False

# Django settings for projects

ALLOWED_HOSTS = ['*']


CORS_ORIGIN_ALLOW_ALL = False
CORS_ORIGIN_WHITELIST = (
    'http://localhost:8080',
)



TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            os.path.join(BASE_DIR, 'templates'),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]


# Database
# https://docs.djangoproject.com/en/4.2/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    }
}

print(BASE_DIR / 'db.sqlite3')


# Static files (CSS, JavaScript, Images)
# STATICFILES_DIRS = [
#     os.path.join(BASE_DIR, 'static'),
# ]
STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')



# other related path
PIPELINES_DIR = os.environ['PIPELINES_DIR']

# main entrance for launch bioinformatics pipeline
PIPELINE_ERNA = os.path.join(BASE_DIR, 'erna_app.py')

# third-party bioinformatics tools
EXTERNALS_DIR = os.environ['EXTERNALS_DIR']

# store data uploaded by user
DATA_DIR = os.environ['DATA_DIR']

# raw data namely fastq
RAW_DATA_DIRS = os.environ['RAW_DATA_DIR'].split(' ')

# analytic results
RESULTS_DIR = os.environ['RESULTS_DIR']

# reference namely genome DNA in fa format
REFERENCES_DIR = os.environ['REFERENCES_DIR']


#celery settings
CELERY_BROKER_URL="redis://127.0.0.1:6379"
CELERY_RESULT_BACKEND = 'django-db'
CELERY_RESULT_EXTENDED = True
CELERY_ACCEPT_CONTENT = ['application/json']
CELERY_RESULT_SERIALIZER = 'json'
CELERY_TASK_SERIALIZER = 'json'
CELERY_TIMEZONE = "America/New_York"
CELERY_ALWAYS_EAGER = True
CELERY_BEAT_SCHEDULER = 'django_celery_beat.schedulers:DatabaseScheduler'
CELERY_BEAT_SCHEDULE = {
    "execute_task": {
        "task": "celery_tasks.tasks.execute_task",
        "schedule": crontab(minute="*/1"),
    },
}
