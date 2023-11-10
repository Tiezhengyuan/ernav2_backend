# ernav2_backend
Back-end of eRNAv2 including Django web end and bioinformatics pipelines

# Development

## Run local version of eRNA_v2.
The local version doesn't provide user authentication. And the default databse is sqllite3 rather than MySQL. The app may not be remotely accessed.

### run back end (Django) locally

Django project known as "erna" could be runned locally.
```
cd erna
source venv/bin/activate
pip install -r requiements.txt
python manage.py runserver
```
Once Django server is running:
RestFull APIs are available at http://localhost:8000/api/.
Django Administration is accessible at http://localhost:8000/admin/.

Start Celery in another terminal showed as the below.
It is recommended to test a demo celery task using http://localhost:8000/celery_tasks/test_async/ in broswer. Task ID would reveal asynchronous tasks execution in eRNAv2 is working.
```
cd erna
source venv/bin/activate
celery -A erna worker --beat --scheduler django -l info --pool=solo
```
