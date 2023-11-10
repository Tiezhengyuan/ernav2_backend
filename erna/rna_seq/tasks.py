from celery import shared_task
from time import sleep
from django.core.mail import send_mail

@shared_task
def add(x, y):
    sleep(10)
    return x + y

@shared_task
def mul(x, y):
    return x * y

@shared_task
def send_email():
    subject ='subject'
    msg = "good morning"
    from_email = "tiezhengyuan@hotmail.com"
    receiver = "tiezhengyuan@hotmail.com"
    email_sent = send_mail(
      subject,
      msg,
      from_email,
      [receiver],
    )
    return email_sent