from django.db import models


class BaseModel(models.Model):
    
    @staticmethod
    def get_last_object():
        return self.objects.last()