import json
from django.core.serializers import serialize
from django.db import models
from django.conf import settings

from commons.models import CustomUser

STATUS_OPTIONS = [
    # ready for edit/execution
    ('active', 'active'),
    # include ready tasks
    ('ready', 'ready'),
    # can't be executed because some tasks are launched
    ('locked', 'locked'),
    ('deleted', 'deleted'),
]

SEQUENCING_OPTIONS = [
    ('mrna-seq', 'mRNA-Seq'),
    ('mirna-seq', 'miRNA-Seq'),
    ('scrna-seq', 'scRNA-Seq'),
    ('other', 'other')
]

class ProjectManager(models.Manager):

    def get_last_project_id(self):
        return self.model.objects.last().project_id

    def get_next_project_id(self):
        last = self.model.objects.last()
        if last:
            next_id = str(int(last.project_id[1:]) + 1)
            return f"P{next_id.zfill(5)}"
        return f"P00001"

    def get_project_by_project_id(self, project_id:str):
        '''
        get projects by project_id
        '''
        return self.model.objects.get(project_id=project_id, status='A')

    def get_project_by_project_name(self, project_name:str):
        '''
        get projects by project_name
        '''
        return self.model.objects.filter(project_name=project_name, status='A')

    def get_projects_by_owner(self, user):
        '''
        get projects by user
        '''
        try:
            obj = self.model.objects.filter(owner=user)
            serialized_data = json.loads(serialize('json', obj))
            return serialized_data
        except Exception as e:
            print(e)
        return []
        
    def insert(self, data:dict):
        '''
        insert a new project
        project_id and owner is automatically generated
        '''
        data['project_id'] = self.model.objects.get_next_project_id()
        Project.objects.create(**data)
        return data['project_id']
   
    def delete(self, project_id):
        return Project.objects.filter(project_id=project_id).update(status='D')


class Project(models.Model):
    project_id = models.CharField(
        primary_key= True,
        max_length= 10,
        verbose_name="project ID"
    )
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        verbose_name="owner identified by user_id"
    )
    genome = models.ForeignKey(
        'rna_seq.Genome',
        on_delete=models.CASCADE,
        null=True,
        blank=True
    )
    status = models.CharField(
        max_length=10,
        default='active',
        choices=STATUS_OPTIONS,
    )
    sequencing = models.CharField(
        max_length=10,
        default='mrna-seq',
        choices=SEQUENCING_OPTIONS,
        verbose_name = "Sequencing technique",
    )
    create_time = models.DateTimeField(auto_now_add=True)
    project_name = models.CharField(
        max_length=50,
        null=True,
        blank=True,
        verbose_name='Project name',
    )
    description = models.CharField(
        max_length=526,
        null=True,
        blank=True,
    )

    objects = ProjectManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ['project_id', 'owner']

    def __str__(self):
        return self.project_id
    

    