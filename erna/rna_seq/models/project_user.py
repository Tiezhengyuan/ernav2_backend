
from django.db import models
from django.conf import settings

from commons.models import CustomUser
from .project import Project

class ProjectUserManager(models.Manager):

    def get_projects_by_user_name(self, user_name)->[Project]:
        projects = []
        users = CustomUser.objects.filter(user_name=user_name)
        for user in users:
            project_users = self.model.objects.filter(user=user)
            for p in project_users:
                if p.project not in projects:
                    projects.append(p.project)
        return projects

    def get_users_by_project_id(self, project_id):
        users = []
        projects = Project.objects.filter(project_id=project_id)
        for project in projects:
            project_users = self.model.objects.filter(project=project)
            for p in project_users:
                if p.user not in users:
                    users.append(p.user)
        return users


class ProjectUser(models.Model):
    project = models.ForeignKey(
        Project,
        on_delete=models.CASCADE,
        verbose_name="Executed project identified by project_id"
    )
    user = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        verbose_name="executor identified by user_id"
    )

    objects = ProjectUserManager()

    class Meta:
        app_label = 'rna_seq'
        ordering = ['project', 'user']

    def __str__(self):
        return self.project.project_id
