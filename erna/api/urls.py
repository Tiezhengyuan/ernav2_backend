

from django.urls import path,include
from rest_framework import routers

from api.views import *

router = routers.DefaultRouter(trailing_slash=True)

# app: commons
router.register('user', CustomUserViewSet, basename='user')

# app: rna_seq
router.register('project', ProjectViewSet, basename='project')
router.register('project_user', ProjectUserViewSet, basename='project_user')
router.register('task', TaskViewSet, basename='task')
router.register('task_tree', TaskTreeViewSet, basename='task_tree')
router.register('task_execution', TaskExecutionViewSet, basename='task_execution')

router.register('tool', ToolViewSet, basename='tool')
router.register('method', MethodViewSet, basename='method')
router.register('method_tool', MethodToolViewSet, basename='method_tool')
router.register('method_relation', MethodRelationViewSet, basename='method_relation')
router.register('pipeline', PipelineViewSet, basename='pipeline')

router.register('raw_data', RawDataViewSet, basename='raw_data')
router.register('sample', SampleViewSet, basename='sample')
router.register('sample_file', SampleFileViewSet, basename='sample_file')
router.register('sample_project', SampleProjectViewSet, basename="sample_project")

router.register('specie', SpecieViewSet, basename='specie')
router.register('genome', GenomeViewSet, basename='genome')
router.register('annotation', AnnotationViewSet, basename='annotation')
router.register('reference', ReferenceViewSet, basename='reference')

#celery tasks
router.register('celery_task_result', TaskResultViewSet, basename='celery_task_result')

urlpatterns = [
    path('', include(router.urls)),
]