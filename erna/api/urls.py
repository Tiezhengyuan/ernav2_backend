

from django.urls import path,include
from rest_framework import routers

from api.views import *

router = routers.DefaultRouter(trailing_slash=True)

router.register('annotation', AnnotationViewSet, basename='annotation')
router.register('celery_task_result', TaskResultViewSet, basename='celery_task_result')
router.register('execution_tree', ExecutionTreeViewSet, basename='execution_tree')
router.register('genome', GenomeViewSet, basename='genome')

router.register('method', MethodViewSet, basename='method')
router.register('method_relation', MethodRelationViewSet, basename='method_relation')
router.register('method_tool', MethodToolViewSet, basename='method_tool')

router.register('pipeline', PipelineViewSet, basename='pipeline')
router.register('project', ProjectViewSet, basename='project')
router.register('project_user', ProjectUserViewSet, basename='project_user')

router.register('raw_data', RawDataViewSet, basename='raw_data')
router.register('reference', ReferenceViewSet, basename='reference')
router.register('rna', RNAViewSet, basename='rna')

router.register('sample', SampleViewSet, basename='sample')
router.register('sample_file', SampleFileViewSet, basename='sample_file')
router.register('sample_project', SampleProjectViewSet, basename="sample_project")
router.register('specie', SpecieViewSet, basename='specie')

router.register('task', TaskViewSet, basename='task')
router.register('task_annotation', TaskAnnotationViewSet, basename='task_annotation')
router.register('task_execution', TaskExecutionViewSet, basename='task_execution')
router.register('task_tree', TaskTreeViewSet, basename='task_tree')
router.register('tool', ToolViewSet, basename='tool')

router.register('user', CustomUserViewSet, basename='user')



urlpatterns = [
    path('', include(router.urls)),
]