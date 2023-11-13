'''
scheduled tasks:
search table Task and detect tasks with 'is_ready'=True
'''

from rna_seq.models import Project, Task

class ScheduleTasks:
    def __init__(self):
        pass
    
    def run_task(self):
        ready_tasks = Task.objects.filter(is_ready=True)
        for task in ready_tasks:
            print(task)