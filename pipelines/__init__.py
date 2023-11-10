import os
import sys
DIR_PROJECT = os.path.dirname(os.path.dirname(__file__))
DIR_PIPELINES = os.path.join(DIR_PROJECT, 'pipelines')

sys.path.append(DIR_PIPELINES)