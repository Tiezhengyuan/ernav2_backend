import os
from .base import *


if os.environ.get('mode') == 'prod':
    from .prod import *
else:
    from .dev import *