def getBase():
    import os
    import sys
    # needed to get Base
    sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
getBase()    
from BulletWrap import *

