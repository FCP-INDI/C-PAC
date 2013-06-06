from behave import *
from hamcrest import assert_that, is_not, greater_than

import numpy as np
import nibabel as nib
import rpy2.robjects as robjects
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr
robjects.conversion.py2ri = numpy2ri

from os import path as op
import sys
curfile  = op.abspath(__file__)
testpath = op.dirname(op.dirname(op.dirname(curfile)))
rpath    = op.join(testpath, "R")
pypath   = op.dirname(testpath)
sys.path.append(pypath)

from cwas import *
from utils import *
