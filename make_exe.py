from distutils.core import setup
from bbfreeze import Freezer
import ctypes



nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)



f = Freezer("e2o_dstools-" + str(nrbits) +"-bit",includes=('numpy','scipy','scipy.special._ufuncs_cxx'))
f.addScript("e2o_dstools/e2o_radiation.py")
f.addScript("e2o_dstools/e2o_getvar.py")
f.addScript("e2o_dstools/e2o_calculateEvaporation.py")
f()    # starts the freezing process
