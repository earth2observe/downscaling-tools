

"""
This script makes a stand-alone 'executable' of the e2o scripts. It is tested using
Anaconda on windows 64 bit and ubuntu xenial 64 bit

supported tagets:
- normal
"""

target ='normal'

from cx_Freeze import setup, Executable, hooks


import ctypes,glob,os,shutil
import matplotlib
import scipy
import sys
import e2o_dstools



pdir = os.path.dirname(sys.executable) + "/"

if sys.platform == 'win32':
    # list comes from: c:\Anaconda\conda-meta\mkl-11.3.3-1.json
    MKL_files= [pdir + "Library/bin/cilkrts20.dll",
        pdir + "Library/bin/ifdlg100.dll",
        pdir + "Library/bin/libchkp.dll",
        pdir + "Library/bin/libicaf.dll",
        pdir + "Library/bin/libifcoremd.dll",
        pdir + "Library/bin/libifcoremdd.dll",
        pdir + "Library/bin/libifcorert.dll",
        pdir + "Library/bin/libifcorertd.dll",
        pdir + "Library/bin/libifportmd.dll",
        pdir + "Library/bin/libimalloc.dll",
        pdir + "Library/bin/libiomp5md.dll",
        pdir + "Library/bin/libiompstubs5md.dll",
        pdir + "Library/bin/libmmd.dll",
        pdir + "Library/bin/libmmdd.dll",
        pdir + "Library/bin/libmpx.dll",
        pdir + "Library/bin/liboffload.dll",
        pdir + "Library/bin/mkl_avx.dll",
        pdir + "Library/bin/mkl_avx2.dll",
        pdir + "Library/bin/mkl_avx512.dll",
        pdir + "Library/bin/mkl_core.dll",
        pdir + "Library/bin/mkl_def.dll",
        pdir + "Library/bin/mkl_intel_thread.dll",
        pdir + "Library/bin/mkl_mc.dll",
        pdir + "Library/bin/mkl_mc3.dll",
        pdir + "Library/bin/mkl_msg.dll",
        pdir + "Library/bin/mkl_rt.dll",
        pdir + "Library/bin/mkl_sequential.dll",
        pdir + "Library/bin/mkl_tbb_thread.dll",
        pdir + "Library/bin/mkl_vml_avx.dll",
        pdir + "Library/bin/mkl_vml_avx2.dll",
        pdir + "Library/bin/mkl_vml_avx512.dll",
        pdir + "Library/bin/mkl_vml_cmpt.dll",
        pdir + "Library/bin/mkl_vml_def.dll",
        pdir + "Library/bin/mkl_vml_mc.dll",
        pdir + "Library/bin/mkl_vml_mc2.dll",
        pdir + "Library/bin/mkl_vml_mc3.dll",
        pdir + "Library/bin/svml_dispmd.dll"]


data_files=[]
scipy_path = os.path.dirname(scipy.__file__)
data_files.append(scipy_path)


def load_scipy_patched(finder, module):
    """the scipy module loads items within itself in a way that causes
        problems without the entire package and a number of other subpackages
        being present."""
    finder.IncludePackage("scipy._lib")  # Changed include from scipy.lib to scipy._lib
    finder.IncludePackage("scipy.misc")

hooks.load_scipy = load_scipy_patched

def mkdatatuples(thelist,destdir="."):
    """
    input list of input files output lis list of tuples including destination
    :param list:
    :return:
    """
    ret = []
    for item in thelist:
        destfile = os.path.join(destdir,os.path.basename(item))
        ret.append((item,destfile))
    return ret

data_files.append('packages.txt')
os.system('conda list' + ">" + os.path.join('packages.txt'))
# matplolib data files


mpl =  matplotlib.get_py2exe_datafiles()

mplfiles = []
for mpldir in mpl:
    ddir = os.path.join('mpl-data',os.path.basename(mpldir[0]))
    data_files.extend(mkdatatuples(mpldir[1],destdir=ddir))

if sys.platform == 'win32':
    # MKL files
    data_files.extend(mkdatatuples(MKL_files,destdir="."))
    # pcraster dll's
    ddir = "c:/pcraster4-64/lib"
    data_files.extend(mkdatatuples(glob.glob(ddir + "/*.dll"),destdir='.'))

# GDAL data files
gdaldata = os.getenv("GDAL_DATA")
data_files.extend(mkdatatuples(glob.glob(gdaldata + "/*.*"),destdir='gdal-data'))

# Add linke and precip data file
data_files.extend(mkdatatuples(glob.glob("e2o_dstools/data/*.*"),destdir='data'))
data_files.extend(mkdatatuples(glob.glob("e2o_dstools/data/Prec/*.*"),destdir='data/Prec'))

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)
#includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','bmi','bmi.wrapper',"pcraster","osgeo.ogr"]

thename = "e2o_downscale-"+e2o_dstools.__version__+'-'+e2o_dstools.__release__+'-'+target+'-'+sys.platform+'-'+nrbits

packages = ["osgeo"]


includes = ['e2o_dstools.e2o_utils']

#  "include_msvcr": True,
options = {"includes": includes, "packages": packages,'include_files': data_files, "build_exe": thename,
            'excludes': ['collections.abc']}
base=None




executables = [
    Executable('e2o_dstools/e2o_getvar.py', base=base),
    Executable('e2o_dstools/e2o_radiation.py', base=base),
    Executable('e2o_dstools/e2o_calculateEvaporation.py', base=base)
]

setup(name='e2o_dstools',
      version=e2o_dstools.__version__+'.'+e2o_dstools.__release__,
      description='e2o_dstools',
      options={"build_exe" : options},
      executables=executables,
      )