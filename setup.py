import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from setuptools import find_packages
here = os.path.abspath(os.path.dirname(__file__))

README = open(os.path.join(here, 'README.rst')).read()

requires = [
    'netCDF4',
    'numpy',
    'matplotlib',
    'pcraster', 'osgeo','scipy']

datadir = os.path.join('e2o_dstools','data')
datafiles = [(d, [os.path.join(d,f) for f in files])
    for d, folders, files in os.walk(datadir)]



# Source dist
setup(name='e2o_dstoools',
      version= "2017.2",
      packages=['e2o_dstools'],
      package_dir={'e2o_dstools': 'e2o_dstools'},
      author='J. Schellekens/F. Sperna Weiland',
      author_email='jaap.schellekens@deltares.nl',
      url='http://www.earth2observe.eu',
      license = "GPL",
      scripts=['e2o_dstools/e2o_getvar.py','e2o_dstools/e2o_calculateEvaporation.py',
               'e2o_dstools/e2o_radiation.py'],
      description='Download and downscaling tools for the earth2observe datasets',
      data_files = datafiles,
      )

