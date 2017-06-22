"""Central place for analysis environment setup
and standard data file locations."""
import os
from os.path import join
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
#from enrico.extern.odict import OrderedDict


#Submission farm name
#Currently supported : LAPP-Annecy, MPIK-HD, CCIN2P3, LOCAL-relaxedtimes
#FARM = os.environ.get('FARM','MPIK')
FARM  = os.environ.get('FARM','LOCAL')
QUEUE = os.environ.get('QUEUE','batch')
TORQUE_RESOURCES = os.environ.get('TORQUE_RESOURCES','')

# Directory names
ENRICO_DIR = os.environ.get('ENRICO_DIR', '')
FERMI_DIR = os.environ.get('FERMI_DIR', '')
FERMI_DATA_DIR = os.environ.get('FERMI_DATA_DIR', '')
CATALOG_DIR = os.environ.get('FERMI_CATALOG_DIR', '')
CATALOG_TEMPLATE_DIR = ''

DIFFUSE_DIR = os.environ.get('FERMI_DIFFUSE_DIR', '')
DOWNLOAD_DIR = os.environ.get('FERMI_DOWNLOAD_DIR', '')
WEEKLY_DIR = ''
if DOWNLOAD_DIR :
  WEEKLY_DIR = join(DOWNLOAD_DIR, 'weekly/photon')
PREPROCESSED_DIR = os.environ.get('FERMI_PREPROCESSED_DIR', '')
CONFIG_DIR = join(os.path.dirname(__file__), 'config')


try :
    from enrico.extern.odict import OrderedDict
    DIRS = OrderedDict(FERMI_DATA_DIR=FERMI_DATA_DIR,
                   FERMI_DIR=FERMI_DIR,
                   CATALOG_DIR=CATALOG_DIR,
                   DIFFUSE_DIR=DIFFUSE_DIR,
                   PREPROCESSED_DIR=PREPROCESSED_DIR,
                   DOWNLOAD_DIR=DOWNLOAD_DIR,
                   WEEKLY_DIR=WEEKLY_DIR,
                   CONFIG_DIR=CONFIG_DIR,
                   ENRICO_DIR=ENRICO_DIR)
except :
    DIRS = {}

# File names
CATALOG_VERSION = '16'
TEMPLATE_VERSION = '15'
if CATALOG_DIR :
  CATALOG_TEMPLATE_DIR = join(CATALOG_DIR, 'Extended_archive_v%s/Templates'% TEMPLATE_VERSION)
  try :
    os.mkdir(join(CATALOG_DIR, 'Extended_archive_v%s'% TEMPLATE_VERSION))
  except:
    pass

CATALOG = 'gll_psc_v%s.fit' % CATALOG_VERSION
DIFFUSE_GAL = 'gll_iem_v06.fits'
DIFFUSE_ISO_SOURCE = 'iso_P8R2_SOURCE_V6_v06.txt'
DIFFUSE_ISO_SOURCE = 'iso_P8R2_SOURCE_V6_v06.txt'
DIFFUSE_ISO_SOURCEPSF0 = 'iso_P8R2_SOURCE_V6_PSF0_v06.txt'
DIFFUSE_ISO_SOURCEPSF1 = 'iso_P8R2_SOURCE_V6_PSF1_v06.txt'
DIFFUSE_ISO_SOURCEPSF2 = 'iso_P8R2_SOURCE_V6_PSF2_v06.txt'
DIFFUSE_ISO_SOURCEPSF3 = 'iso_P8R2_SOURCE_V6_PSF3_v06.txt'
DIFFUSE_ISO_SOURCEEDISP0 = 'iso_P8R2_SOURCE_V6_EDISP0_v06.txt'
DIFFUSE_ISO_SOURCEEDISP1 = 'iso_P8R2_SOURCE_V6_EDISP1_v06.txt'
DIFFUSE_ISO_SOURCEEDISP2 = 'iso_P8R2_SOURCE_V6_EDISP2_v06.txt'
DIFFUSE_ISO_SOURCEEDISP3 = 'iso_P8R2_SOURCE_V6_EDISP3_v06.txt'
DIFFUSE_ISO_CLEAN = 'iso_P8R2_CLEAN_V6_v06.txt'
DIFFUSE_ISO_CLEANPSF0 = 'iso_P8R2_CLEAN_V6_PSF0_v06.txt'
DIFFUSE_ISO_CLEANPSF1 = 'iso_P8R2_CLEAN_V6_PSF1_v06.txt'
DIFFUSE_ISO_CLEANPSF2 = 'iso_P8R2_CLEAN_V6_PSF2_v06.txt'
DIFFUSE_ISO_CLEANPSF3 = 'iso_P8R2_CLEAN_V6_PSF3_v06.txt'
DIFFUSE_ISO_CLEANEDISP0 = 'iso_P8R2_CLEAN_V6_EDISP0_v06.txt'
DIFFUSE_ISO_CLEANEDISP1 = 'iso_P8R2_CLEAN_V6_EDISP1_v06.txt'
DIFFUSE_ISO_CLEANEDISP2 = 'iso_P8R2_CLEAN_V6_EDISP2_v06.txt'
DIFFUSE_ISO_CLEANEDISP3 = 'iso_P8R2_CLEAN_V6_EDISP3_v06.txt'
SPACECRAFT = 'lat_spacecraft_merged.fits'

def check_command_line_tools():
    """Check command line tool availability"""
    from subprocess import Popen, PIPE
    print('*** COMMAND LINE TOOLS ***')
    for tool in ['python', 'ipython', 'gtlike', 'enrico_setupcheck']:
        location = Popen(['which', tool],
                         stdout=PIPE).communicate()[0]
        print('{0:.<20} {1}'.format(tool, location.strip() or 'MISSING'))


def check_python_modules():
    """Check python package availability"""
    # @todo: Use this fast method to check for python module availability:
    # http://stackoverflow.com/questions/2617704/
    # checking-for-module-availability-programmatically-in-python
    ASTRO_PACKAGES = ['pyfits', 'kapteyn']
    # @todo: Here it's enough to try one of the Fermi python modules
    # and to show where they are located.
    FERMI_PACKAGES = ['gt_apps', 'UnbinnedAnalysis','BinnedAnalysis']
    PACKAGES = ['enrico', 'IPython'] + ASTRO_PACKAGES + FERMI_PACKAGES
    print('*** PYTHON PACKAGES ***')
    for package in PACKAGES:
        try:
            exec('import %s' % package)
            filename = eval('%s.__file__' % package)
            print('{0:.<20} {1}'.format(package, filename))
        except ImportError:
            print('{0:.<20} {1}'.format(package, 'MISSING'))


def print_farm():
   """Print the name of the submission farm"""
   print('*** FARM ***')
   if FARM=='':
     print('{0:.<20} {1}'.format("FARM", 'MISSING'))
   else:
     print('{0:.<20} {1}'.format("FARM", FARM))

