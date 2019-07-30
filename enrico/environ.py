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
WEEKLY_SC_DIR = ''
if DOWNLOAD_DIR :
  WEEKLY_DIR = join(DOWNLOAD_DIR, 'weekly/photon')
  WEEKLY_SC_DIR = join(DOWNLOAD_DIR, 'weekly/spacecraft')
PREPROCESSED_DIR = os.environ.get('FERMI_PREPROCESSED_DIR', '')
CONFIG_DIR = join(os.path.dirname(__file__), 'config')
USE_FULLMISSION_SPACECRAFT = bool(os.environ.get('USE_FULLMISSION_SPACECRAFT','False')=='True')

try :
    from enrico.extern.odict import OrderedDict
    DIRS = OrderedDict(FERMI_DATA_DIR=FERMI_DATA_DIR,
                   FERMI_DIR=FERMI_DIR,
                   CATALOG_DIR=CATALOG_DIR,
                   DIFFUSE_DIR=DIFFUSE_DIR,
                   PREPROCESSED_DIR=PREPROCESSED_DIR,
                   DOWNLOAD_DIR=DOWNLOAD_DIR,
                   WEEKLY_DIR=WEEKLY_DIR,
                   WEEKLY_SC_DIR=WEEKLY_SC_DIR,
                   CONFIG_DIR=CONFIG_DIR,
                   ENRICO_DIR=ENRICO_DIR)
except :
    DIRS = {}

# File names
CATALOG_VERSION = '19'
CATALOG_8yr_VERSION = '6'
TEMPLATE_VERSION = '8years' #v15
if CATALOG_DIR :
  CATALOG_TEMPLATE_DIR = join(CATALOG_DIR, 'Extended_archive_%s/Templates'% TEMPLATE_VERSION)
  try :
    os.mkdir(join(CATALOG_DIR, 'Extended_archive_%s'% TEMPLATE_VERSION))
  except:
    pass

CATALOG = 'gll_psc_v%s.fit' % CATALOG_VERSION
CATALOG_8yr = 'gll_psc_8year_v%s.fit' % CATALOG_8yr_VERSION
DIFFUSE_GAL = 'gll_iem_v07.fits'
DIFFUSE_ISO_SOURCE = 'iso_P8R3_SOURCE_V2_v1.txt'
DIFFUSE_ISO_SOURCE = 'iso_P8R3_SOURCE_V2_v1.txt'
#DIFFUSE_ISO_SOURCEFRONT = 'iso_P8R3_SOURCE_V2_FRONT_v1.txt'
#DIFFUSE_ISO_SOURCEBACK = 'iso_P8R3_SOURCE_V2_BACK_v1.txt'
DIFFUSE_ISO_SOURCEPSF0 = 'iso_P8R3_SOURCE_V2_PSF0_v1.txt'
DIFFUSE_ISO_SOURCEPSF1 = 'iso_P8R3_SOURCE_V2_PSF1_v1.txt'
DIFFUSE_ISO_SOURCEPSF2 = 'iso_P8R3_SOURCE_V2_PSF2_v1.txt'
DIFFUSE_ISO_SOURCEPSF3 = 'iso_P8R3_SOURCE_V2_PSF3_v1.txt'
DIFFUSE_ISO_SOURCEEDISP0 = 'iso_P8R3_SOURCE_V2_EDISP0_v1.txt'
DIFFUSE_ISO_SOURCEEDISP1 = 'iso_P8R3_SOURCE_V2_EDISP1_v1.txt'
DIFFUSE_ISO_SOURCEEDISP2 = 'iso_P8R3_SOURCE_V2_EDISP2_v1.txt'
DIFFUSE_ISO_SOURCEEDISP3 = 'iso_P8R3_SOURCE_V2_EDISP3_v1.txt'
DIFFUSE_ISO_CLEAN = 'iso_P8R3_CLEAN_V2_v1.txt'
DIFFUSE_ISO_CLEANPSF0 = 'iso_P8R3_CLEAN_V2_PSF0_v1.txt'
DIFFUSE_ISO_CLEANPSF1 = 'iso_P8R3_CLEAN_V2_PSF1_v1.txt'
DIFFUSE_ISO_CLEANPSF2 = 'iso_P8R3_CLEAN_V2_PSF2_v1.txt'
DIFFUSE_ISO_CLEANPSF3 = 'iso_P8R3_CLEAN_V2_PSF3_v1.txt'
DIFFUSE_ISO_CLEANEDISP0 = 'iso_P8R3_CLEAN_V2_EDISP0_v1.txt'
DIFFUSE_ISO_CLEANEDISP1 = 'iso_P8R3_CLEAN_V2_EDISP1_v1.txt'
DIFFUSE_ISO_CLEANEDISP2 = 'iso_P8R3_CLEAN_V2_EDISP2_v1.txt'
DIFFUSE_ISO_CLEANEDISP3 = 'iso_P8R3_CLEAN_V2_EDISP3_v1.txt'
#DIFFUSE_ISO_ULTRACLEAN = 'iso_P8R3_ULTRACLEAN_V2_v1.txt'
#DIFFUSE_ISO_ULTRACLEANPSF0 = 'iso_P8R3_ULTRACLEAN_V2_PSF0_v1.txt'
#DIFFUSE_ISO_ULTRACLEANPSF1 = 'iso_P8R3_ULTRACLEAN_V2_PSF1_v1.txt'
#DIFFUSE_ISO_ULTRACLEANPSF2 = 'iso_P8R3_ULTRACLEAN_V2_PSF2_v1.txt'
#DIFFUSE_ISO_ULTRACLEANPSF3 = 'iso_P8R3_ULTRACLEAN_V2_PSF3_v1.txt'
#DIFFUSE_ISO_ULTRACLEANEDISP0 = 'iso_P8R3_ULTRACLEAN_V2_EDISP0_v1.txt'
#DIFFUSE_ISO_ULTRACLEANEDISP1 = 'iso_P8R3_ULTRACLEAN_V2_EDISP1_v1.txt'
#DIFFUSE_ISO_ULTRACLEANEDISP2 = 'iso_P8R3_ULTRACLEAN_V2_EDISP2_v1.txt'
#DIFFUSE_ISO_ULTRACLEANEDISP3 = 'iso_P8R3_ULTRACLEAN_V2_EDISP3_v1.txt'
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
    ASTRO_PACKAGES = ['astropy', 'kapteyn']
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

