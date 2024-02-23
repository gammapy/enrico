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
CATALOG_TEMPLATE_DIR =  join(CATALOG_DIR,'Templates')

DIFFUSE_DIR = os.environ.get('FERMI_DIFFUSE_DIR', '')
DOWNLOAD_DIR = os.environ.get('FERMI_DOWNLOAD_DIR', '')
WEEKLY_DIR = ''
WEEKLY_SC_DIR = ''
if DOWNLOAD_DIR :
  WEEKLY_DIR = join(DOWNLOAD_DIR, 'weekly/photon')
  WEEKLY_SC_DIR = join(DOWNLOAD_DIR, 'weekly/spacecraft')
PREPROCESSED_DIR = os.environ.get('FERMI_PREPROCESSED_DIR', '')
CONFIG_DIR = join(os.path.dirname(__file__), 'config')
USE_FULLMISSION_SPACECRAFT = os.environ.get('USE_FULLMISSION_SPACECRAFT','True')#bool(os.environ.get('USE_FULLMISSION_SPACECRAFT','False')=='True')
COMPRESS_WEEKLY_FILES = os.environ.get('COMPRESS_WEEKLY_FILES','False')


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

CATALOG_VERSION = '30'
TEMPLATE_VERSION = '12years'

#if CATALOG_DIR :
#  CATALOG_TEMPLATE_DIR = join(CATALOG_DIR, 'Extended_archive_%s/Templates'% TEMPLATE_VERSION)
#  try :
#    os.mkdir(join(CATALOG_DIR, 'Extended_archive_%s'% TEMPLATE_VERSION))
#  except:
#    pass

IRF_TAG = "P8R3"
ISO_MAJOR = "V3"
ISO_MINOR = "v1"
TAG_ISO="{0}_{1}".format(ISO_MAJOR,ISO_MINOR)
CATALOG = 'gll_psc_v%s.fit' % CATALOG_VERSION
DIFFUSE_GAL = 'gll_iem_v07.fits'
DIFFUSE_ISO_SOURCE = 'iso_{2}_SOURCE_{0}_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEFRONT = 'iso_{2}_SOURCE_{0}_FRONT_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEBACK = 'iso_{2}_SOURCE_{0}_BACK_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEPSF0 = 'iso_{2}_SOURCE_{0}_PSF0_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEPSF1 = 'iso_{2}_SOURCE_{0}_PSF1_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEPSF2 = 'iso_{2}_SOURCE_{0}_PSF2_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEPSF3 = 'iso_{2}_SOURCE_{0}_PSF3_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEEDISP0 = 'iso_{2}_SOURCE_{0}_EDISP0_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEEDISP1 = 'iso_{2}_SOURCE_{0}_EDISP1_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEEDISP2 = 'iso_{2}_SOURCE_{0}_EDISP2_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_SOURCEEDISP3 = 'iso_{2}_SOURCE_{0}_EDISP3_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEAN = 'iso_{2}_CLEAN_{0}_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANFRONT = 'iso_{2}_CLEAN_{0}_FRONT_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANBACK = 'iso_{2}_CLEAN_{0}_BACK_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANPSF0 = 'iso_{2}_CLEAN_{0}_PSF0_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANPSF1 = 'iso_{2}_CLEAN_{0}_PSF1_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANPSF2 = 'iso_{2}_CLEAN_{0}_PSF2_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANPSF3 = 'iso_{2}_CLEAN_{0}_PSF3_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANEDISP0 = 'iso_{2}_CLEAN_{0}_EDISP0_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANEDISP1 = 'iso_{2}_CLEAN_{0}_EDISP1_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANEDISP2 = 'iso_{2}_CLEAN_{0}_EDISP2_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
DIFFUSE_ISO_CLEANEDISP3 = 'iso_{2}_CLEAN_{0}_EDISP3_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEAN = 'iso_{2}_ULTRACLEAN_{0}_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANPSF0 = 'iso_{2}_ULTRACLEAN_{0}_PSF0_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANPSF1 = 'iso_{2}_ULTRACLEAN_{0}_PSF1_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANPSF2 = 'iso_{2}_ULTRACLEAN_{0}_PSF2_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANPSF3 = 'iso_{2}_ULTRACLEAN_{0}_PSF3_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANEDISP0 = 'iso_{2}_ULTRACLEAN_{0}_EDISP0_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANEDISP1 = 'iso_{2}_ULTRACLEAN_{0}_EDISP1_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANEDISP2 = 'iso_{2}_ULTRACLEAN_{0}_EDISP2_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
#DIFFUSE_ISO_ULTRACLEANEDISP3 = 'iso_{2}_ULTRACLEAN_{0}_EDISP3_{1}.txt'.format(ISO_MAJOR,ISO_MINOR,IRF_TAG)
SPACECRAFT = 'lat_spacecraft_merged.fits'

def check_command_line_tools():
    """Check command line tool availability"""
    from subprocess import Popen, PIPE
    print('*** COMMAND LINE TOOLS ***')
    for tool in ['python', 'ipython', 'gtlike', 'enrico_setupcheck']:
        location = Popen(['which', tool],
                         stdout=PIPE).communicate()[0]
        print(('{0:.<20} {1}'.format(tool, location.strip() or 'MISSING')))


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
            print(('{0:.<20} {1}'.format(package, filename)))
        except ImportError:
            print(('{0:.<20} {1}'.format(package, 'MISSING')))


def print_farm():
   """Print the name of the submission farm"""
   print('*** FARM ***')
   if FARM=='':
     print(('{0:.<20} {1}'.format("FARM", 'MISSING')))
   else:
     print(('{0:.<20} {1}'.format("FARM", FARM)))

