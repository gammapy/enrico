"""Central place for analysis environment setup
and standard data file locations."""
import os
from os.path import join
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
from extern.odict import OrderedDict as dict


# Directory names
ENRICO_DIR = os.environ.get('ENRICO_DIR', '')
FERMI_DATA_DIR = os.environ.get('FERMI_DATA_DIR', '')
CATALOG_DIR = os.environ.get('FERMI_CATALOG_DIR', '')
CATALOG_TEMPLATE_DIR = join(CATALOG_DIR, 'Templates') if CATALOG_DIR else ''
DIFFUSE_DIR = os.environ.get('FERMI_DIFFUSE_DIR', '')
DOWNLOAD_DIR = os.environ.get('FERMI_DOWNLOAD_DIR', '')
WEEKLY_DIR = join(DOWNLOAD_DIR, 'weekly/photon') if DOWNLOAD_DIR else ''
PREPROCESSED_DIR = os.environ.get('FERMI_PREPROCESSED_DIR', '')
CONFIG_DIR = join(os.path.dirname(__file__), 'data', 'config')
XML_DIR = join(os.path.dirname(__file__), 'data', 'xml')

DIRS = dict(FERMI_DATA_DIR=FERMI_DATA_DIR,
            CATALOG_DIR=CATALOG_DIR,
            DIFFUSE_DIR=DIFFUSE_DIR,
            PREPROCESSED_DIR=PREPROCESSED_DIR,
            DOWNLOAD_DIR=DOWNLOAD_DIR,
            WEEKLY_DIR=WEEKLY_DIR,
            CONFIG_DIR= CONFIG_DIR,
            XML_DIR=XML_DIR,
            ENRICO_DIR=ENRICO_DIR)

# File names
CATALOG_VERSION = '06'
CATALOG = 'gll_psc_v%s.fit' % CATALOG_VERSION
DIFFUSE_GAL = 'gal_2yearp7v6_v0.fits'
DIFFUSE_ISO_SOURCE = 'iso_p7v6source.txt'
DIFFUSE_ISO_CLEAN = 'iso_p7v6clean.txt'
SPACECRAFT = 'lat_spacecraft_merged.fits'

# Download URLs
FSSC_URL = 'http://fermi.gsfc.nasa.gov/ssc'
FSSC_FTP_URL = 'ftp://legacy.gsfc.nasa.gov/fermi/data/lat'
CATALOG_URL = join(FSSC_URL, 'data/access/lat/2yr_catalog')
DIFFUSE_URL = join(FSSC_URL, 'data/analysis/software/aux')
WEEKLY_URL = join(FSSC_FTP_URL, 'weekly/photon')
SPACECRAFT_URL = join(FSSC_FTP_URL, 'mission/spacecraft',
                      SPACECRAFT)

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
    # http://stackoverflow.com/questions/2617704/checking-for-module-availability-programmatically-in-python
    ASTRO_PACKAGES = ['pyfits', 'kapteyn']
    # @todo: Here it's enough to try one of the Fermi python modules 
    # and to show where they are located.
    FERMI_PACKAGES = ['gt_apps', 'UnbinnedAnalysis']
    PACKAGES = ['enrico', 'IPython'] + ASTRO_PACKAGES + FERMI_PACKAGES
    print('*** PYTHON PACKAGES ***')
    for package in PACKAGES:
        try:
            exec('import %s' % package)
            filename = eval('%s.__file__' % package)
            print('{0:.<20} {1}'.format(package, filename))
        except ImportError:
            print('{0:.<20} {1}'.format(package, 'MISSING'))
