"""Central place for analysis environment setup."""
import os
from os.path import join
#import logging
#logging.basicConfig(level=logging.INFO)
#log = logging.getLogger(__name__)
import subprocess
import urllib

# @todo: Optionally use config file for directories!

# Directory names
FERMI_DIR = os.environ.get('FERMI_DIR', '')
CATALOG_DIR = os.environ.get('FERMI_CATALOG_DIR', '')
CATALOG_TEMPLATE_DIR = join(CATALOG_DIR, 'Templates') if CATALOG_DIR else ''
# @note: At the moment the Fermi tools don't include all diffuse models,
# but if a future version does we could simply use those:
#DIFFUSE_DIR_FERMI = join(FERMI_DIR, 'refdata/fermi/diffuseModels') if FERMI_DIR else ''
#DIFFUSE_DIR = os.environ.get('FERMI_DIFFUSE_DIR', DIFFUSE_DIR_FERMI)
DIFFUSE_DIR = os.environ.get('FERMI_DIFFUSE_DIR', '')
DOWNLOAD_DIR = os.environ.get('FERMI_DOWNLOAD_DIR', '')
DATA_DIR = os.environ.get('FERMI_DATA_DIR', '')
CONFIG_DIR = join(os.path.dirname(__file__), 'data', 'config')
XML_DIR = join(os.path.dirname(__file__), 'data', 'xml')
#       (_tag,          _dir)
DIRS = [('FERMI_DIR',   FERMI_DIR),
        ('CATALOG_DIR', CATALOG_DIR),
        ('DIFFUSE_DIR', DIFFUSE_DIR),
        ('DATA_DIR',    DATA_DIR),
        ('DOWNLOAD_DIR',DOWNLOAD_DIR),
        ('CONFIG_DIR',  CONFIG_DIR),
        ('XML_DIR',     XML_DIR)]

# File names
#CATALOG = join(CATALOG_DIR, 'gll_psc_v05.fit')
#DIFFUSE_GAL = join(DIFFUSE_DIR, 'gal_2yearp7v6_v0.fits')
#DIFFUSE_ISO_SOURCE = join(DIFFUSE_DIR, 'iso_p7v6source.txt')
#DIFFUSE_ISO_CLEAN = join(DIFFUSE_DIR, 'iso_p7v6clean.txt')
CATALOG = 'gll_psc_v05.fit'
DIFFUSE_GAL = 'gal_2yearp7v6_v0.fits'
DIFFUSE_ISO_SOURCE = 'iso_p7v6source.txt'
DIFFUSE_ISO_CLEAN = 'iso_p7v6clean.txt'
SPACECRAFT = 'lat_spacecraft_merged.fits'

# Download URLs
CATALOG_URL = 'http://fermi.gsfc.nasa.gov/ssc/data/access/lat/2yr_catalog'
DIFFUSE_URL = 'http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux'
DATA_URL = 'TODO'

#        (_tag,                 _url,        _dir,        _file)
FILES = [('CATALOG',            CATALOG_URL, CATALOG_DIR, CATALOG),
         ('DIFFUSE_GAL',        DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_GAL),
         ('DIFFUSE_ISO_SOURCE', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCE),
         ('DIFFUSE_ISO_CLEAN',  DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEAN)]

# Data selections
#             (tag,     #weeks)
SELECTIONS = [('week',       1),
              ('month',      4),
              ('11month',   48),
              ('24month',  105),
              ('36month',  156),
              ('all',     1000)]

def check_dirs():
    """Check directory availability"""
    print('*** DIRECTORIES ***')
    for key, value in DIRS:
        print('{0:.<20} {1}'.format(key, value or 'MISSING'))

def check_files():
    """Check file availability"""
    print('*** FILES ***')
    for _tag, _url, _dir, _file in FILES:
        path = join(_dir, _file)
        status = path if os.path.isfile(path) else 'MISSING'
        print('{0:.<20} {1}'.format(_tag, status))

def check_catalog_templates():
    """Check catalog template availability"""
    return False
    
def check_tools():
    """Check command line tool availability"""
    print('*** COMMAND LINE TOOLS ***')
    for tool in ['python', 'ipython', 'gtlike', 'enrico_setup']:
        location = subprocess.Popen(['which', tool], stdout=subprocess.PIPE).communicate()[0]
        print('{0:.<20} {1}'.format(tool, location.strip() or 'MISSING'))

def check_python():
    """Check python package availability"""
    # @todo: Use this fast method to check for python module availability:
    # http://stackoverflow.com/questions/2617704/checking-for-module-availability-programmatically-in-python
    ASTRO_PACKAGES = ['pyfits', 'kapteyn']
    # @todo: Here it's enough to try one of the Fermi python modules and to show where they are located.
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

def check():
    """Check everything"""
    check_dirs()
    check_files()
    check_tools()
    check_python()
    print('''
       If something is missing, please install it or adjust your environment
       variables (PATH, PYTHONPATH, FERMI_...) as described in the setup tutorial.
    ''')

def download_aux():
    """Download missing diffuse model and catalog files"""
    # Catalog and diffuse files
    for _tag, _url, _dir, _file in FILES:
        path = join(_dir, _file)
        url = join(_url, _file)
        if _dir:
            if not os.path.isfile(path):
                print('Downloading %s' % path)
                urllib.urlretrieve(url, path)
            else:
                print('%s exists already. Not downloading.' % _file)
        else:
            print('Set DIRECTORIES before downloading the files.')

    # Diffuse emission templates
    if CATALOG_DIR:
        if not os.path.isdir(CATALOG_TEMPLATE_DIR):
            print('Creating directory %s' % CATALOG_TEMPLATE_DIR)
            os.mkdir(CATALOG_TEMPLATE_DIR)
        # Use any one of the templates to check if they are installed
        filename = join(CATALOG_TEMPLATE_DIR, 'CenALobes.fits')
        if not os.path.isfile(filename):
            # Now we know that we want to download
            filename = 'gll_psc_v05_templates.tgz'
            url = join(CATALOG_URL, filename)
            path = join(CATALOG_DIR, filename)
            print('Downloading %s' % path)
            urllib.urlretrieve(url, path)
            print('Unpacking')
            subprocess.call(['tar', 'zxvf', path, '-C', CATALOG_DIR])
        else:
            print('Diffuse emission templates exist already. Not downloading.')
    else:
        print('Set CATALOG_DIR before downloading the files.')
            

def download_data(spacecraft=True, photon=True):
    """Download spacecraft file and weekly photon files.

    Uses wget to download new data incrementally from 
    the Fermi LAT data server.
    For info on wget options see
    http://www.gnu.org/software/wget/manual/wget.html"""
    os.chdir(DOWNLOAD_DIR)
    if spacecraft:
        # -m --mirror
        cmd = ('wget -N ftp://legacy.gsfc.nasa.gov/FTP/fermi/data/'
               'lat/mission/spacecraft/lat_spacecraft_merged.fits')
        print(cmd)
        os.system(cmd)
    if photon:
        # -m --mirror
        # -P --directory-prefix, in this case put in sub-directory weekly
        # -nH --no-host-directories
        # -np --no-parent
        print(cmd)
        cmd = ('wget -m -P weekly -nH --cut-dirs=4 -np '
               'ftp://legacy.gsfc.nasa.gov/fermi/data/lat/weekly/photon/')
        os.system(cmd)

def prepare_data():
