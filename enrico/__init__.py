"""Enrico helps you with your Fermi data analysis."""
import os
from os.path import join
import logging
logger = logging.getLogger(__name__)
import urllib

# @todo: How to make this easy to set up for the user?
DATA_DIR = '.'
DIFFUSE_DIR = '.'
CATALOG_DIR = '.'

try:
    FERMI_DIR = os.environ['FERMI_DIR']
except:
    FERMI_DIR = None

def check_installation():
    """Check if software, data, diffuse models,
    catalog are available"""
    print('FERMI_DIR: %s' % FERMI_DIR)

def download_diffuse():
    """Download missing diffuse model files"""
    filenames = ['gal_2yearp7v6_v0.fits', 'iso_p7v6source.txt', 'iso_p7v6clean.txt']
    for filename in filenames[1:]:
        path = join(DIFFUSE_DIR, filename)
        url = 'http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/' + filename
        if not os.path.isfile(path):
            logger.info('Downloading %s to %s' % (url, path))
            urllib.urlretrieve(url, path)