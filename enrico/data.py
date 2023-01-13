"""Utilities to download / preprocess data"""
import os
from os.path import join
import glob
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

# @todo: use config file for this?
from enrico.environ import CATALOG_DIR, CATALOG, DIFFUSE_DIR, DIFFUSE_GAL, DIFFUSE_ISO_SOURCE, DIFFUSE_ISO_CLEAN
from enrico.environ import DIFFUSE_ISO_SOURCEFRONT, DIFFUSE_ISO_SOURCEBACK, DIFFUSE_ISO_CLEANFRONT, DIFFUSE_ISO_CLEANBACK
from enrico.environ import DIFFUSE_ISO_SOURCEPSF0, DIFFUSE_ISO_SOURCEPSF1, DIFFUSE_ISO_SOURCEPSF2, DIFFUSE_ISO_SOURCEPSF3
from enrico.environ import DIFFUSE_ISO_SOURCEEDISP0, DIFFUSE_ISO_SOURCEEDISP1, DIFFUSE_ISO_SOURCEEDISP2, DIFFUSE_ISO_SOURCEEDISP3
from enrico.environ import DIFFUSE_ISO_CLEANPSF0, DIFFUSE_ISO_CLEANPSF1, DIFFUSE_ISO_CLEANPSF2, DIFFUSE_ISO_CLEANPSF3
from enrico.environ import DIFFUSE_ISO_CLEANEDISP0, DIFFUSE_ISO_CLEANEDISP1, DIFFUSE_ISO_CLEANEDISP2, DIFFUSE_ISO_CLEANEDISP3
from enrico.environ import DIRS, DOWNLOAD_DIR, CATALOG_TEMPLATE_DIR, TEMPLATE_VERSION, PREPROCESSED_DIR
from enrico.environ import WEEKLY_DIR, SPACECRAFT, USE_FULLMISSION_SPACECRAFT,COMPRESS_WEEKLY_FILES,WEEKLY_SC_DIR
#check platform
from sys import platform 

#TODO: Read from default config
default_filter = 'DATA_QUAL==1&&LAT_CONFIG==1&&ABS(ROCK_ANGLE)<52'

# Download URLs
FSSC_URL = 'http://fermi.gsfc.nasa.gov/ssc'
#FSSC_FTP_URL = 'ftp://legacy.gsfc.nasa.gov/fermi/data/lat'
FSSC_FTP_URL = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/'
HEASARC_FTP = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat'
WEEKLY_DIFFUSE_URL = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat'
CATALOG_URL = join(FSSC_URL, 'data/access/lat/12yr_catalog')
CATALOG_URL_EXTENDED = join(FSSC_URL, 'data/access/lat/12yr_catalog')
# CATALOG_URL_8yr = join(FSSC_URL, 'data/access/lat/fl8y/')
DIFFUSE_URL = join(FSSC_URL, 'data/analysis/software/aux/')
WEEKLY_URL = join(FSSC_FTP_URL, 'weekly/photon')
WEEKLY_SC_URL = join(FSSC_FTP_URL, 'weekly/spacecraft')
WEEKLY_DIFFRSP_URL = join(HEASARC_FTP, 'weekly/diffuse')
SPACECRAFT_URL = join(FSSC_FTP_URL, 'mission/spacecraft', SPACECRAFT)

#        (_tag,                 _url,        _dir,        _file)
FILES = [('CATALOG',            CATALOG_URL, CATALOG_DIR, CATALOG),
         # ('CATALOG',            CATALOG_URL_8yr, CATALOG_DIR, CATALOG_8yr),
         ('DIFFUSE_GAL',        DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_GAL),
         ('DIFFUSE_ISO_SOURCE', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCE),
         ('DIFFUSE_ISO_SOURCEFRONT', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEFRONT),
         ('DIFFUSE_ISO_SOURCEBACK', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEBACK),
         ('DIFFUSE_ISO_SOURCEPSF0', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEPSF0),
         ('DIFFUSE_ISO_SOURCEPSF1', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEPSF1),
         ('DIFFUSE_ISO_SOURCEPSF2', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEPSF2),
         ('DIFFUSE_ISO_SOURCEPSF3', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEPSF3),
         ('DIFFUSE_ISO_CLEAN', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEAN),
         #('DIFFUSE_ISO_CLEANFRONT', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANFRONT),
         #('DIFFUSE_ISO_CLEANBACK', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANBACK),
         ('DIFFUSE_ISO_CLEANPSF0', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANPSF0),
         ('DIFFUSE_ISO_CLEANPSF1', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANPSF1),
         ('DIFFUSE_ISO_CLEANPSF2', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANPSF2),
         ('DIFFUSE_ISO_CLEANPSF3', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANPSF3),
         #('DIFFUSE_ISO_SOURCEEDISP0', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEEDISP0),
         #('DIFFUSE_ISO_SOURCEEDISP1', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEEDISP1),
         #('DIFFUSE_ISO_SOURCEEDISP2', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEEDISP2),
         #('DIFFUSE_ISO_SOURCEEDISP3', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_SOURCEEDISP3),
         #('DIFFUSE_ISO_CLEANEDISP0', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANEDISP0),
         #('DIFFUSE_ISO_CLEANEDISP1', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANEDISP1),
         #('DIFFUSE_ISO_CLEANEDISP2', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANEDISP2),
         #('DIFFUSE_ISO_CLEANEDISP3', DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEANEDISP3),
         #('DIFFUSE_ISO_CLEAN',  DIFFUSE_URL, DIFFUSE_DIR, DIFFUSE_ISO_CLEAN)
]


# fermievtypes = {\
#     'FRONT': 1,
#     'BACK': 2,
#     'PSF0': 4,
#     'PSF1': 8,
#     'PSF2': 16,
#     'PSF3': 32,
#     'EDISP0': 64,
#     'EDISP1': 128,
#     'EDISP2': 256,
#     'EDISP3': 512,
#     }


def check_dirs():
    """Check directory availability"""
    print('*** DIRECTORIES ***')
    for tag in list(DIRS.keys()):
        dir = DIRS.get(tag, 'MISSING')
        status = 'YES' if os.path.isdir(dir) else 'NO'
        print(('{tag:.<20} {status:.<10} {dir}'.format(**locals())))


def check_files():
    """Check file availability"""
    print('*** FILES ***')
    for _tag, _url, _dir, _file in FILES:
        path = join(_dir, _file)
        status = path if os.path.isfile(path) else 'MISSING'
        print(('{0:.<20} {1}'.format(_tag, status)))


def check_catalog_templates():
    """Check catalog template availability"""
    return False


class Data(object):
    """Standard data directories and files"""
    # Event classes pass 7 are still hierarchichal

    # P8R3 IRF name   Event Class (evclass)   Class Hierarchy Photon File Extended File
    # P8R3_SOURCEVETO_V3  2048    Standard    X   X
    # P8R3_ULTRACLEANVETO_V3  1024    Standard    X   X
    # P8R3_ULTRACLEAN_V3  512 Standard    X   X
    # P8R3_CLEAN_V3   256 Standard    X   X
    # P8R3_SOURCE_V3  128 Standard    X   X
    # P8R3_TRANSIENT010_V3    64  Standard        X
    # P8R3_TRANSIENT020_V3    16  Standard        X
    # P8R3_TRANSIENT010E_V3   32  Extended        X
    # P8R3_TRANSIENT020E_V3   8   Extended        X
    EVENT_CLASSES = dict(source=128,
                         clean=256)

    # Time selections
    SELECTIONS = dict(week=1,
                      month=4,
                      lastmonth=-4,
                      eleven_months=48,
                      two_year=105,
                      three_year=156,
                      all=1000)

    # Energy selections
    EMINS = [100, 1000, 10000, 100000]

    # Available preprocessing steps
    STEPS = ['gtselect', 'gtmktime', 'gtltcube']

    def __init__(self, chatter=4, clobber='no', debug='no'):
        self.chatter = chatter
        self.clobber = clobber
        self.debug = debug

    def download(self, spacecraft=True, photon=True):
        """Download spacecraft file and weekly photon files.

        Uses wget to download new data incrementally from
        the Fermi LAT data server.
        For info on wget options see
        http://www.gnu.org/software/wget/manual/wget.html"""
        
        os.system("mkdir -p "+DOWNLOAD_DIR)
        os.chdir(DOWNLOAD_DIR)
        if spacecraft:
            # -m --mirror
            if USE_FULLMISSION_SPACECRAFT=="True":
                if platform == "linux" or platform == "linux2": # linux 
                    cmd = 'wget -N ' + SPACECRAFT_URL
                elif platform == "darwin": 
                    cmd = 'curl -O -N ' + SPACECRAFT_URL
            else:
                list_spacecrafts = sorted(glob.glob("weekly/spacecraft/lat_spacecraft*.fits*"))
                list_of_files=[F.replace('weekly/spacecraft/','').replace('.gz','') for F in list_spacecrafts]
                if platform == "linux" or platform == "linux2": # linux 
                    cmd = 'wget -c -N -R '+','.join(list_of_files[:-1])+' -r -l 1 -A fits --random-wait -e robots=off -nH --cut-dirs=4 -np ' + WEEKLY_SC_URL
                elif platform == "darwin": 
                    cmd = 'curl -O  -c  -r -l 1 -A fits --random-wait -e robots=off -nH --cut-dirs=4 -np ' + WEEKLY_SC_URL
            
                if COMPRESS_WEEKLY_FILES=="True":
                    # always leave out the last file uncompressed, so that wget can resume download
                    list_of_files=sorted(glob.glob("weekly/spacecraft/lat_spacecraft*.fits"))[:-1]
                    for f in list_of_files:
                        cmd = 'gzip -f '+f
                        print(('EXEC: ', cmd))
                        os.system(cmd)
                
            print(('EXEC: ', cmd))
            os.system(cmd)
        if photon:
            # -m --mirror
            # -P --directory-prefix, in this case put in sub-directory weekly
            # -nH --no-host-directories
            # -np --no-parent
            print(WEEKLY_URL)
            list_photons = sorted(glob.glob("weekly/photon/lat_photon_weekly*.fits*"))
            list_of_files=[F.replace('weekly/photon/','').replace('.gz','') for F in list_photons]
            cmd = 'wget -c -N -r -l 1 -A fits --random-wait -e robots=off -R '+','.join(list_of_files[:-1])+' -nH --cut-dirs=4 -np ' + WEEKLY_URL
            print(('EXEC: ', cmd))
            os.system(cmd)
            if COMPRESS_WEEKLY_FILES=="True":
                # always leave out the last file uncompressed, so that wget can resume download
                list_of_files=sorted(glob.glob("weekly/photon/lat_photon_weekly*.fits"))[:-1]
                for f in list_of_files:
                    cmd = 'gzip -f '+f
                    print(('EXEC: ', cmd))
                    os.system(cmd)


    def download_aux(self):
        """Download missing diffuse model and catalog files"""
        from urllib.request import urlretrieve
        from subprocess import call
        # Catalog and diffuse files
        for _tag, _url, _dir, _file in FILES:
            os.system("mkdir -p "+_dir)
            path = join(_dir, _file)
            url = join(_url, _file)
            if _dir:
                if not os.path.isfile(path):
                    print(('Downloading %s -> %s' %(url,path)))
                    urlretrieve(url, path)
                else:
                    print(('%s exists already. Not downloading.' % _file))
            else:
                print('Set DIRECTORIES before downloading the files.')
        # Diffuse emission templates
        if CATALOG_DIR:
            if not os.path.isdir(CATALOG_TEMPLATE_DIR):
                print(('Creating directory %s' % CATALOG_TEMPLATE_DIR))
                os.mkdir(CATALOG_TEMPLATE_DIR)
            # Use any one of the templates to check if they are installed
            filename = join(CATALOG_DIR, CATALOG)
            if not os.path.isfile(filename):
                url = join(CATALOG_URL, filename)
                path = join(CATALOG_DIR, filename)
                print(('Downloading %s -> %s' %(url,path)))
                urlretrieve(url, path)
                
            filename = join(CATALOG_TEMPLATE_DIR, 'CenALobes.fits')
            if not os.path.isfile(filename):
                # Now we know that we want to download
                filename = 'LAT_extended_sources_%s.tgz' % TEMPLATE_VERSION
                url = join(CATALOG_URL_EXTENDED, filename)
                path = join(CATALOG_DIR, filename)
                print(('Downloading %s -> %s' %(url,path)))
                urlretrieve(url, path)
                print('Unpacking')
                call(['tar', 'zxvf', path, '-C', CATALOG_DIR])
            else:
                print('Diffuse emission templates exist already. '
                      'Not downloading.')
        else:
            print('Set CATALOG_DIR before downloading the files.')

    def preprocess(self, steps=None, event_classes=None,
                   selections=None, emins=None):
        """Preprocess data (run gtselect, gtmktim, gtltcube)"""
        if steps == None:
            steps = self.STEPS
        if event_classes == None:
            event_classes = list(self.EVENT_CLASSES.keys())
        if selections == None:
            selections = list(self.SELECTIONS.keys())
        if emins == None:
            emins = self.EMINS
        log.info('steps: %s' % steps)
        log.info('event_classes: %s' % event_classes)
        log.info('selections: %s' % selections)
        log.info('emins: %s' % emins)
        for event_class in event_classes:
            evclass = self.EVENT_CLASSES[event_class]
            for selection in selections:
                for emin in emins:
                    WORKDIR = join(PREPROCESSED_DIR, event_class,
                                   selection, 'emin_%06d' % emin)
                    try:
                        os.makedirs(WORKDIR)
                        log.info('MKDIR: %s' % WORKDIR)
                    except OSError:
                        log.info('%s exists already' % WORKDIR)
                    os.chdir(WORKDIR)
                    # @note: This is not the most efficient way,
                    # gtselect for higher emin could be rerun on
                    # the event file for lower emin instead of
                    # re-processing all weekly files.
                    # Also gtltcube is energy-independent, so
                    # it doesn't have to be re-run for each energy band.
                    # But this only runs once, so what the heck...
                    # also Fermi tools have been known to get confused
                    # by re-running tools in non-standard order.
                    self._preprocess_list(selection)
                    if 'gtselect' in steps:
                        self._preprocess_gtselect(evclass, emin)
                    if 'gtmktime' in steps:
                        self._preprocess_gtmktime()
                    if 'gtltcube' in steps:
                        self._preprocess_gtltcube()

    def _preprocess_list(self, selection):
        """Produce lists of weekly event list files."""
        files = os.listdir(WEEKLY_DIR)
        # Select only fits files
        files = [join(WEEKLY_DIR, _) + '\n'
                 for _ in files if _.endswith('.fits')]
        # Sort chronologically
        files.sort()
        # Select only a subset of weeks
        weeks = self.SELECTIONS[selection]
        if weeks <0:
            files = files[weeks:]
        else:
            files = files[:weeks]
        log.debug('Writing weeks.lis with %04d lines.' % len(files))
        open('weeks.lis', 'w').writelines(files)
        if USE_FULLMISSION_SPACECRAFT=="False":
            """Produce lists of weekly spacecraft list files."""
            files = os.listdir(WEEKLY_SC_DIR)
            # Select only fits files
            files = [join(WEEKLY_SC_DIR, _) + '\n'
                     for _ in files if _.endswith('.fits')]
            # Sort chronologically
            files.sort()
            # Select only a subset of weeks
            # weeks = self.SELECTIONS[selection]
            if weeks <0:
                files = files[weeks:]
            else:
                files = files[:weeks]
            # files = files[:weeks]
            log.debug('Writing weeks.lis with %04d lines.' % len(files))
            open('weeks_sc.lis', 'w').writelines(files)

    def _preprocess_gtselect(self, evclass, emin):
        """Run gtselect"""
        from gt_apps import filter as tool
        self._set_common_tool_options(tool)
        tool['infile'] = '@weeks.lis'
        tool['outfile'] = 'gtselect.fits'
        tool['ra'] = 'INDEF'
        tool['dec'] = 'INDEF'
        tool['rad'] = 'INDEF'
        tool['tmin'] = 'INDEF'
        tool['tmax'] = 'INDEF'
        tool['emin'] = emin
        tool['emax'] = int(1e10)
        tool['zmax'] = 100
        tool['evclass'] = evclass
        tool['convtype'] = -1
        tool['phasemin'] = 0
        tool['phasemax'] = 1
        tool['evtable'] = 'EVENTS'
        tool['clobber'] = 'yes'
        tool.run()

    def _preprocess_gtmktime(self):
        """Run gtmktime"""
        from gt_apps import maketime as tool
        self._set_common_tool_options(tool)
        if USE_FULLMISSION_SPACECRAFT=="True":
            tool['scfile'] = join(DOWNLOAD_DIR, SPACECRAFT)
        else:
            tool['scfile'] =  '@weeks_sc.lis'
        tool['sctable'] = 'SC_DATA'
        tool['filter'] = default_filter
        tool['roicut'] = 'no'
        tool['evfile'] = 'gtselect.fits'
        tool['evtable'] = 'EVENTS'
        tool['outfile'] = 'gtmktime.fits'
        tool['apply_filter'] = 'yes'
        tool['overwrite'] = 'no'
        tool['header_obstimes'] = 'yes'
        tool['gtifile'] = 'default'
        tool.run()

    def _preprocess_gtltcube(self):
        """Run gtltcube"""
        from gt_apps import GtApp
        tool = GtApp('gtltcube')
        self._set_common_tool_options(tool)
        tool['evfile'] = 'gtmktime.fits'
        tool['evtable'] = 'EVENTS'
        if USE_FULLMISSION_SPACECRAFT=="True":
            tool['scfile'] = join(DOWNLOAD_DIR, SPACECRAFT)
        else:
            tool['scfile'] = '@weeks_sc.lis'
        tool['sctable'] = 'SC_DATA'
        tool['outfile'] = 'gtltcube.fits'
        tool['dcostheta'] = 0.025
        tool['binsz'] = 1
        tool['phibins'] = 0
        tool['tmin'] = 0
        tool['tmax'] = 0
        tool['file_version'] = 1
        tool['zmax'] = 180
        tool.run()

    def _set_common_tool_options(self, tool):
        tool['chatter'] = self.chatter
        tool['clobber'] = self.clobber
        tool['debug'] = self.debug
        tool['gui'] = 'no'
        tool['mode'] = 'ql'
