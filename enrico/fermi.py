"""Helper functions to call Fermi tools
@todo: Old stuff, can be removed?"""
import os.path
import logging
from config import get_default_config
config = get_default_config()
SCFILE = config['file']['spacecraft']
LTCUBE = config['file']['ltcube']

ZMAX = config['analysis']['zmax']
EVCLSMIN = config['analysis']['evclsmin']
EVCLSMAX = config['analysis']['evclsmax']
FILTER = config['analysis']['filter']
IRFS = config['analysis']['irfs']
BINSZ = config['space']['binsz']
NXPIX = config['space']['nxpix']
NYPIX = config['space']['nypix']
COORDSYS = config['space']['coordsys']
XREF = config['space']['xref']
YREF = config['space']['yref']
PROJ = config['space']['proj']
EMIN = config['energy']['emin']
EMAX = config['energy']['emax']
ENUMBINS = config['energy']['nenergies']


# -------------------------------------------
# General utility functions
# -------------------------------------------

def call(cmd, echo=True):
    """Call a command line tool"""
    from os import system
    cmd = [str(_) for _ in cmd]
    if echo:
        print ' '.join(cmd)
    system(' '.join(cmd))

def mkdir(path):
    """Make a dir if it doesn't exist already"""
    if not os.path.exists(path):
        logging.info('Creating dir: {0}'.format(path))
        os.mkdir(dir)
    else:
        logging.debug('Existing dir: {0}'.format(path))

def quote(text, quote_char="'"):
    """Add quotes to a string"""
    return quote_char + text + quote_char

def etag(emin, emax, fmt='%07d'):
    return ('emin_%s_emax_%s' % (fmt, fmt)) % (emin, emax)

#------------------------------------------------------------
# Helper functions to run Fermi tools
#------------------------------------------------------------

def execute_command(cmd, logfile=None, add_common_commands=True):
    """A helper function to call tools. Some common options
    are set via global options in this module"""
    if add_common_commands:
        cmd += ' chatter=%d' % config['tool']['chatter']
        cmd += ' clobber=%s' % config['tool']['clobber']
        cmd += ' debug=%s' % config['tool']['debug']
        cmd += ' gui=%s' % config['tool']['gui']
        cmd += ' mode=%s' % config['tool']['mode']
    if logfile:
        cmd += ' 2>&1 > ' + logfile    
    print cmd
    from os import system
    system(cmd)

def run_select(infile, outfile,
               ra='INDEF', dec='INDEF', rad='INDEF',
               tmin='INDEF', tmax='INDEF',
               emin=EMIN, emax=EMAX, zmax=ZMAX,
               evclsmin=EVCLSMIN, evclsmax=EVCLSMAX,
               convtype=-1, phasemin=0, phasemax=1,
               logfile=None):
    """Run gtselect"""
    cmd = 'gtselect '
    cmd += ' infile=%s' % infile
    cmd += ' outfile=%s' % outfile
    cmd += ' ra=%s' % ra
    cmd += ' dec=%s' % dec
    cmd += ' rad=%s' % rad
    cmd += ' tmin=%s' % tmin
    cmd += ' tmax=%s' % tmax
    cmd += ' emin=%g' % emin
    cmd += ' emax=%g' % emax
    cmd += ' zmax=%d' % zmax
    cmd += ' evclsmin=%s' % evclsmin
    cmd += ' evclsmax=%s' % evclsmax
    cmd += ' evclass=INDEF'
    cmd += ' convtype=%s' % convtype
    cmd += ' phasemin=%s' % phasemin
    cmd += ' phasemax=%s' % phasemax    
    execute_command(cmd, logfile)

def run_mktime(scfile, evfile, outfile, 
               filter=FILTER, roicut='no', header_obstimes='yes',
               logfile=None):
    """Run gtmktime"""
    cmd = 'gtmktime '
    cmd += ' scfile=%s' % scfile
    cmd += ' filter="%s"' % filter
    cmd += ' roicut=%s' % roicut
    cmd += ' header_obstimes=%s' % header_obstimes
    cmd += ' evfile=%s' % evfile
    cmd += ' outfile=%s' % outfile
    execute_command(cmd, logfile)
    
def run_ltcube(evfile, scfile, outfile, 
               dcostheta=0.025, binsz=1, zmax=ZMAX,
               tmin=0, tmax=0,
               logfile=None):
    """Run gtltcube"""
    cmd = 'gtltcube'
    cmd += ' evfile=%s' %evfile
    cmd += ' scfile=%s' % scfile
    cmd += ' outfile=%s' % outfile
    cmd += ' dcostheta=%s' % dcostheta
    cmd += ' binsz=%s' % binsz
    cmd += ' zmax=%s' % zmax
    cmd += ' tmin=%s' % tmin
    cmd += ' tmax=%s' % tmax
    execute_command(cmd, logfile)

def run_bin(evfile, outfile, scfile=SCFILE,
            binsz=BINSZ, nxpix=NXPIX, nypix=NYPIX,
            coordsys=COORDSYS, xref=XREF, yref=YREF, proj=PROJ,
            emin=EMIN, emax=EMAX, enumbins=ENUMBINS,
            logfile=None):
    """Run gtbin"""
    cmd = 'gtbin '
    cmd += ' algorithm=CCUBE'
    cmd += ' ebinalg=LOG'
    cmd += ' evfile=%s' % evfile
    cmd += ' outfile=%s' % outfile
    cmd += ' scfile=%s' % scfile
    cmd += ' nxpix=%d' % nxpix
    cmd += ' nypix=%d' % nypix
    cmd += ' binsz=%g' % binsz
    cmd += ' coordsys=%s' % coordsys
    cmd += ' xref=%s' % xref
    cmd += ' yref=%s' % yref
    cmd += ' axisrot=0'
    cmd += ' proj=%s' % proj
    cmd += ' emin=%g' % emin
    cmd += ' emax=%g' % emax
    cmd += ' enumbins=%d' % enumbins
    execute_command(cmd, logfile)
    
def run_expcube(infile, evfile, cmfile, outfile, irfs=IRFS,
                nxpix=NXPIX, nypix=NYPIX, pixscale=BINSZ, 
                coordsys=COORDSYS, xref=XREF, yref=YREF, proj=PROJ,
                emin=EMIN, emax=EMAX, enumbins=ENUMBINS, 
                bincalc='EDGE', logfile=None):
    """Run gtexpcube"""
    cmd = 'gtexpcube '
    cmd += ' infile=%s' % infile
    cmd += ' evfile=%s' % evfile
    cmd += ' cmfile=%s' % cmfile
    cmd += ' outfile=%s' % outfile
    cmd += ' irfs=%s' % irfs
    cmd += ' nxpix=%d' % nxpix
    cmd += ' nypix=%d' % nypix
    cmd += ' pixscale=%g' % pixscale
    cmd += ' coordsys=%s' % coordsys
    cmd += ' xref=%s' % xref
    cmd += ' yref=%s' % yref
    cmd += ' axisrot=0'
    cmd += ' proj=%s' % proj
    cmd += ' emin=%g' % emin
    cmd += ' emax=%g' % emax
    cmd += ' enumbins=%d' % enumbins
    cmd += ' bincalc=%s' % bincalc
    execute_command(cmd, logfile)

def run_srcmaps(cmap, srcmdl, bexpmap, outfile,
                scfile=SCFILE, expcube=LTCUBE, irfs=IRFS,
                convol='yes', resample='yes', rfactor=2,
                ptsrc='yes', psfcorr='yes', logfile=None):
    """Run gtsrcmaps"""
    cmd = 'gtsrcmaps '
    cmd += ' scfile=%s' % scfile 
    cmd += ' expcube=%s' % expcube  
    cmd += ' cmap=%s' % cmap 
    cmd += ' srcmdl=%s' % srcmdl 
    cmd += ' bexpmap=%s' % bexpmap 
    cmd += ' outfile=%s' % outfile 
    cmd += ' irfs=%s' % irfs 
    cmd += ' convol=%s' % convol 
    cmd += ' resample=%s' % resample
    cmd += ' rfactor=%s' % rfactor
    cmd += ' ptsrc=%s' % ptsrc
    cmd += ' psfcorr=%s' % psfcorr
    execute_command(cmd, logfile)

def run_model(srcmaps, srcmdl, outfile,
              irfs=IRFS, expcube=LTCUBE, bexpmap='none',
              convol='yes', resample='yes', rfactor=2,
              logfile=None):
    """Run gtmodel"""
    cmd = 'gtmodel '
    cmd += ' srcmaps=%s' % srcmaps 
    cmd += ' srcmdl=%s' % srcmdl 
    cmd += ' outfile=%s' % outfile
    cmd += ' irfs=%s' % irfs  
    cmd += ' expcube=%s' % expcube  
    cmd += ' bexpmap=%s' % bexpmap 
    cmd += ' convol=%s' % convol 
    cmd += ' resample=%s' % resample
    cmd += ' rfactor=%s' % rfactor
    execute_command(cmd, logfile)

def run_obssim(infile, srclist, simtime, scfile='none', evroot='mc',
               ltfrac=1.0, tstart=239557417.,
               nevents="no", maxtime=3e8, startdate='2001-01-01',
               offset=0, use_ac='no', ra=0, dec=0, radius=20,
               emin=EMIN, emax=EMAX, edisp='no', irfs=IRFS,
               area=1, maxrows=int(1e8), seed=0, logfile=None):
    """Run gtobssim"""
    cmd = 'gtobssim '
    cmd += ' infile=%s' % infile
    cmd += ' srclist=%s' % srclist
    cmd += ' scfile=%s' % scfile
    cmd += ' evroot=%s' % evroot
    cmd += ' simtime=%s' % simtime
    cmd += ' ltfrac=%s' % ltfrac
    cmd += ' tstart=%s' % tstart
    cmd += ' nevents=%s' % nevents
    cmd += ' maxtime=%s' % maxtime
    cmd += ' startdate=%s' % startdate
    cmd += ' offset=%s' % offset
    cmd += ' use_ac=%s' % use_ac
    cmd += ' ra=%s' % ra
    cmd += ' dec=%s' % dec
    cmd += ' radius=%s' % radius
    cmd += ' emin=%s' % emin
    cmd += ' emax=%s' % emax
    cmd += ' edisp=%s' % edisp
    cmd += ' irfs=%s' % irfs
    cmd += ' area=%s' % area
    cmd += ' maxrows=%s' % maxrows
    cmd += ' seed=%s' % seed
    execute_command(cmd, logfile)

#------------------------------------------------------------
# Helper functions to run FTOOLS
#------------------------------------------------------------

def run_fgauss(infile, outfile, sigma,
               ratio=1, theta=0, nsigma=3,
               boundary='nearest', constant=0,
               datatype='E', nullval=0, 
               copyprime='yes', copyall='no',
               logfile=None):
    """Run fgauss"""
    cmd = 'fgauss '
    cmd += ' infile=%s' % infile 
    cmd += ' outfile=%s' % outfile
    cmd += ' sigma=%s' % sigma 
    cmd += ' ratio=%s' % ratio 
    cmd += ' theta=%s' % theta
    cmd += ' nsigma=%s' % nsigma
    cmd += ' boundary=%s' % boundary
    cmd += ' constant=%s' % constant
    cmd += ' datatype=%s' % datatype
    cmd += ' nullval=%s' % nullval
    cmd += ' copyprime=%s' % copyprime
    cmd += ' copyall=%s' % copyall
    execute_command(cmd, logfile, add_common_commands=False)
