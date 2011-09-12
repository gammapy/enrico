""""xml model utility functions.

At the moment this is a thin wrapper of user_contrib.make2FGLxml,
but maybe some day we'll write our own."""
from os.path import join
from enrico.user_contrib.make2FGLxml import addSrcs2
import enrico.environ as env 

class SourceList:
    """Create and xml model file from the 2FGL catalog"""
    def __init__(self, config):
        """Set all options"""
        self.out = join(config['out'], config['file']['xml'])
        self.clobber = config['tool']['clobber']
        self.ra = config['target']['ra']
        self.dec = config['target']['dec']
        self.rad = config['space']['rad']
        self.srcrad = config['space']['srcrad']
        self.psF = config['model']['point_only']
        self.radLim = config['model']['max_radius']
        self.roi = (self.ra, self.dec, self.rad)
        self.srcs = join(env.CATALOG_DIR, env.CATALOG)
        self.extD = env.CATALOG_TEMPLATE_DIR        
        self.GDname = config['model']['diffuse_gal']
        self.GDfile = join(env.DIFFUSE_DIR, self.GDname + '.fits')
        self.ISOname = config['model']['diffuse_iso']
        self.ISOfile = join(env.DIFFUSE_DIR, self.ISOname + '.txt')
        self.signif = config['model']['min_significance']
    def makeModel(self):
        """Construct the source list and write to self.out"""
        print('Writing %s' % self.out)
        addSrcs2(self,
                 self.GDfile, 
                 self.GDname,
                 self.ISOfile, 
                 self.ISOname, 
                 self.signif)
