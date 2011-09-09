#!/usr/bin/env python
"""
Runs the example analysis from the Fermi ScienceTools python tutorial:
http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/python_tutorial.html

See also
http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/python_usage_notes.html

- Can be used to test the installation after a software update.
- Note that all output goes to the current working directory.

@note: The sections on making the TS map and the butterfly are not implemented.
"""
import sys; sys.path.append('/Users/deil/work/workspace/galpop/src/scripts/fermi/user_contrib')
import logging
from os.path import join
import numpy as np
import gt_apps
import pyLikelihood # @note If this is removed I get crashes!
from UnbinnedAnalysis import UnbinnedObs, UnbinnedAnalysis

def check_limits(like, verbose=True):
    """Check if some of the parameters are almost at
    the allowed limit.
    This is useful to identify fit convergence problems."""
    n_pars_at_limit = 0
    for src in like.sourceNames():
        for name in like.model[src].funcs['Spectrum'].paramNames:
            bounds = like.model[src].funcs['Spectrum'].getParam(name).getBounds()
            value = like.model[src].funcs['Spectrum'].getParam(name).value()
            for bound in bounds:
                if( abs(bound - value) < 0.01):
                    n_pars_at_limit += 1
                    if verbose:
                        print src,name,value,bound
    return n_pars_at_limit


class UnbinnedAnalysisRunner:
    """The ScienceTools user-contributed likeSED requires that
    the the event list, exposure files and xml model have already
    been produced. You can do this with this class."""

    def __init__(self, options={}):
        self.options = options

    def setup(self, tool):
        """This method is used to set all parameters the tool accepts
        from self.options. This is convenient, because now we don't
        have to pass every single parameter to every single tool
        explicitly. If something is set to a value it shouldn't we
        can always set it explicity after before calling run()."""
        for key, value in self.options.items():
            try:
                tool[key] = value
            except KeyError:
                pass

    def filter(self):
        logging.info('filter')
        self.setup(gt_apps.filter)
        gt_apps.filter['infile'] = self.options['EVENTS_UNFILTERED']
        gt_apps.filter['outfile'] = self.options['EVENTS_FILTERED']
        gt_apps.filter.run()

    def maketime(self):
        logging.info('maketime')
        self.setup(gt_apps.maketime)
        gt_apps.maketime['evfile'] = self.options['EVENTS_FILTERED']
        gt_apps.maketime['outfile'] = self.options['EVENTS']
        gt_apps.maketime.run()

    def expCube(self):
        logging.info('expCube')
        self.setup(gt_apps.expCube)
        gt_apps.expCube['evfile'] = self.options['EVENTS']
        gt_apps.expCube['outfile'] = self.options['expcube']
        gt_apps.expCube['dcostheta'] = 0.025
        gt_apps.expCube['binsz'] = 1
        gt_apps.expCube.run()

    def expMap(self):
        logging.info('expMap')
        self.setup(gt_apps.expMap)
        gt_apps.expMap['evfile'] = self.options['EVENTS']
        gt_apps.expMap['outfile'] = self.options['expmap']
        gt_apps.expMap.run()
    
    def make_xml(self):
        """Use the user-contributed make2FGLxml script and 2FGL
        catalog to make an xml model file.
        The same argument names as in make2FGLxml are used."""
        logging.info('make_xml')
        import make2FGLxml
        CATALOG_DIR = '/Users/deil/work/workspace/survey/catalogs/fermi/'
        EXTENDED_DIR = join(CATALOG_DIR, 'Templates')
        DIFFUSE_DIR = '/Users/deil/bin/fermi/ScienceTools-v9r23p1-fssc-20110726/external/diffuseModels/'
        # Make a srcList object
        sources = join(CATALOG_DIR, 'gll_psc_v04.fit')
        ft1 = self.options['EVENTS']
        out = self.options['xmlmodel']
        source_list = make2FGLxml.srcList(sources, ft1, out)
        # Call the srcList.makeModel method
        GDfile = join(DIFFUSE_DIR, 'gal_2yearp7v6_v0.fits')
        GDname = 'gal_2yearp7v6_v0' 
        ISOfile = join(DIFFUSE_DIR, 'iso_p7v6source.txt') 
        ISOname = 'iso_p7v6source'
        extDir = EXTENDED_DIR
        radLim = 5 # ROI center offset beyond which all spectral parameters are fixed (-1 meas = ROI radius) 
        signif = 10 # Sources with lower significance will have spectral parameters fixed
        psForce = False # Treat exteded sources as point sources? 
        source_list.makeModel(GDfile, GDname, ISOfile, ISOname,
                              extDir, radLim, signif, psForce)
        # @todo Implement more flexible criteria which parameters
        # to fix based on source offset, flux and significance!!!
    
    def load_obs(self):
        logging.info('load_obs')
        eventFile = self.options['EVENTS']
        scFile = self.options['scfile']
        expMap = self.options['expmap']
        expCube = self.options['expcube']
        irfs = self.options['irfs']
        self.obs = UnbinnedObs(eventFile, scFile, expMap, expCube, irfs)
    
    def run_fit1(self, verbosity=10):
        logging.info('run_fit1')
        srcModel = self.options['xmlmodel']
        optimizer = 'DRMNFB'
        logging.info('Creating like1 ... this takes a few minutes')
        self.like1 = UnbinnedAnalysis(self.obs, srcModel, optimizer)
        like = self.like1
        like.tol = 1e-1
        logging.info('Fitting ... this takes a few minutes')
        like.fit(verbosity=verbosity)
        print like.model['_2FGLJ1555.7+1111'] 
        like.logLike.writeXml(self.options['xmlfit1'])
        #like.plot()
        
    def run_fit2(self, verbosity=10):
        logging.info('run_fit2')
        srcModel = self.options['xmlfit1']
        optimizer = 'NewMinuit'
        logging.info('Creating like2 ... this takes about two minutes')
        self.like2 = UnbinnedAnalysis(self.obs, srcModel, optimizer)
        like = self.like2
        like.tol = 1e-2
        self.like2obj = pyLikelihood.Minuit(like.logLike)
        logging.info('Fitting ... this takes a few minutes')
        like.fit(verbosity=verbosity,covar=True,optObject=self.like2obj)
        print 'fit 2 return code:', self.like2obj.getRetCode()
        like.tolType = 0
        like.fit(verbosity=verbosity,covar=True,optObject=self.like2obj)
        print 'fit 2 return code second try:', self.like2obj.getRetCode()
        like.logLike.writeXml(self.options['xmlfit2'])
        print like.model['_2FGLJ1555.7+1111']
        emin, emax = self.options['emin'], self.options['emax']
        print('Integral flux: %s +- %s', 
              (like.flux('_2FGLJ1555.7+1111', emin, emax),
               like.fluxError('_2FGLJ1555.7+1111', emin, emax)))
        ts =  like.Ts('_2FGLJ1555.7+1111')
        print 'TS:', ts
        print 'Significance:', np.sqrt(ts) 
        
    def run_likeSED(self):
        logging.info('run_likeSED')
        from likeSED import likeInput, likeSED
        logging.info('Creating like3 ... this takes about two minutes')
        srcModel = self.options['xmlfit2']
        optimizer = 'NewMinuit'
        like = UnbinnedAnalysis(self.obs, srcModel, optimizer)
        SrcName = '_2FGLJ1555.7+1111'
        model = 'PG1553_fit2.xml'
        low_edges = [200.,427.69,914.61,1955.87,4182.56,8944.27,19127.05,40902.61]
        high_edges = [427.69,914.61,1955.87,4182.56,8944.27,19127.05,40902.61,187049.69]
        nbins = len(low_edges)
        inputs = likeInput(like, SrcName,
                           model, nbins)
        logging.info('customBins ... this will take about 10 minutes')
        inputs.customBins(low_edges, high_edges)
        # inputs.plotBins() # @todo What does this do?
        logging.info('fullFit')
        inputs.fullFit()
        
        sed=likeSED(inputs)
        sed.getECent()
        sed.fitBands()
        sed.Plot()

    def check_results(self):
        logging.info('check_results')
        pass

def download_data(dir='.'):
    url = 'http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/data/pyLikelihood/'
    files = ['L110728114426E0D2F37E85_PH00.fits',
             'L110728114426E0D2F37E85_PH01.fits',
             'L110728114426E0D2F37E85_SC00.fits']
    logging.info('Trying to download data files from %s' % url)
    import urllib
    raise NotImplementedError


def run_python_tutorial_example():
    # Configuration parameters
    DATA_DIR = '/Users/deil/work/_Data/fermi/python_tutorial/'
    files = dict(EVENTS_UNFILTERED = '@' + join(DATA_DIR, 'PG1553_events.list'),
                 scfile = join(DATA_DIR, 'PG1553_SC.fits'),
                 EVENTS_FILTERED = 'PG1553_filtered.fits',
                 EVENTS = 'PG1553_filtered_gti.fits',
                 expcube = 'PG1553_ltCube.fits',
                 expmap = 'PG1553_expMap.fits',
                 xmlmodel = 'PG1553_model.xml',
                 xmlfit1 = 'PG1553_fit1.xml',
                 xmlfit2 = 'PG1553_fit2.xml')
    parameters = dict(evclass=2, ra=238.929, dec=11.1901,
                      rad=10, emin=200, emax=400000,
                      zmax=100, roicut='yes',
                      srcrad=20, nlong=120, nlat=120,
                      nenergies=20,
                      tmin=239557417, tmax=256970880,
                      #tmin=239557417, tmax=239557417 + 7*24*3600,
                      filter='(DATA_QUAL==1)&&(LAT_CONFIG==1)&&ABS(ROCK_ANGLE)<52',
                      irfs='P7SOURCE_V6')
    options = dict(chatter=4, clobber='yes')
    options.update(files)
    options.update(parameters)
    
    # Run the analysis
    #download_data(DATA_DIR)
    runner = UnbinnedAnalysisRunner(options)
    #runner.filter()
    #runner.maketime()
    #runner.expCube()
    #runner.expMap()
    runner.make_xml()
    
    runner.load_obs()
    runner.run_fit1()
    runner.run_fit2()
    runner.run_likeSED()
    return runner.like2

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.info('main')    
    like = run_python_tutorial_example()
