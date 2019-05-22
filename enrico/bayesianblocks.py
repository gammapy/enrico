import os
from math import sqrt
import numpy as np
import scipy.optimize
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 15})
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt
from enrico import utils
from enrico import plotting
from enrico import environ
from enrico import lightcurve
from enrico.config import get_config
from enrico.constants import LightcurvePath,FoldedLCPath
from enrico.submit import call
from enrico.RunGTlike import run, GenAnalysisObjects
from enrico import Loggin
from enrico.plotting import plot_errorbar_withuls

pol0 = lambda x,p1: p1
pol1 = lambda x,p1,p2: p1+p2*x


class BayesianBlocks(lightcurve.LightCurve):
    """Class to calculate light curves and variability indexes."""
    def __init__(self, config, parent_filename=""):
        super().__init__(config, parent_filename)

        # Convert time if necessary
        if self.config['time']['type']=='MJD':
             self.config['time']['tmin'] = utils.MJD_to_met(self.config['time']['tmin'])
             self.config['time']['tmax'] = utils.MJD_to_met(self.config['time']['tmax'])
             self.config['time']['type']=='MET'
        elif self.config['time']['type']=='JD':
             self.config['time']['tmin'] = utils.JD_to_met(self.config['time']['tmin'])
             self.config['time']['tmax'] = utils.JD_to_met(self.config['time']['tmax'])
             self.config['time']['type']=='MET'
        self.tmin = self.config['time']['tmin']
        self.tmax = self.config['time']['tmax']

        self.p0 = self.config['BayesianBlocks']['p0']
        self.config['Spectrum']['FitsGeneration'] = self.config['BayesianBlocks']['FitsGeneration']
        self.config['Spectrum']['FrozenSpectralIndex'] = self.config['BayesianBlocks']['SpectralIndex']
        self.config['UpperLimit']['TSlimit'] = self.config['BayesianBlocks']['TSLightCurve']

        # Check apperture light curve have benn run first
        self._CheckAppertureLightCurveFile()

    def _CheckAppertureLightCurveFile(self):
        ''' Check the existance of apperture light curve file (the selected evt list with good gti)'''
        import os.path
        evtfile = str("%s/AppertureLightCurve/%s_%s_MkTime.fits"%(self.folder,self.srcname,self.Tag))
        if not os.path.isfile(evtfile):
            raise Exception('The apperture photometry events list doesn\'t exist\nPlease run enrico_applc first')

    def _MakeTimeBins(self):
        from astropy.stats import bayesian_blocks
        from astropy.table import Table

        evtfile = str("%s/AppertureLightCurve/%s_%s_MkTime.fits"%(self.folder,self.srcname,self.Tag))
        evtlist = Table.read(evtfile, hdu='EVENTS')

        meanRate = float(len(evtlist))/float(self.tmax-self.tmin)
        print("Mean photon rate %s s^-1" %meanRate)
        print("Mean photon rate %s day^-1" %(meanRate*3600.*24))

        edges = bayesian_blocks(evtlist['TIME'], fitness='events', p0=self.p0)

        self.Nbin = len(edges)-1
        self.time_array = np.zeros(self.Nbin*2)
        self.gtifile = []
        for i in xrange(self.Nbin):
            self.time_array[2*i] = edges[i]
            self.time_array[2*i+1]= edges[i+1]

        self.info("Running LC with "+str(self.Nbin)+" bins")
        for i in xrange(self.Nbin):
            print "Bin ",i," Start=",self.time_array[2*i]," Stop=",self.time_array[2*i+1]
        print

    def _ManageFolder(self,path):
        """   All files will be stored in a subfolder name path + NLCbin
        Create a subfolder"""
        self.LCfolder =  self.folder+"/BayesianBlocks/"
        utils.mkdir_p(self.LCfolder)
        self.config['out'] = self.LCfolder

    def _MakeLC(self,Path=LightcurvePath) :
        #import gc
        import os
        #gc.enable()
        '''Main function of the Lightcurve script. Read the config file and run the gtlike analysis'''
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')

        self.PrepareLC(self.config['BayesianBlocks']['MakeConfFile'])#Get the config file

        for i in xrange(self.Nbin):
            #gc.collect()
            cmd = str("enrico_sed %s && enrico_plot_bayesianblocks %s" %(self.configfile[i], self.parent_filename))
            if self.submit == 'yes':
                scriptname = self.LCfolder+"LC_Script_"+str(i)+".sh"
                JobLog = self.LCfolder+"LC_Job_"+str(i)+".log"
                JobName = (self.config['target']['name'] + "_" +
                       self.config['analysis']['likelihood'] +
                       "_LC_" + self.config['file']['tag'])+"_"+str(i)+".log"

                call(cmd,enricodir,fermidir,scriptname,JobLog,JobName)#Submit the job
            else :
                os.system(cmd)

                #run(self.configfile[i])#run in command line
