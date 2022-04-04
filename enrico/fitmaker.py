"""
fitmaker.py written by David Sanchez : david.sanchez@lapp.in2p3.fr
Collection of functions to run the fit (gtlike)
the class FitMaker will call the function of the observation class (gtfunction.py) and prepare the fit
by computing the fits file.
it can distinguish between the binned and unbinnned analysis
begun September 2011
"""
#import logging
#logging.basicConfig(level=logging.INFO)
#log = logging.getLogger(__name__)
import numpy as np
import string
try:
    import astropy.io.fits as fits
except ImportError:
    import pyfits as fits
from UnbinnedAnalysis import UnbinnedAnalysis, UnbinnedObs
from BinnedAnalysis import BinnedAnalysis, BinnedObs, BinnedConfig
from enrico import utils
from enrico import Loggin
from enrico import environ
import pyLikelihood as pyLike

class FitMaker(Loggin.Message):
    """Collection of functions to prepare/run the GTLIKE fit
     and compute an upper limit is needed"""
    def __init__(self, obs, config):
        super(FitMaker,self).__init__()
        Loggin.Message.__init__(self)
        self.obs = obs
        self.config = config
        self.task_number = 1
        self.log_like = 0

    def _log(self, task='', description=''):
        print()
        print(("\033[34m"+'# ' + '*' * 60))
        if task:
            task = '%10s --- ' % task
        print(("\033[34m"+'# *** %3d %s%s' %
              (self.task_number, task, description)))
        print(("\033[34m"+'# ' + '*' * 60+"\033[0m"))
        self.task_number += 1

    def FirstSelection(self,config=None):
        """Make a coarse selection of events from original file"""
        self._log('gtselect', 'Select data from library, coarse cut')#run gtselect
        if config!=None:
            self.obs.Configuration = config
            self.obs.LoadConfiguration()
            self.obs.FirstCut()
            self.obs.Configuration = self.config
            self.obs.LoadConfiguration()
        else:
            self.obs.FirstCut()

    def GenerateFits(self):
        """Run the different ST tools and compute the fits files
           First it runs the tools that are common to the binned
           and unbinned analysis chain then it run the specific
           tools following the choise of the user"""

        #Run the tools common to binned and unbinned chain
        self._log('gtselect', 'Select data from library, fine cut')#run gtselect
        self.obs.SelectEvents()
        self._log('gtmktime', 'Update the GTI and cut data based on ROI')#run gtmktime
        self.obs.MkTime()
        if (self.config["analysis"]["ComputeDiffrsp"] == "yes" and self.config["analysis"]["likelihood"] == "unbinned"):
            self._log('gtdiffrsp', 'Compute Diffuse response')
            self.obs.DiffResps()#run gtdiffresp
        self._log('gtbin', 'Create count maps (square fully embed in the ROI circle)')
        self.obs.Gtbin()
        self._log('gtltcube', 'Make live time cube')#run gtltcube
        self.obs.ExpCube()

        self.info('Compute the psf')#run gtpsf
        self.obs.GtPSF()


        #Choose between the binned of the unbinned analysis
        if self.config['analysis']['likelihood'] == 'binned': #binned analysis chain
            self._log('gtbin', 'Make count map CCUBE')#run gtbin
            self.obs.GtCcube()
            self._log('gtexpcube2', 'Make binned exposure cube')#run gtexpcube2
            self.obs.GtBinnedMap()
            self._log('gtsrcmap', 'Make a source map')#run gtsrcmap
            self.obs.SrcMap()
            self.info('Compute the psf')#run gtpsf
            self.obs.GtPSF()
            self.info('Compute Energy Migration Matrices')#run gtpsf
            self.obs.GtDRM()
            # the model should be generated for each component after optimizing.
            #self._log('gtmodel', 'Make a model map')#run gtmodel
            #self.obs.ModelMap(self.config["file"]["xml"])
            self.info('Compute Energy Migration Matrices')#run gtdrm
            self.obs.GtDRM()


        if self.config['analysis']['likelihood'] == 'unbinned': #unbinned analysis chain
            self._log('gtexpmap', 'Make an exposure map')
            self.obs.ExpMap()
    #the function ends here. It does not run gtlike

    def CreateLikeObject(self):
        """Create an UnbinnedAnalysis or a BinnedAnalysis and retrun it."""
        #create binnedAnalysis object
        if self.config['analysis']['likelihood'] == 'binned':
            use_edisp  = self.config['analysis']['EnergyDispersion'] == 'yes'
            edisp_bins = -2 if use_edisp==True else 0
            Obs = BinnedObs(srcMaps=self.obs.srcMap,
                            expCube=self.obs.Cubename,
                            binnedExpMap=self.obs.BinnedMapfile,
                            irfs=self.obs.irfs)
            Cfg = BinnedConfig(use_edisp=use_edisp, 
                               edisp_bins=edisp_bins)
            Fit = BinnedAnalysis(Obs, self.obs.xmlfile, config=Cfg,
                                 optimizer=self.config['fitting']['optimizer'])
            
            Fit.setEnergyRange(self.obs.Emin,self.obs.Emax)
            #print(("Is edisp enabled? {0}".format(str(Fit.logLike.use_edisp()))))
            
        #create a unbinnedAnalysis object
        if self.config['analysis']['likelihood'] == 'unbinned':
            Obs = UnbinnedObs(self.obs.mktimefile,
                              self.obs.ft2,
                              expMap=self.obs.Mapname,
                              expCube=self.obs.Cubename,
                              irfs=self.obs.irfs)
            Fit = UnbinnedAnalysis(Obs, self.obs.xmlfile,
                                   optimizer=self.config['fitting']['optimizer'])
        
        # Fix this, EBL absorbed models use LogParabola with b=0 instead of PowerLaw, 
        # we may want to allow fixed shape for that case
        if float(self.config['Spectrum']['FrozenSpectralIndex']!=0): 
            parameters = dict()
            parameters['Index']  = -float(self.config['Spectrum']['FrozenSpectralIndex'])
            parameters['alpha']  = +float(self.config['Spectrum']['FrozenSpectralIndex'])
            parameters['Index1'] = -float(self.config['Spectrum']['FrozenSpectralIndex'])
            parameters['beta']   = 0
            parameters['Index2'] = 2.
            parameters['Cutoff'] = 30000. # set the cutoff to be high

            for key in list(parameters.keys()):
                IdGamma = utils.getParamIndx(Fit, self.obs.srcname, key)
                if (IdGamma == -1):
                    continue
                else:
                    self.info("Freezing %s = %s" %(str(key),str(parameters[key])))
                    Fit[IdGamma] = parameters[key] # set the parameter
                    Fit[IdGamma].setFree(False)#the variable index is frozen to compute the UL
                     
        return Fit #return the BinnedAnalysis or UnbinnedAnalysis object.

    def GetOptObject(self,Fit,method):
        if method == 'MINUIT':
            optObject = pyLike.Minuit(Fit.logLike)
        elif method == 'NEWMINUIT':
            optObject = pyLike.NewMinuit(Fit.logLike)
        else:
            optFactory = pyLike.OptimizerFactory_instance()
            optObject = optFactory.create(str(method), Fit.logLike)
        return optObject

    def GetStatus(self,Fit,optObject):
        status = optObject.getRetCode()
        if isinstance(optObject, pyLike.Minuit) == 'MINUIT':
            edm     = optObject.getDistance()
            quality = optObject.getQuality() if status==0 else 2
        elif isinstance(optObject, pyLike.NewMinuit):
            edm     = optObject.getDistance()
            check   = np.bitwise_and(np.bitwise_not(156), status)
            quality = 3 if check==0 else 0
        else:
            edm     = optObject.getDistance()
            quality = 3 if status==0 else 0
        
        loglike = optObject.stat().value()
        
        return(edm,quality,loglike)

    def _PerformFit(self, Fit, methods=["MINUIT"], covar=False):
        """Perform a fit with the selected methods and extract the edm, quality of fit, etc."""
        for method in methods:
            try:
                optObject = self.GetOptObject(Fit,self.config['fitting']['optimizer'])
                self.log_like = Fit.fit(verbosity=0, covar=True, optimizer=method, optObject=optObject)
                edm,quality,loglike = self.GetStatus(Fit,optObject)
                if self.config['verbose'] == 'yes' :
                    print(('Fit output with {1}: {0} [quality: {2}]'.format(
                        self.log_like,self.config['fitting']['optimizer'],quality))) 
            except Exception as exc:
                self.warning('Exception while running {0}, {1}'.format(method,str(exc)))
            else:
                break

        return(Fit)

    def PerformFit(self, Fit, writeXml = True):
        """Run gtlile tool. First it run gtlike with the DRNMGB optimizer
        and the user optimizer after. A dictionnay is return with all
        the releveant results"""

        self._log('gtlike', 'Run likelihood analysis')
        # TODO fix this part

        # Change the fit tolerance to the one given by the user
        Fit.ftol = float(self.config['fitting']['ftol'])
        
        # Fit with DRMNGB/DRMNFB (as recommended in gtlike fhelp) optimizer to 
        # obtain initial parameters close to the minimum. Then switch optimizer.
        list_of_methods = ['DRMNFB','DRMNGB']
        Fit = self._PerformFit(Fit, list_of_methods, covar=False)
        
        # 2nd precise fit with the user optimizer and ask gtlike to compute the covariance matrix
        if self.config['fitting']['optimizer'] in ["MINUIT","NEWMINUIT"]:
            list_of_methods = ["NEWMINUIT","MINUIT"]
        else:
            list_of_methods = [self.config['fitting']['optimizer']]
       
        Fit = self._PerformFit(Fit, list_of_methods, covar=True) 

        # remove source with TS<min_source_TS (default=1)
        # to improve convergence and re-fit
        try:             self.config['fitting']['min_source_TS']
        except KeyError: self.config['fitting']['min_source_TS'] = 1.

        Fit = self.RemoveWeakSources(Fit, self.config['fitting']['min_source_TS'])
        self._log('Re-optimize', '')
        self._PerformFit(Fit, list_of_methods, covar=True) 
        
        self.outXml = None
        if writeXml :
            self.outXml = utils._dump_xml(self.config)
            Fit.writeXml(self.outXml)

        self.success("Fit with gtlike performed")

    def RemoveWeakSources(self,Fit,minTS=1.0):
        """Remove weak sources and re-optimize to get a better minimum."""
        self._log('','Remove all the weak (TS<%.2f) sources' %minTS)
        SrcTsDict = dict()

        for src in Fit.model.srcNames:
            SrcTsDict[src] = Fit.Ts(src)
            if (SrcTsDict[src]<minTS) and not(src == self.obs.srcname):
                for comp in Fit.components:
                    if comp.logLike.getSource(src).getType() == 'Point':
                        if self.config['verbose'] == 'yes' :
                            self.info("deleting source "+src+" with TS = "+\
                                      str(SrcTsDict[src])+" from the model")
                        comp.deleteSource(src)

        for src in Fit.model.srcNames:
            if (SrcTsDict[src]>=minTS):
                if self.config['verbose'] == 'yes' :
                    self.info('keeping source {0} with TS={1:.2e}'.format(src,SrcTsDict[src]))
        

        return Fit

    def GetAndPrintResults(self, Fit):
        """Get and print some useful results. Also contruct a dictonnary and fill it with results"""
        if self.config['verbose'] == 'yes' :
            self._log('Results', 'Print results of the fit')
        Result = {}

        if self.config['verbose'] == 'yes' :
            print((Fit.model,"\n"))
            self.info("Results for the Fit")
            # Print src name, Npred and TS for source with TS > 5
            print("Source Name\tNpred\tTS")
            #TODO
        #for src in Fit.model.srcNames:
        #if Fit.Ts(src) > 5:
        #   print src, "\t%2.3f\t%2.3f" % (Fit.NpredValue(src), Fit.Ts(src))

        # fill the dictonnary with some values
        Result['Optimizer'] = self.config['fitting']['optimizer']
        Result['Npred'] = Fit.NpredValue(self.obs.srcname)
        Result['TS'] = Fit.Ts(self.obs.srcname)
        if self.config['verbose'] == 'yes' :
            print("Values and (MINOS) errors for " + self.obs.srcname)
            print("TS : ", Fit.Ts(self.obs.srcname))

        # Get the python object 'Spectrum' for the source of interest
        spectrum = Fit[self.obs.srcname].funcs['Spectrum']
        # Get the names of the parameters for the source of interest
        ParName = spectrum.paramNames
        #Get the model type and fill the dictonnary
        stype = Fit.model.srcs[self.obs.srcname].spectrum().genericName()
        Result['ModelType'] = stype
        Result['log_like'] = Fit.logLike.value()
        #Add the energy information to the result dictionnary
        Result['Emin'] = self.obs.Emin
        Result['Emax'] = self.obs.Emax
        #Add the time information to the result dictionnary
        Result['tmin'] = self.config['time']['tmin']
        Result['tmax'] = self.config['time']['tmax']
        Result['SrcName'] = self.obs.srcname
        Result['Flux'] = Fit.flux(self.obs.srcname,self.obs.Emin,self.obs.Emax)
        Result['dFlux'] = Fit.fluxError(self.obs.srcname,self.obs.Emin,self.obs.Emax)
        Result['EFlux'] = Fit.energyFlux(self.obs.srcname,self.obs.Emin,self.obs.Emax)
        Result['dEFlux'] = Fit.energyFluxError(self.obs.srcname,self.obs.Emin,self.obs.Emax)
        
        for par in ParName : #Loop over the parameters and get value, error and scale
            ParValue = spectrum.getParam(par).value()
            ParError = spectrum.getParam(par).error()
            Scale    = spectrum.getParam(par).getScale()
            Result[par] = ParValue * Scale
            Result['d'+par] = ParError * Scale
            if ParError>0: # Compute MINOS errors for relevent parameters  Fit.Ts(self.obs.srcname) > 5 and
                try:
                    MinosErrors = Fit.minosError(self.obs.srcname, par)
                    if self.config['verbose'] == 'yes' :
                       print((par+" :  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ] %2.0e" %
                          (ParValue, ParError, MinosErrors[0], MinosErrors[1], Scale)))
                       Result.update({'d'+par+'-': MinosErrors[0] * Scale})
                       Result.update({'d'+par+'+': MinosErrors[1] * Scale})
                except:
                    if self.config['verbose'] == 'yes' :
                        print((par+" :  %2.2f +/-  %2.2f  %2.0e" %
                          (ParValue, ParError, Scale)))
            else:
                if self.config['verbose'] == 'yes' :
                    print((par+" :  %2.2f   %2.0e" %
                      (ParValue, Scale)))

        try: # get covariance matrix
            if self.config['verbose'] == 'yes' :
                utils.GetCovar(self.obs.srcname, Fit)
        except:
            pass #if the covariance matrix has not been computed

        if self.config['verbose'] == 'yes' :
            utils.GetFluxes(Fit,self.obs.Emin,self.obs.Emax) #print the flux of all the sources

        #Compute an UL if the source is too faint
        if float(self.config['UpperLimit']['TSlimit']) > Fit.Ts(self.obs.srcname):
            if self.config['UpperLimit']['envelope'] == 'yes':
                self.EnvelopeUL(Fit)
            else:
                Ulval = self.ComputeUL(Fit)
                Result['Ulvalue'] = Ulval

        return Result   #Return the dictionnary


    def PoissonUL(self,Fit):
        """ Compute UL using Feldman-cousin poisson stat"""
        self.info('Compute the exposure')#run gtexposure

        try :
            spfile = fits.open(self.obs.lcfile)
        except:
            self.obs.GtLCbin(dt = self.config['time']['tmax']-self.config['time']['tmin'])
            #spfile = pyfits.open(self.obs.lcfile)

        try:
            self.obs.Configuration['AppLC']['index'] = self.config['UpperLimit']['SpectralIndex']
        except:
            try:
                self.obs.Configuration['AppLC']['index']
            except:
                self.info("Cannot find the spectral index")
                self.obs.Configuration['AppLC']['index'] = 1.5

        self.info("Assuming spectral index of %s" %self.info("Assuming default index of 1.5"))
        self.obs.GtExposure()
        Exposure = np.sum(spfile[1].data.field("EXPOSURE"))

        ### Run it in any case (useful for instance for the DL3)
        #self.info('Compute the psf')#run gtpsf
        #self.obs.GtPSF()

        ccube = fits.open(self.obs.ccube)
        psfres = fits.open(self.obs.psf)

        #read psf and get the 68% containement radius
        theta = (psfres[2].data["Theta"])
        theta68 = np.zeros(len(psfres[1].data["psf"]))
        for i in range(len(psfres[1].data["psf"])):
            integral = np.trapz(psfres[1].data["psf"][i],theta)
            for j in range(psfres[1].data["psf"][i].size):
               if np.trapz(psfres[1].data["psf"][i][:j],theta[:j])/integral>0.68:
                   theta68[i] = theta[j]
                   break

        #read the CCUBE
        x = np.arange(-abs(int(ccube[0].header["CRPIX1"])*ccube[0].header["CDELT1"]),\
                abs(ccube[0].header["CRPIX1"]*ccube[0].header["CDELT1"]),abs(ccube[0].header["CDELT1"]))
        y = np.arange(-abs(ccube[0].header["CRPIX2"]*ccube[0].header["CDELT2"]),\
                abs(int(ccube[0].header["CRPIX2"])*ccube[0].header["CDELT2"]),abs(ccube[0].header["CDELT2"]))

        xx, yy = np.meshgrid(x, y)
        dist = np.sqrt(xx**2 + yy**2)
        Obsevt = 0 #compute the number of events within the PSF radius
        for i in range(len(psfres[1].data["psf"])):
            maps = ccube[0].data[i]
            Obsevt += sum(maps[dist<theta68[i]])/0.68

        nbg = max(0,int(Obsevt-Fit.NpredValue(self.obs.srcname)))
        Obsevt = int(Fit.NpredValue(self.obs.srcname)+nbg)
        if Obsevt> 20:
            self.warning("Observed Numbers of event too high (>20)\n abort and return -1")
            return -1

        cl = str(int(float(self.config['UpperLimit']['cl'])*100))
        try :
            ullookup = np.genfromtxt(environ.ENRICO_DIR+'/enrico/extern/UL_poisson_'+cl+'.dat',unpack=True)
        except:
            self.warning("cannot find the file "+environ.ENRICO_DIR+'/enrico/extern/UL_poisson_'+cl+'.dat')
        bkglookup = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0])
        measurement = ullookup[0] #the first row is the measurement

        uls = ullookup[2:-1:2] #keep only 1 row over 2 since we don't care about LL

        self.info("Found "+str(Obsevt)+" events for "+str(nbg)+" background event ")
#        print uls[bkglookup.searchsorted(nbg)][measurement.searchsorted(Obsevt)]," ",Exposure
        return uls[bkglookup.searchsorted(nbg)][measurement.searchsorted(Obsevt)]/Exposure

    def ComputeUL(self, Fit):
        """Compute an Upper Limit using either the profil or integral method
        See the ST cicerone for more information on the 2 method"""

        self._log('UpperLimit', 'Compute upper Limit')
        #Index given by the user
        self.info("Assumed index is "+str(self.config['UpperLimit']['SpectralIndex']))

        parameters = dict()
        parameters['Index']  = -float(self.config['UpperLimit']['SpectralIndex'])
        parameters['alpha']  = +float(self.config['UpperLimit']['SpectralIndex'])
        parameters['Index1'] = -float(self.config['UpperLimit']['SpectralIndex'])
        parameters['beta']   = 0
        parameters['Index2'] = 2.
        parameters['Cutoff'] = 30000. # set the cutoff to be high

        for key in list(parameters.keys()):
            try:
                utils.FreezeParams(Fit,self.obs.srcname, key, parameters[key])
            except:
                continue

        import scipy.stats
        cl = float(self.config['UpperLimit']['cl'])
        delta = 0.5*scipy.stats.chi2.isf(1-2*(cl-0.5), 1)


        if self.config['UpperLimit']['Method'] == "Profile": #The method is Profile
            if Fit.Ts(self.obs.srcname)<2 :
                self.warning("TS of the source is very low, better to use Integral method")
            try:
                import UpperLimits
                ulobject = UpperLimits.UpperLimits(Fit)
                ul, _ = ulobject[self.obs.srcname].compute(emin=self.obs.Emin,
                                          emax=self.obs.Emax,delta=delta)
                                          #delta=2.71 / 2)
                self.info("Upper limit using Profile method: ")
                #print ulobject[self.obs.srcname].results
                self.warning("Be sure to have enough photons to validate the gaussian assumption")
            except RuntimeError:
                self.warning("ST UpperLimits returned RuntimeError, trying Integral")
                self.config['UpperLimit']['Method'] = 'Integral'
        if self.config['UpperLimit']['Method'] == "Integral": #The method is Integral
            import IntegralUpperLimit
            try:
                ul, _ = IntegralUpperLimit.calc_int(Fit, self.obs.srcname, cl=cl,
                                                    verbosity=0,emin=self.obs.Emin,
                                                    emax=self.obs.Emax)
            except RuntimeError:
                self.warning("ST UpperLimits returned RuntimeError, trying Poisson")
                self.config['UpperLimit']['Method'] = 'Poisson'
                
            print(("Upper limit using Integral method: ", ul))
            self.warning("Be sure to have enough photons to validate the gaussian assumption")

        if self.config['UpperLimit']['Method'] == "Poisson": #The method is Poisson
            ul = self.PoissonUL(Fit)
            print(("Upper limit using Poisson statistic: ", ul))

        print("This is an ul on the integral flux in ph/cm2/s")
        return ul #Return the result. This is an ul on the integral flux in ph/cm2/s

    def EnvelopeUL(self, Fit):
        """Compute the envelope UL. An UL is computed for different index and the maximum is taken at each energy.
        This is usefull when the source index is not know or can not be constrain by theoritical argument
        The index range form 1.5 to 2.5"""
        import IntegralUpperLimit
        import UpperLimits
        self._log('EnvelopeUL', 'Compute upper limit envelope')
        PhIndex = Fit.par_index(self.obs.srcname, 'Index')
        Nbp = 20 #Make Nbp computations
        Npgraph = 100#The graph has Npgraph points
        ener = np.logspace(np.log10(self.obs.Emin),
                           np.log10(self.obs.Emax), Npgraph)#the array containing the energy
        Ulenv = np.array(Npgraph * [0.])#the array containing the UL value

        for i in range(Nbp):
            indx = -1.5 - i / (Nbp - 1.)
            utils.FreezeParams(Fit,self.srcname,PhIndex,indx)
            #Use either the profile or the integral method
            self.info("Method used: "+self.config['UpperLimit']['Method'])
            if self.config['UpperLimit']['Method'] == "Profile":
                ul = UpperLimits.UpperLimits(Fit)
                source_ul = ul[self.obs.srcname]
                ul_val, _ = source_ul.compute(emin=self.obs.Emin,
                                              emax=self.obs.Emax,
                                              delta=2.71 / 2)
            if self.config['UpperLimit']['Method'] == "Integral":
                ul_val, _ = IntegralUpperLimit.calc_int(Fit, self.obs.srcname,
                                                        verbosity=0)
            self.success("Upper Limits calculated")
            print(("Index = ", indx, " UL = ", ul_val))  #small print
            for j in range(Npgraph):
                model_name = Fit.model.srcs[self.obs.srcname].spectrum().genericName()
                #compute the DNDE value. The computation change is
                #the model is PowerLaw or PowerLaw2
                #Note : Other model are not taken into account
                #and no UL will be computed
                if model_name == 'PowerLaw2':
                    newUl = ul_val * (indx + 1) * pow(ener[j], indx + 2) \
                        / (pow(self.obs.Emax, indx + 1) - pow(self.obs.Emin, indx + 1))
                elif model_name == 'PowerLaw':
                    IdEScale = utils.getParamIndx(Fit, self.obs.srcname, 'Scale')
                    Escale = Fit[IdEScale].value()
                    newUl = ul_val * pow(ener[j] / Escale, indx + 2)*Escale**2*1.6022e-6
                Ulenv[j] = max(Ulenv[j], newUl)

        print()
        self.info("Result of the UL envelope")
        for j in range(Npgraph):
            print((ener[j], " ", Ulenv[j]))

    def ComputeSED(self, Fit, dump=False):
        """compute the SED with the butterfly for all the model and save it into an ascii file"""
        self._log('PlotSED', 'Generate SED plot')
        from enrico import plotting #plotting is the dedicated library
        from enrico.constants import SpectrumPath
        # filename = self.config['out'] + '/'+SpectrumPath+\
        #         '/SED_' + self.obs.srcname +'_'+ self.config['target']['spectrum']
        filename = utils._SpecFileName(self.config)
        Param = plotting.Params(self.obs.srcname, Emin=self.obs.Emin,
                              Emax=self.obs.Emax, PlotName=filename)
        result = plotting.Result(Fit, Param)
        result.GetDecorrelationEnergy(Param)

        if self.config['Spectrum']['WriteCovMatrix'] == 'yes':
            result._WriteCovMatrix(Param)

        if (dump):
            result._DumpSED(Param)

        return(result)

    def ModelMap(self, xml=None):
        """Make a model Map. Valid only if the statistic is binned"""
        if self.config['analysis']['likelihood'] == 'binned':
            self._log('gtmodel', 'Make model map')#run gtmodel
            if xml==None:
                try:
                    xml = self.outXml 
                except AttributeError:
                    xml = utils._dump_xml(self.config)

            self.obs.ModelMap(xml)
