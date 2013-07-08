# fitmake.py written by David Sanchez : david.sanchez@mpi-hd.mpg.de
# Collection of functions to run the fit (gtlike)
# the class Makefit will call the function of the observation class (gtfunction.py) and prepare the fit
# by computing the fits file.
# it can distinguish between the binned and unbinnned analysis
# begun September 2011

import numpy as np
from UnbinnedAnalysis import UnbinnedAnalysis, UnbinnedObs
from BinnedAnalysis import BinnedAnalysis, BinnedObs
import utils
import os
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

class FitMaker(object):
    """Collection of functions to prepare/run the GTLIKE fit
     and compute an upper limit is needed"""
    def __init__(self, obs, config):
        self.obs = obs
        self.config = config
        self.task_number = 1
        self.log_like = 0

    def _log(self, task='', description=''):
        print
        print('# ' + '*' * 60)
        if task:
            task = '%10s --- ' % task
        print('# *** %3d %s%s' %
              (self.task_number, task, description))
        print '# ' + '*' * 60
        self.task_number += 1

    def GenerateFits(self):
        """Run the different ST tools and compute the fits files
           First it runs the tools that are common to the binned 
           and unbinned analysis chain then it run the specific
           tools following the choise of the user"""

        #Run the tools common to binned and unbinned chain
        self._log('gtselect', 'Select data from library')#run gtselect
        self.obs.FirstCut()
        self._log('gtmktime', 'Update the GTI and cut data based on ROI')#run gtdiffresp
        self.obs.MkTime()
        if self.config["analysis"]["ComputeDiffrsp"] == "yes":
            self._log('gtdiffrsp', 'Compute Diffuse response')
            self.obs.DiffResps()#run gtbin
        self._log('gtbin', 'Create a count map')
        self.obs.Gtbin()
        self._log('gtltcube', 'Make live time cube')#run gtexpcube
        self.obs.ExpCube()

        #Choose between the binned of the unbinned analysis
        if self.config['analysis']['likelihood'] == 'binned': #binned analysis chain
            self._log('gtbin', 'Make count map CCUBE')#run gtbin
            self.obs.GtCcube()
            self._log('gtexpcube2', 'Make binned exposure cube')#run gtexpcube2
            self.obs.GtBinnedMap()
            self._log('gtsrcmap', 'Make a source map')#run gtsrcmap
            self.obs.SrcMap()

        if self.config['analysis']['likelihood'] == 'unbinned': #unbinned analysis chain
            self._log('gtexpmap', 'Make an exposure map')
            self.obs.ExpMap()
    #the function ends here. It does not run gtlike

    def CreateLikeObject(self):
        """Create an UnbinnedAnalysis or a BinnedAnalysis and retrun it.
        By default, the optimizer is DRMNGB """

        #create binnedAnalysis object
        if self.config['analysis']['likelihood'] == 'binned':
            Obs = BinnedObs(srcMaps=self.obs.scrMap,
                            expCube=self.obs.Cubename,
                            binnedExpMap=self.obs.BinnedMapfile,
                            irfs=self.obs.irfs)
            Fit = BinnedAnalysis(Obs, self.obs.xmlfile,
                                 optimizer='DRMNGB')

        #create a unbinnedAnalysis object
        if self.config['analysis']['likelihood'] == 'unbinned':
            Obs = UnbinnedObs(self.obs.eventfile,
                              self.obs.ft2,
                              expMap=self.obs.Mapname,
                              expCube=self.obs.Cubename,
                              irfs=self.obs.irfs)
            Fit = UnbinnedAnalysis(Obs, self.obs.xmlfile,
                                   optimizer='DRMNGB')

        if float(self.config['Spectrum']['FrozenSpectralIndex']) > 0:
            if Fit.model.srcs[self.obs.srcname].spectrum().genericName()=="PowerLaw" or Fit.model.srcs[self.obs.srcname].spectrum().genericName()=="PowerLaw2":
                PhIndex = Fit.par_index(self.obs.srcname, 'Index')
                Fit[PhIndex] = -float(self.config['Spectrum']['FrozenSpectralIndex'])
                Fit.freeze(PhIndex)
                print "Freezing spectral index at ",-float(self.config['Spectrum']['FrozenSpectralIndex'])
            else:
              log.warning("The model is not a PowerLaw. Cannot freeze the index.")
        return Fit #return the BinnedAnalysis or UnbinnedAnalysis object.

    def PerformFit(self, Fit):
        """Run gtlile tool. First it run gtlike with the DRNMGB optimizer
        and the user optimizer after. A dictionnay is return with all
        the releveant results"""

        self._log('gtlike', 'Run likelihood analysis')
        try:
            Fit.fit(0) #first try to run gtlike to approche the minimum
        except:
            pass

        # Now the precise fit will be done
        #change the fit tolerance to the one given by the user
        Fit.ftol = float(self.config['fitting']['ftol'])
        #fit with the user optimizer and ask gtlike to compute the covariance matrix 
        self.log_like = Fit.fit(0,covar=True, optimizer=self.config['fitting']['optimizer']) #fit with the user optimizer and ask gtlike to compute the covariance matrix 

        self.RemoveWeakSources(Fit)#remove source with TS<1 to be sure that MINUIT will converge

        Fit.writeXml(utils._dump_xml(self.config))

    def RemoveWeakSources(self,Fit):
        """Remove the weak source after a fit and reoptimized
         weak mens TS<1"""
        self._log('','Remove all the weak (TS<1) sources')
        NoWeakSrcLeft = False
        while not(NoWeakSrcLeft):
            NoWeakSrcLeft = True
            for src in Fit.model.srcNames:
                ts = Fit.Ts(src)
                if  ts< 1 and not(src == self.obs.srcname) and Fit.logLike.getSource(src).getType() == 'Point':
                    print "delete source : ", src," with TS = ",ts
                    NoWeakSrcLeft = False
                    Fit.deleteSource(src)
            if not(NoWeakSrcLeft):
                self._log('Re-optimize', False)
                Fit.fit(0,covar=True, optimizer=self.config['fitting']['optimizer'])
            print
        return Fit

    def GetAndPrintResults(self, Fit):
        """Get and print some useful results. Also contruct a dictonnary and fill it with results"""
        if self.config['verbose'] == 'yes' :
            self._log('Results', 'Print results of the fit')
        Result = {}

        if self.config['verbose'] == 'yes' :
            print Fit.model
            print
            # Print src name, Npred and TS for source with TS > 5
            print "Source Name\tNpred\tTS"
            for src in Fit.model.srcNames:
                if Fit.Ts(src) > 5:
                    print src, "\t%2.3f\t%2.3f" % (Fit.NpredValue(src), Fit.Ts(src))
            print '\n# ' + '*' * 60 +'\n'

        # fill the dictonnary with some values
        Result['Optimizer'] = self.config['fitting']['optimizer']
        Result['Npred'] = Fit.NpredValue(self.obs.srcname)
        Result['TS'] = Fit.Ts(self.obs.srcname)
        if self.config['verbose'] == 'yes' :
            print "Values and (MINOS) errors for " + self.obs.srcname
            print "TS : ", Fit.Ts(self.obs.srcname)

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

        for par in ParName : #Loop over the parameters and get value, error and scale
            ParValue = spectrum.getParam(par).value()
            ParError = spectrum.getParam(par).error()
            Scale = spectrum.getParam(par).getScale()
            Result[par] = ParValue * Scale
            Result['d'+par] = ParError * Scale
            if ParError>0: # Compute MINOS errors for relevent parameters  Fit.Ts(self.obs.srcname) > 5 and
                try:
                    MinosErrors = Fit.minosError(self.obs.srcname, par)
                    if self.config['verbose'] == 'yes' :
                       print(par+" :  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ] %2.0e" %
                          (ParValue, ParError, MinosErrors[0], MinosErrors[1], Scale))
                       Result.update({'d'+par+'-': MinosErrors[0] * Scale})
                       Result.update({'d'+par+'+': MinosErrors[1] * Scale})
                except:
                    if self.config['verbose'] == 'yes' :
                        print(par+" :  %2.2f +/-  %2.2f  %2.0e" %
                          (ParValue, ParError, Scale))
            else:
                if self.config['verbose'] == 'yes' :
                    print(par+" :  %2.2f +/-  %2.2f  %2.0e" %
                      (ParValue, ParError, Scale))

        try: # get covariance matrix
            if self.config['verbose'] == 'yes' :
                utils.GetCovar(self.obs.srcname, Fit)
        except:
            pass #if the covariance matrix has not been computed

        #Compute an UL if the source is too faint
        if float(self.config['UpperLimit']['TSlimit']) > Fit.Ts(self.obs.srcname):
            if self.config['UpperLimit']['envelope'] == 'yes':
                self.EnvelopeUL(Fit)
            else:
                Ulval = self.ComputeUL(Fit)
                Result['Ulvalue'] = Ulval

        return Result   #Return the dictionnary

    def ComputeUL(self, Fit):
        """Compute an Upper Limit using either the profil or integral method
        See the ST cicerone for more information on the 2 method"""

        self._log('UpperLimit', 'Compute upper Limit')
        #Index given by the user  
        print "Assumed index is ", self.config['UpperLimit']['SpectralIndex']

        IdGamma = utils.getParamIndx(Fit, self.obs.srcname, 'Index')
        Fit[IdGamma] = -self.config['UpperLimit']['SpectralIndex']#set the index
        Fit[IdGamma].setFree(0)#the variable index is frozen to compute the UL

        import scipy.stats
        cl = float(self.config['UpperLimit']['cl'])
        delta = 0.5*scipy.stats.chi2.isf(1-2*(cl-0.5), 1)
        print cl,' ',delta

        if self.config['UpperLimit']['Method'] == "Profile": #The method is Profile
            import UpperLimits
            ulobject = UpperLimits.UpperLimits(Fit)
            ul, _ = ulobject[self.obs.srcname].compute(emin=self.obs.Emin,
                                      emax=self.obs.Emax,delta=delta)
                                      #delta=2.71 / 2)
            print "Upper limit using Profile method: "
            print ulobject[self.obs.srcname].results
        if self.config['UpperLimit']['Method'] == "Integral": #The method is Integral
            import IntegralUpperLimit
            ul, _ = IntegralUpperLimit.calc_int(Fit, self.obs.srcname, cl=cl,
                                                verbosity=0)
            print "Upper limit using Integral method: ", ul
        return ul #Return the result. This is an ul on the integral flux in ph/cm2/s 

    def EnvelopeUL(self, Fit):
        """Compute the envelope UL. An UL is computed for different index and the maximum is taken at each energy.
        This is usefull when the source index is not know or can not be constrain by theoritical argument
        The index range form 1.5 to 2.5"""
        import IntegralUpperLimit
        self._log('EnvelopeUL', 'Compute upper limit envelope')
        PhIndex = Fit.par_index(self.obs.srcname, 'Index')
        Nbp = 20 #Make Nbp computations
        Npgraph = 100#The graph has Npgraph points
        ener = np.logspace(np.log10(self.obs.Emin),
                           np.log10(self.obs.Emax), Npgraph)#the array containing the energy
        Ulenv = np.array(Npgraph * [0.])#the array containing the UL value

        for i in xrange(Nbp):
            indx = -1.5 - i / (Nbp - 1.)
            Fit[PhIndex] = indx
            Fit.freeze(PhIndex)#Freeze the index
            #Use either the profile or the integral method
            if self.config['UpperLimit']['Method'] == "Profile":
                ul = UpperLimits.UpperLimits(Fit)
                source_ul = ul[self.obs.srcname]
                ul_val, _ = source_ul.compute(emin=self.obs.Emin,
                                              emax=self.obs.Emax,
                                              delta=2.71 / 2)
            if self.config['UpperLimit']['Method'] == "Integral":
                ul_val, _ = IntegralUpperLimit.calc_int(Fit, self.obs.srcname,
                                                        verbosity=0)
            print "Index = ", indx, " UL = ", ul_val  #small print
            for j in xrange(Npgraph):
                model_name = Fit.model.srcs[self.obs.srcname].spectrum().genericName()
                #compute the DNDE value. The computation change is 
                #the model is PowerLaw or PowerLaw2
                #Note : Other model are not taken into account 
                #and no UL will be computed
                if model_name == 'PowerLaw2':
                    newUl = ul_val * (indx + 1) * pow(ener[j], indx + 2) / (pow(self.obs.Emax, indx + 1) - pow(self.obs.Emin, indx + 1))
                elif model_name == 'PowerLaw':
                    IdEScale = utils.getParamIndx(Fit, self.obs.srcname, 'Scale')
                    Escale = Fit[IdEScale].value()
                    newUl = ul_val * pow(ener[j] / Escale, indx + 2)
                Ulenv[j] = max(Ulenv[j], newUl)

        print
        print "Result of the UL envelope"
        for j in xrange(Npgraph):
            print ener[j], " ", Ulenv[j]

    def ComputeSED(self, Fit):
        """compute the SED with the butterfly for all the model and save it into an ascii file"""
        self._log('PlotSED', 'Generate SED plot')
        import plotting#plotting is the dedicated library
        filename = self.config['out'] + '/Spectrum/SED_' + self.obs.srcname +'_'+ Fit[self.obs.srcname].funcs['Spectrum'].genericName()
        Param = plotting.Params(self.obs.srcname, Emin=self.obs.Emin, 
                              Emax=self.obs.Emax, PlotName=filename)
        result = plotting.Result(Fit, Param)
        result._DumpSED(Param)

    def ModelMap(self, xml):
        """Make a model Map. Valid only if the statistic is binned"""
        if self.config['analysis']['likelihood'] == 'binned':
            self._log('gtmodel', 'Make model map')#run gtmodel
            self.obs.ModelMaps(xml)
