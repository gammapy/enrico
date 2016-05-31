import os
import scipy.stats
import SummedLikelihood
try:
    import cPickle as pickle
except:
    print("cPickle not found, using pickle instead")
from enrico.fitmaker import FitMaker
from enrico.gtfunction import Observation
from enrico import Loggin

class ModelTester(Loggin.Message):
    """Class to est several models to check
        which is statistically prefered."""
    def __init__(self, config):
        super(ModelTester,self).__init__()
        Loggin.Message.__init__(self)
        self.config = config
        self.folder = self.config['out']
        os.system("mkdir -p "+self.folder+"/TestModel")
        self.modellist = ["PowerLaw","LogParabola","PLSuperExpCutoff"]
        
        '''
        try:
            with open(self.folder+"/TestModel/Fit.pickle","r") as pfile:
                print("Retrieving previous Fit from %s" \
                    %(self.folder+"/TestModel/Fit.pickle"))
                self.FitRunner = pickle.load(pfile)
                self.Fit = self.FitRunner.CreateLikeObject()
        except:
            self._GenFit()
            self.FitRunner.PerformFit(self.Fit, False)

            with open(self.folder+"/TestModel/Fit.pickle","w") as pfile:
                print("Saving current Fit to %s" \
                    %(self.folder+"/TestModel/Fit.pickle"))
                pickle.dump(self.FitRunner,pfile)
        '''
        self._GenFit()
        self.FitRunner.PerformFit(self.Fit, False)

        # Store the results in a dictionnary
        self.Results = {}

    def _GenFit(self):
         try :
             del self.Fit
         except :
             pass

         if self.config['Spectrum']['SummedLike'] == 'yes':
             Obs1 = Observation(self.folder, self.config, tag="FRONT")
             Obs2 = Observation(self.folder, self.config, tag="BACK")
             self.FitRunnerfront = FitMaker(Obs1, self.config)
             self.FitRunnerback = FitMaker(Obs2, self.config)
             self.FitRunnerfront.CreateLikeObject()
             self.FitRunnerback.CreateLikeObject()
             self.Fit = SummedLikelihood.SummedLikelihood()
         else:
             Obs = Observation(self.folder, self.config, tag="")
             self.FitRunner = FitMaker(Obs, self.config)##Class
             self.Fit = self.FitRunner.CreateLikeObject()

    def _printResults(self):
        print 
        self.info("Summary of the results")
        for key in self.modellist:
           if key == "PowerLaw":
               print key," Log(Like) = ",self.Results[key]
               llpl = self.Results[key]
           else :
               TS = 2*(self.Results[key]-llpl)
               prob = 1 - scipy.stats.chi2.cdf(TS, 1)
               print key," Log(Like) = ",self.Results[key]," TS = ",TS," Pvalue = ",prob

    def TestModel(self):
        """ actually test the models """
        Dumpfile = open(self.folder+"/TestModel/TestModel.results","w")
        for key in self.modellist:
            self.Results[key] = self.RunAFit(self.config["target"]["name"],key)
            Dumpfile.write(key + '\t' + str(self.Results[key]) + '\n')
        Dumpfile.close()
        self._printResults()

    def TestModelFromFile(self,inputfile):
        """ Set model and pars from file (test only a custom model)
            This function allow us to test several custom models by calculating 
            their likelihood """
        
        Inputfile = open(inputfile,'r')
        with open(self.folder+"/TestModel/TestModel.results","w") as empty: pass
        Dumpfile  = open(self.folder+"/TestModel/TestModel.results","a+")
        
        for currentmodel in Inputfile.readlines():
            content = currentmodel.rstrip().split(",")
            model = content[0]
            pars = content[1:]

            if model not in self.modellist:
                print("WARNING: given model %s not in the valid range %s" \
                    %(str(model), str(self.modellist)))
                continue

            for k in xrange(len(pars)):
                try: pars[k] = float(pars[k])
                except: pars[k]=None

            # Reduce the list of possible models to the current one
            self.modellist = [model]

            print("Using model %s with parameters %s" %(str(model), str(pars)))
            self.Results[model] = self.RunAFit(self.config["target"]["name"],model,pars)
            _sep_= ', '
            Dumpfile.write(model + _sep_ + _sep_.join([str(k) for k in pars]) \
                + _sep_ + str(self.Results[model]) + '\n')
            
            print("%s Log(Like) = %s" %(model,self.Results[model]))
        
        Inputfile.close()
        Dumpfile.close()

    def RunAFit(self,srcname,model,pars=None):
        # Compute the loglike for the current model for the given parameter set.
        self.info("Computing loglike value for "+model)
        self.Fit.logLike.getSource(srcname).setSpectrum(model)
        
        if (pars==None):
            pars=[None,None,None,None]
        else:
            # Release diffuse models (background)
            self.Fit.thaw(self.Fit.par_index("IsoDiffModel", 'Normalization'))
            self.Fit.thaw(self.Fit.par_index("GalDiffModel", 'Value'))

        if model=="PowerLaw":
            self._setPowerLaw(srcname,pars)
        if model=="LogParabola":
            self._setLogParabola(srcname,pars)
        if model=="PLSuperExpCutoff":
            self._setPLSuperExpCutoff(srcname,pars)
        
        #change the fit tolerance to the one given by the user
        #self.Fit.ftol = float(self.config['fitting']['ftol'])
        try :
          self.Fit.fit(0,covar=True,optimizer=self.config["fitting"]["optimizer"])
          spectrum = self.Fit[self.config['target']['name']].funcs['Spectrum']
          # Get the names of the parameters for the source of interest
          print "Loglike Value for ",model,": ",self.Fit.logLike.value()
          self.Fit.writeXml(self.folder+"/TestModel/TestModel"+model+".xml")
          return self.Fit.logLike.value()
        except :
          self.warning("No convergence for model : "+model+" ??")
          return 0

    def _setPowerLaw(self,name,pars=None):
        SrcSpectrum = self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum']
        SrcSpectrum.getParam('Prefactor').setBounds(1e-7,1e7)
        SrcSpectrum.getParam('Prefactor').setScale(1e-11)
        SrcSpectrum.getParam('Prefactor').setValue(1.)
        SrcSpectrum.getParam('Prefactor').setFree(1)

        SrcSpectrum.getParam('Index').setBounds(-5,0)
        SrcSpectrum.getParam('Index').setValue(-2)
        SrcSpectrum.getParam('Index').setFree(1)

        SrcSpectrum.getParam('Scale').setValue(self.config['energy']['emin'])
        SrcSpectrum.getParam('Scale').setBounds(20,3e6)
        
        # Set each non-None parameter to the wanted value and fix it.
        if pars[0]!=None:
            print("Fixing Prefactor")
            SrcSpectrum.getParam('Prefactor').setFree(0)
            SrcSpectrum.getParam('Prefactor').setValue(pars[0]/1e-11)
            par = self.Fit.par_index(name, 'Prefactor')
            self.Fit.freeze(par)
        if pars[1]!=None:
            print("Fixing Index")
            SrcSpectrum.getParam('Index').setScale(pars[1])
            SrcSpectrum.getParam('Index').setFree(0)
            par = self.Fit.par_index(name, 'Prefactor')
            self.Fit-freeze(par)


    def _setLogParabola(self,name,pars=None):
        SrcSpectrum = self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum']
        SrcSpectrum.getParam('norm').setBounds(1e-7,1e7)
        SrcSpectrum.getParam('norm').setScale(1e-11)
        SrcSpectrum.getParam('norm').setValue(1.)
        SrcSpectrum.getParam('norm').setFree(1)

        SrcSpectrum.getParam('alpha').setBounds(0,5)
        SrcSpectrum.getParam('alpha').setValue(2)
        SrcSpectrum.getParam('alpha').setFree(1)

        SrcSpectrum.getParam('beta').setBounds(0.01,10)
        SrcSpectrum.getParam('beta').setValue(0.5)
        SrcSpectrum.getParam('beta').setFree(1)

        SrcSpectrum.getParam('Eb').setValue(self.config['energy']['emin'])
        SrcSpectrum.getParam('Eb').setFree(0)
        SrcSpectrum.getParam('Eb').setBounds(20,3e6) 
        
        # Set each non-None parameter to the wanted value and fix it.
        if pars[0]!=None:
            print("Fixing norm")
            SrcSpectrum.getParam('norm').setFree(0)
            SrcSpectrum.getParam('norm').setValue(pars[0]/1e-11)
            par = self.Fit.par_index(name, 'norm')
            self.Fit.freeze(par)
        if pars[1]!=None:
            print("Fixing alpha")
            SrcSpectrum.getParam('alpha').setFree(0)
            SrcSpectrum.getParam('alpha').setScale(pars[1])
            par = self.Fit.par_index(name, 'alpha')
            self.Fit.freeze(par)
        if pars[2]!=None:
            print("Fixing beta")
            SrcSpectrum.getParam('beta').setFree(0)
            SrcSpectrum.getParam('beta').setScale(pars[2])
            par = self.Fit.par_index(name, 'beta')
            self.Fit.freeze(par)



    def _setPLSuperExpCutoff(self,name,pars=None):
        SrcSpectrum = self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum']
        SrcSpectrum.getParam('Prefactor').setBounds(1e-7,1e7)
        SrcSpectrum.getParam('Prefactor').setScale(1e-11)
        SrcSpectrum.getParam('Prefactor').setValue(1.)
        SrcSpectrum.getParam('Prefactor').setFree(1)

        SrcSpectrum.getParam('Index1').setBounds(-5,0)
        SrcSpectrum.getParam('Index1').setValue(-2)
        SrcSpectrum.getParam('Index1').setFree(1)

        SrcSpectrum.getParam('Index2').setValue(-1)
        SrcSpectrum.getParam('Index2').setBounds(-5,-0.05)
        SrcSpectrum.getParam('Index2').setFree(0)

        SrcSpectrum.getParam('Cutoff').setBounds(20,3e6)
        SrcSpectrum.getParam('Cutoff').setValue(1e4)
        SrcSpectrum.getParam('Cutoff').setFree(1)

        SrcSpectrum.getParam('Scale').setValue(self.config['energy']['emin'])
        SrcSpectrum.getParam('Scale').setBounds(20,3e6)
        
        # Set each non-None parameter to the wanted value and fix it.
        if pars[0]!=None:
            print("Fixing Prefactor")
            SrcSpectrum.getParam('Prefactor').setFree(0)
            SrcSpectrum.getParam('Prefactor').setValue(pars[0]/1e-11)
            par = self.Fit.par_index(name, 'Prefactor')
            self.Fit.freeze(par)
        if pars[1]!=None:
            print("Fixing Index1")
            SrcSpectrum.getParam('Index1').setScale(pars[1])
            SrcSpectrum.getParam('Index1').setFree(0)
            par = self.Fit.par_index(name, 'Index1')
            self.Fit.freeze(par)
        if pars[2]!=None:
            print("Fixing Index2")
            SrcSpectrum.getParam('Index2').setScale(pars[2])
            SrcSpectrum.getParam('Index2').setFree(0)
            par = self.Fit.par_index(name, 'Index2')
            self.Fit.freeze(par)
        if pars[3]!=None:
            print("Fixing Cutoff")
            SrcSpectrum.getParam('Cutoff').setScale(pars[3])
            SrcSpectrum.getParam('Cutoff').setFree(0)
            par = self.Fit.par_index(name, 'Cutoff')
            self.Fit.freeze(par)

