import os
import scipy.stats
import SummedLikelihood
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
        """ Set model and pars from file (test only a custom model) """
        with open(inputfile,'r') as r:
            content = r.read().rstrip().split(",")
            model = content[0]
            pars = content[1:]
               
            Dumpfile = open(self.folder+"/TestModel/TestModel.results","w")
            if model not in self.modellist:
                print("WARNING: given model %s not in the valid range %s" %(str(model), str(self.modellist)))
                for key in self.modellist:
                    self.Results[key] = self.RunAFit(self.config["target"]["name"],key)
                    Dumpfile.write(key + '\t' + str(self.Results[key]) + '\n')
            else:
                for k in xrange(len(pars)):
                    try: pars[k] = float(pars[k])
                    except: pars[k]=None
                
                # Reduce the list of possible models to the current one
                self.modellist = [model]
                
                print("Using model %s with parameters %s" %(str(model), str(pars)))
                self.Results[model] = self.RunAFit(self.config["target"]["name"],model,pars)
                Dumpfile.write(model + '\t' + str(self.Results[model]) + '\n')
                Dumpfile.close()
                print("%s Log(Like) = %s" %(model,self.Results[model]))
                #self._printResults()

    def RunAFit(self,srcname,model,pars=None):
        self.info("Computing loglike value for "+model)
#        self._GenFit()
        self.Fit.logLike.getSource(srcname).setSpectrum(model)
        
        if model=="PowerLaw":
            self._setPowerLaw(srcname,pars)
        if model=="LogParabola":
            self._setLogParabola(srcname,pars)
        if model=="PLSuperExpCutoff":
            self._setPLSuperExpCutoff(srcname,pars)

        #change the fit tolerance to the one given by the user
#        self.Fit.ftol = float(self.config['fitting']['ftol'])
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

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-7,1e7)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-11)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(1.)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setBounds(-5,0)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setValue(-2)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setFree(1)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setValue(self.config['energy']['emin'])
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(20,3e6)
        
        if pars!=None:
            if pars[0]!=None:
                print("Fixing Prefactor")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-11)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(pars[0]/1e-11)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(pars[0],pars[0])
            if pars[1]!=None:
                print("Fixing Index")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setScale(pars[1])
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setFree(0)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setBounds(pars[1],pars[1])


    def _setLogParabola(self,name,pars=None):
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('norm').setBounds(1e-7,1e7)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('norm').setScale(1e-11)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('norm').setValue(1.)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('alpha').setBounds(0,5)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('alpha').setValue(2)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('alpha').setFree(1)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('beta').setBounds(0.01,10)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('beta').setValue(0.5)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('beta').setFree(1)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Eb').setValue(self.config['energy']['emin'])
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Eb').setFree(0)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Eb').setBounds(20,3e6)
        
        if pars!=None:
            if pars[0]!=None:
                print("Fixing norm")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('norm').setScale(1e-11)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('norm').setValue(pars[0]/1e-11)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('norm').setBounds(pars[0]/1e-11,pars[0]/1e-11)
            if pars[1]!=None:
                print("Fixing alpha")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('alpha').setScale(pars[1])
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('alpha').setFree(0)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('alpha').setBounds(pars[1],pars[1])
            if pars[2]!=None:
                print("Fixing beta")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('beta').setScale(pars[2])
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('beta').setFree(0)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('beta').setBounds(pars[2],pars[2])



    def _setPLSuperExpCutoff(self,name,pars=None):
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-7,1e7)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-11)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(1.)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index1').setBounds(-5,0)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index1').setValue(-2)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index1').setFree(1)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index2').setValue(-1)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index2').setBounds(-5,-0.05)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index2').setFree(0)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Cutoff').setBounds(20,3e6)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Cutoff').setValue(1e4)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Cutoff').setFree(1)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setValue(self.config['energy']['emin'])
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(20,3e6)
        
        if pars!=None:
            if pars[0]!=None:
                print("Fixing Prefactor")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-11)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(pars[0]/1e-11)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(pars[0]/1e-11,pars[0]/1e-11)
            if pars[1]!=None:
                print("Fixing Index1")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index1').setScale(pars[1])
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index1').setFree(0)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index1').setBounds(pars[1],pars[1])
            if pars[2]!=None:
                print("Fixing Index2")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index2').setScale(pars[2])
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index2').setFree(0)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index2').setBounds(pars[2],pars[2])
            if pars[3]!=None:
                print("Fixing Cutoff")
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Cutoff').setScale(pars[3])
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Cutoff').setFree(0)
                self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Cutoff').setBounds(pars[3],pars[3])

