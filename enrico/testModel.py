from fitmaker import FitMaker
from gtfunction import Observation
from enrico.config import get_config
import os

class ModelTester:
    """Class to est several models to check
        which is statistically prefered."""
    def __init__(self, infile):
         self.config = get_config(infile)
         self.folder = self.config['out']
         os.system("mkdir -p "+self.folder+"/TestModel")
         convtype = self.config['analysis']['convtype']

         Obs = Observation(self.folder, self.config, convtype, tag="")
         FitRunner = FitMaker(Obs, self.config)##Class
         self.Fit = FitRunner.CreateLikeObject()
         # Store the results in a dictionnary
         self.Results = {}
         self.Results["PowerLaw"] = 0
         self.Results["LogParabola"] = 0
         self.Results["PLSuperExpCutoff"] = 0

    def TestModel(self):
        """ actually test the models """
        Dumpfile = open(self.folder+"/TestModel/TestModel.results","w")
        for key in self.Results.iterkeys():
            self.Results[key] = self.RunAFit(self.config["target"]["name"],key)
            Dumpfile.write(key + '\t' + str(self.Results[key]) + '\n')
        Dumpfile.close()

    def RunAFit(self,srcname,model):
        print "Computing loglike value for: ",model
        self.Fit.logLike.getSource(srcname).setSpectrum(model)
        if model=="PowerLaw":
            self._setPowerLaw(srcname)
        if model=="LogParabola":
            self._setLogParabola(srcname)
        if model=="PLSuperExpCutoff":
            self._setPLSuperExpCutoff(srcname)

        self.Fit.fit(0,optimizer=self.config["fitting"]["optimizer"])
        return self.Fit.logLike.value()


    def _setPowerLaw(self,name):
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-7,1e7)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-11)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(1.)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setBounds(-5,0)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setValue(-2)
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setFree(1)

        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setValue(self.config['energy']['emin'])
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(20,3e6)

    def _setLogParabola(self,name):
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
        self.Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Eb').setBounds(20,3e6)



    def _setPLSuperExpCutoff(self,name):
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

