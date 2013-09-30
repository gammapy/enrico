import os
from math import sqrt
import numpy as np
import ROOT
from enrico import utils
from enrico import root_style
from enrico import plotting
from enrico import environ
from enrico.config import get_config
from enrico.constants import LightcurvePath,FoldedLCPath
from enrico.submit import call
from enrico.RunGTlike import run, GenAnalysisObjects

class LightCurve:
    """Class to calculate light curves and variability indexes."""
    def __init__(self, config):
        ROOT.gROOT.SetBatch(ROOT.kTRUE) #Batch mode
        self.config = config

        #Read the config
        self.srcname = self.config['target']['name'] #src name
        self.Tag = self.config['file']['tag']
        self.tmin = self.config['time']['tmin']
        self.tmax = self.config['time']['tmax']

        self.submit = self.config['Submit']
        # One point of the LC will be computed as a spectrum plot.
        # enrico_sed will be used
        # Do fits files will be generated
        self.config['Spectrum']['FitsGeneration'] = self.config['LightCurve']['FitsGeneration']
        #Froze the Spectral index at a value of self.config['LightCurve']['SpectralIndex'] (no effect if 0)
        self.config['Spectrum']['FrozenSpectralIndex'] = self.config['LightCurve']['SpectralIndex']
        #TS limit. Compute an UL if the TS is below TSLightCurve
        self.config['UpperLimit']['TSlimit'] = self.config['LightCurve']['TSLightCurve']

        self.folder = self.config['out']

        #No plot, no bin in energy, Normal UL
        self.config['Spectrum']['ResultPlots'] = 'no'
        self.config['Ebin']['NumEnergyBins'] = 0
        self.config['UpperLimit']['envelope'] = 'no'
        #No submition. Submission will be directly handle by this soft
        self.config['Submit'] = 'no'
#        self.config['verbose'] ='no' #Be quiet

        self.configfile = []#All the config file in the disk are stored in a list

    def _MakeTimeBins(self):
        self.time_array = np.zeros(0)
        self.Nbin = 0
        self.gtifile = []
        if self.config['time']['file'] != '':
            print "use ",self.config['time']['file'] 
            self.gtifile.append(self.config['time']['file'])
            times = np.genfromtxt(self.gtifile[0],dtype="float",unpack=True)
            self.Nbin = times.size/2
            self.time_array=np.reshape(times,times.size,'F')
            if self.config['time']['type']=='MJD':
                 self.time_array = utils.MJD_to_met(self.time_array)
            elif config['time']['type']=='JD':
                 self.time_array = utils.JD_to_met(self.time_array)
        else:
            self.Nbin = self.config['LightCurve']['NLCbin']
            self.time_array = np.zeros(self.Nbin*2)
#            self.dt = (self.tmax - self.tmin) / self.Nbin
            t = np.arange(self.tmin,self.tmax+0.000001,(self.tmax - self.tmin) / self.Nbin)
            for i in xrange(self.Nbin):
                self.time_array[2*i] = t[i]
                self.time_array[2*i+1]= t[i+1]

        print "Running LC with ",self.Nbin," bins"
        for i in xrange(self.Nbin):
            print "Bin ",i," Start=",self.time_array[2*i]," Stop=",self.time_array[2*i+1]
        print 

    def _errorReading(self,message,i):
        print "WARNING : "+message+" : ",self.configfile[i]
        print "Job Number : ",i
        print "Please have a look at this job log file"


    def _ManageFolder(self,path):
        """   All files will be stored in a subfolder name path + NLCbin
        Create a subfolder"""
        self.LCfolder =  self.folder+"/"+path+"_"+str(self.Nbin)+"bins/"
        os.system("mkdir -p "+self.LCfolder)
        self.config['out'] = self.LCfolder


    def PrepareLC(self,write = 'no'):
        """Simple function to prepare the LC generation : generate and write the config files"""

        for i in xrange(self.Nbin):
            self.config['time']['tmin'] = self.time_array[2*i]
            self.config['time']['tmax'] = self.time_array[2*i+1]
            self.config['file']['tag'] = self.Tag + '_LC_' + str(i)
            filename = (self.config['out'] + "Config_" + str(i) + "_" +
                    str(self.config['time']['tmin']) + "_" +
                    str(self.config['time']['tmax']))#Name of the config file

            if len(self.gtifile)==1:
                self.config['time']['file']=self.gtifile[0]
            elif len(self.gtifile)>1:
                print 'Time selection file for bin {0} = {1}'.format(i,self.gtifile[i])
                self.config['time']['file']=self.gtifile[i]

            if write == 'yes':
                self.config.write(open(filename, 'w'))

            self.configfile.append(filename)


    def _MakeLC(self,Path=LightcurvePath) :
        '''Main function of the Lightcurve script. Read the config file and run the gtlike analysis'''
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')

        self.PrepareLC(self.config['LightCurve']['MakeConfFile'])#Get the config file

        for i in xrange(self.Nbin):
            if self.submit == 'yes':
                cmd = "enrico_sed "+self.configfile[i]
                scriptname = self.LCfolder+"LC_Script_"+str(i)+".sh"
                JobLog = self.LCfolder+"LC_Job_"+str(i)+".log"
                JobName = (self.config['target']['name'] + "_" +
                       self.config['analysis']['likelihood'] +
                       "_LC_" + self.config['file']['tag'])+"_"+str(i)+".log"

                call(cmd,enricodir,fermidir,scriptname,JobLog,JobName)#Submit the job
            else :
                run(self.configfile[i])#run in command line


    def _MakePhasebin(self):
        """document me """
        self.time_array = np.zeros(self.Nbin*2)
        self.config['time']['type'] = 'MJD'
        T0 = self.config['FoldedLC']['epoch']
        Period = self.config['FoldedLC']['Period']
        t1 = utils.met_to_MJD(self.config['time']['tmin'])
        t2 = utils.met_to_MJD(self.config['time']['tmax'])

        if T0==0:
            T0=t1
        elif t1 < T0:
           T0 -= np.ceil((T0-t1)/Period)*Period
           # find orbit numbers covered by range (t1,t2)
        norbt1 = int(np.floor((t1-T0)/Period))
        norbt2 = int(np.ceil((t2-T0)/Period))

        phase = np.linspace(0,1,self.Nbin+1)

        self.gtifile=[] #reset gtifiles
        for i in range(self.Nbin):
            self.time_array[2*i] = self.config['time']['tmin']
            self.time_array[2*i+1] = self.config['time']['tmax']
            gtifn = os.path.join(self.LCfolder,"TimeSelection_{0:02.0f}.dat".format(i))
            ints=[]

            for norb in range(norbt1,norbt2+1):
                ints.append(((T0+(norb+phase[i])*Period),(T0+(norb+phase[i+1])*Period)))
            tsel = np.array(ints)
            np.savetxt(gtifn,tsel)
            self.gtifile.append(gtifn)


    def MakeLC(self) :
        """Run a std lc """
        self._MakeTimeBins()
        self._ManageFolder(LightcurvePath)
        self._MakeLC()

    def MakeFoldedLC(self):
        """run a folded lc """
        self.Nbin = self.config['FoldedLC']['NLCbin']
        self._ManageFolder(FoldedLCPath)
        self._MakePhasebin()
        self._MakeLC()


    def PlotLC(self):
        '''Plot a lightcurve which have been generated previously'''
        self._MakeTimeBins()
        self._ManageFolder(LightcurvePath)
        self.PrepareLC()#Get the config file
        self._PlotLC()

    def PlotFoldedLC(self):
        '''Plot a lightcurve which have been generated previously'''
        self.Nbin = self.config['FoldedLC']['NLCbin']
        self._ManageFolder(FoldedLCPath)
        self._MakePhasebin()
        self.PrepareLC()#Get the config file
        self._PlotLC(True)


    def _PlotLC(self,folded=False):
        root_style.RootStyle()#Nice plot style

        print "Reading files produced by enrico"
        LcOutPath = self.LCfolder + self.config['target']['name']

        #Result are stored into list. This allow to get rid of the bin which failled
        Time = []
        TimeErr = []
        Flux = []
        FluxErr = []
        FluxForNpred = []
        FluxErrForNpred = []
        Npred = []
        Npred_detected_indices = []
        TS = []

        for i in xrange(self.Nbin):
            CurConfig = get_config(self.configfile[i])
            #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
            try :
                ResultDic = utils.ReadResult(CurConfig)
            except :
                self._errorReading("fail reading config file",i)
                continue

            #Update the time and time error array
            Time.append((ResultDic.get("tmax")+ResultDic.get("tmin"))/2.)
            TimeErr.append((ResultDic.get("tmax")-ResultDic.get("tmin"))/2.)
            #Check is an ul have been computed. The error is set to zero for the TGraph.
            if ResultDic.has_key('Ulvalue') :
                Flux.append(ResultDic.get("Ulvalue"))
                FluxErr.append(0)
            else :
                Flux.append(ResultDic.get("Flux"))
                FluxErr.append(ResultDic.get("dFlux"))
            FluxErrForNpred.append(ResultDic.get("dFlux"))
            FluxForNpred.append(ResultDic.get("Flux"))
            #Get the Npred and TS values
            Npred.append(ResultDic.get("Npred"))
            TS.append(ResultDic.get("TS"))
            if (CurConfig['LightCurve']['TSLightCurve']<float(ResultDic.get("TS"))):
                Npred_detected_indices.append(i)

        #change the list into np array
        TS = np.array(TS)
        Npred = np.array(Npred)
        Npred_detected = Npred[Npred_detected_indices]
        Time = np.array(Time)
        TimeErr = np.array(TimeErr)
        Flux = np.array(Flux)
        FluxErr = np.array(FluxErr)
        FluxForNpred = np.array(FluxForNpred)
        FluxErrForNpred = np.array(FluxErrForNpred)

        #Plots the diagnostic plots is asked
        # Plots are : Npred vs flux
        #             TS vs Time
        if self.config['LightCurve']['DiagnosticPlots'] == 'yes':
            fittedFunc = self.CheckNpred(Npred,FluxForNpred,FluxErrForNpred,Npred_detected_indices)#check the errors calculation
            gTHNpred,TgrNpred = plotting.PlotNpred(Npred,FluxForNpred,FluxErrForNpred)
            CanvNpred = _GetCanvas()
            gTHNpred.Draw()
            TgrNpred.Draw('zP')

            _,TgrNpred_detected = plotting.PlotNpred(Npred_detected,Flux[Npred_detected_indices],FluxErrForNpred[Npred_detected_indices])
            TgrNpred_detected.SetLineColor(2)
            TgrNpred_detected.SetMarkerColor(2)
            TgrNpred_detected.Draw('zP')
            fittedFunc.Draw("SAME")

            CanvNpred.Print(LcOutPath+"_Npred.eps")
            CanvNpred.Print(LcOutPath+"_Npred.C")

            gTHTS,TgrTS = plotting.PlotTS(Time,TimeErr,TS)
            CanvTS = _GetCanvas()
            gTHTS.Draw()
            TgrTS.Draw('zP')
            CanvTS.Print(LcOutPath+'_TS.eps')
            CanvTS.Print(LcOutPath+'_TS.C')


#    Plot the LC itself. This function return a TH2F for a nice plot
#    a TGraph and a list of TArrow for the ULs
        if folded:
            phase = np.linspace(0,1,self.Nbin+1)
            Time = (phase[1:]+phase[:-1])/2.
            TimeErr = (phase[1:]-phase[:-1])/2.
            gTHLC,TgrLC,ArrowLC = plotting.PlotFoldedLC(Time,TimeErr,Flux,FluxErr)
        else :
            gTHLC,TgrLC,ArrowLC = plotting.PlotLC(Time,TimeErr,Flux,FluxErr)
        CanvLC = ROOT.TCanvas()
        gTHLC.Draw()
        TgrLC.Draw('zP')

        #plot the ul as arrow
        for i in xrange(len(ArrowLC)):
            ArrowLC[i].Draw()

        # compute Fvar and probability of being cst
        self.FitWithCst(TgrLC)
        self.Fvar(Flux,FluxErr)

        #Save the canvas in the LightCurve subfolder
        CanvLC.Print(LcOutPath+'_LC.eps')
        CanvLC.Print(LcOutPath+'_LC.C')

        #Dump into ascii
        lcfilename = LcOutPath+"_results.dat"
        print "Write to Ascii file : ",lcfilename
        WriteToAscii(Time,TimeErr,Flux,FluxErr,TS,Npred,lcfilename)

        if self.config["LightCurve"]['ComputeVarIndex'] == 'yes':
             self.VariabilityIndex()

    def CheckNpred(self,Npred,Flux,FluxErr,detected_indices):
        """check if the errors are well computed using the Npred/sqrt(Npred) vs Flux/FluxErr relation
           and print corresponding point which failled"""
        _,TgrNpred = plotting.PlotNpred(Npred[detected_indices],Flux[detected_indices],FluxErr[detected_indices])
        func = ROOT.TF1("func","pol1",np.min(np.sqrt(Npred)),np.max(np.sqrt(Npred)))
        TgrNpred.Fit(func,"Q")
        for i in xrange(len(Flux)):
            if Flux[i]/FluxErr[i]>2*func.Eval(sqrt(Npred[i])):
                self._errorReading("problem in errors calculation for",i)
                print "Flux +/- error = ",Flux[i]," +/- ",FluxErr[i]
                print "V(Npred) = ",sqrt(Npred[i])
                print 
        func.SetLineColor(15)
        func.SetLineStyle(2)
        return func

    def Fvar(self,Flux,FluxErr):
        """Compute the Fvar as defined in Vaughan et al."""
        moy=np.average(Flux)
        var=np.var(Flux)
        expvar=np.average((np.array(FluxErr))**2)
        intvar=var-expvar #Correct for errors
        print 
        try :
            fvar=sqrt(intvar)/moy
            err_fvar = sqrt( ( 1./sqrt(2*len(Flux))*expvar/moy**2/fvar)**2 + (sqrt(expvar/len(Flux))*1./moy)**2)
            print "\t Fvar = ",fvar," +/- ",err_fvar
        except :
            print  "\t Fvar is negative, Fvar**2 = %2.2e +/- %2.2e"%(intvar/(moy*moy), ((1./sqrt(2*len(Flux))*expvar/moy**2)**2/(intvar/(moy*moy)) + (sqrt(expvar/len(Flux))*1./moy)**2))
        print 

    def FitWithCst(self,tgrFlux):
        """Fit the LC with a constant function an
           print the chi2 and proba"""
        func = ROOT.TF1('func','pol0')
        func.SetLineColor(15)
        func.SetLineStyle(3)
        tgrFlux.Fit('func','Q')
        print
        print '\tChi2 = ',func.GetChisquare()," NDF = ",func.GetNDF()
        print '\tprobability of being cst = ',func.GetProb()
        print

    def VariabilityIndex(self):
        """Compute the variability index as in the 2FLG catalogue. (see Nolan et al, 2012)"""
        LcOutPath = self.LCfolder + self.config['target']['name']

        utils._log('Computing Variability index ')

        self.config['Spectrum']['FitsGeneration'] = 'no'
#        ValueDC = self.GetDCValue()
        ResultDicDC = utils.ReadResult(self.config)
        LogL1 = []
        LogL0 = []
        Time = []
        for i in xrange(self.Nbin):
            CurConfig = get_config(self.configfile[i])
            #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
            try :
                ResultDic = utils.ReadResult(CurConfig)
            except :
                print "WARNING : fail reading the config file : ",CurConfig
                print "Job Number : ",i
                print "Please have a look at this job log file"
                continue

#            LogL1.append(ResultDic.get("log_like"))
            #Update the time and time error array
            Time.append((ResultDic.get("tmax")+ResultDic.get("tmin"))/2.)

            ##############################################################
            #   Compute the loglike value using the DC flux or prefactor
            ##############################################################
            # Create one obs instance
            CurConfig['Spectrum']['FitsGeneration'] = 'no'
            _,Fit = GenAnalysisObjects(CurConfig,verbose=0)#be quiet
            Fit.ftol = float(self.config['fitting']['ftol'])

            #Spectral index management!
            utils.FreezeParams(Fit, self.srcname, 'Index', -2)
            LogL1.append(-Fit.fit(0,optimizer=CurConfig['fitting']['optimizer']))

            Model_type = Fit.model.srcs[self.srcname].spectrum().genericName()
            if (Model_type == 'PowerLaw') :
                utils.FreezeParams(Fit, self.srcname, 'Prefactor', utils.fluxNorm(ResultDicDC['Prefactor']))
            if (Model_type == 'PowerLaw2') :
                utils.FreezeParams(Fit, self.srcname, 'Integral', utils.fluxNorm(ResultDicDC['Integral']))
            LogL0.append(-Fit.fit(0,optimizer=CurConfig['fitting']['optimizer']))

            del Fit #Clean memory

        Can = _GetCanvas()
        TgrDC = ROOT.TGraph(len(Time),np.array(Time),np.array(LogL0))
        TgrDC.Draw("ALP*")
        TgrDC = ROOT.TGraph(len(Time),np.array(Time),np.array(LogL0))
        TgrDC.SetMarkerColor(2)
        TgrDC.Draw("PL*")
        #Save the canvas in the LightCurve subfolder
        Can.Print(LcOutPath+'_VarIndex.eps')
        Can.Print(LcOutPath+'_VarIndex.C')
        print 
        print "\t TSvar = ",2*(sum(LogL1)-sum(LogL0))
        print "\t NDF = ",len(LogL0)-1
        print "\t Chi2 prob = ",ROOT.TMath.Prob(2*(sum(LogL1)-sum(LogL0)),len(LogL0)-1)
        print 
def _GetCanvas():
    Canv = ROOT.TCanvas()
    Canv.SetGridx()
    Canv.SetGridy()
    return Canv

def WriteToAscii(Time, TimeErr, Flux, FluxErr, TS, Npred, filename):
    """Write the results of the LC in a Ascii file"""
    flc = open(filename, 'w')
    flc.write('# Time (MET) Delta_Time Flux(ph cm-2 s-1) '
              'Delta_Flux TS Npred\n')
    for i in xrange(len(Time)):
        flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                  str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                  str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    flc.close()

