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
from enrico import Loggin

class LightCurve(Loggin.Message):
    """Class to calculate light curves and variability indexes."""
    def __init__(self, config):
        super(LightCurve,self).__init__()
        Loggin.Message.__init__(self)
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
            elif self.config['time']['type']=='JD':
                 self.time_array = utils.JD_to_met(self.time_array)
        else:
            self.Nbin = self.config['LightCurve']['NLCbin']
            self.time_array = np.zeros(self.Nbin*2)
#            self.dt = (self.tmax - self.tmin) / self.Nbin
            t = np.arange(self.tmin,self.tmax+0.000001,(self.tmax - self.tmin) / self.Nbin)
            for i in xrange(self.Nbin):
                self.time_array[2*i] = t[i]
                self.time_array[2*i+1]= t[i+1]

        self.info("Running LC with "+str(self.Nbin)+" bins")
        for i in xrange(self.Nbin):
            print "Bin ",i," Start=",self.time_array[2*i]," Stop=",self.time_array[2*i+1]
        print 

    def _errorReading(self,message,i):
        self.warning(message+" : "+self.configfile[i])
        print "Job Number : ",i
        self.warning("Please have a look at this job log file")


    def _ManageFolder(self,path):
        """   All files will be stored in a subfolder name path + tag + NLCbin
        Create a subfolder"""
        self.LCfolder =  self.folder+"/"+path+"_"+self.config['file']['tag']+"_"+str(self.Nbin)+"bins/"
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


    def _MakeLC(self) :
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

        self.info("Reading files produced by enrico")
        LcOutPath = self.LCfolder + self.config['target']['name']

        #Result are stored into list. This allow to get rid of the bin which failled
        Time = []
        TimeErr = []
        Flux = []
        FluxErr = []
        Index = []
        IndexErr = []
        Cutoff = []
        CutoffErr = []
        FluxForNpred = []
        FluxErrForNpred = []
        Npred = []
        Npred_detected_indices = []
        TS = []

        # Find name used for index parameter
        if (self.config['target']['spectrum'] == 'PowerLaw' or
                self.config['target']['spectrum'] == 'PowerLaw2'):
            IndexName = 'Index'
            CutoffName = None
        elif (self.config['target']['spectrum'] == 'PLExpCutoff' or
                self.config['target']['spectrum'] == 'PLSuperExpCutoff'):
            IndexName = 'Index1'
            CutoffName = 'Cutoff'
            CutoffErrName = 'dCutoff'
        IndexErrName = 'd' + IndexName

        Nfail = 0
        for i in xrange(self.Nbin):
            CurConfig = get_config(self.configfile[i])
            #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
            try :
                ResultDic = utils.ReadResult(CurConfig)
            except :
                self._errorReading("Fail reading config file",i)
                Nfail+=1
                continue

            #Update the time and time error array
            Time.append((ResultDic.get("tmax")+ResultDic.get("tmin"))/2.)
            TimeErr.append((ResultDic.get("tmax")-ResultDic.get("tmin"))/2.)
            #Check is an ul have been computed. The error is set to zero for the TGraph.
            if ResultDic.has_key('Ulvalue') :
                Flux.append(ResultDic.get("Ulvalue"))
                FluxErr.append(0)
                Index.append(ResultDic.get(IndexName))
                IndexErr.append(0)
            else :
                Flux.append(ResultDic.get("Flux"))
                FluxErr.append(ResultDic.get("dFlux"))
                Index.append(ResultDic.get(IndexName))
                IndexErr.append(ResultDic.get(IndexErrName))
                if CutoffName is not None:
                    Cutoff.append(ResultDic.get(CutoffName))
                    CutoffErr.append(ResultDic.get(CutoffErrName))
            FluxErrForNpred.append(ResultDic.get("dFlux"))
            FluxForNpred.append(ResultDic.get("Flux"))
            #Get the Npred and TS values
            Npred.append(ResultDic.get("Npred"))
            TS.append(ResultDic.get("TS"))
            if (CurConfig['LightCurve']['TSLightCurve']<float(ResultDic.get("TS"))):
                Npred_detected_indices.append(i-Nfail)

        #change the list into np array
        TS = np.array(TS)
        Npred = np.array(Npred)
        Npred_detected = Npred[Npred_detected_indices]
        Time = np.array(Time)
        TimeErr = np.array(TimeErr)
        Flux = np.array(Flux)
        FluxErr = np.array(FluxErr)
        Index = np.array(Index)
        IndexErr = np.array(IndexErr)
        Cutoff = np.array(Cutoff)
        CutoffErr = np.array(CutoffErr)
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

            CanvNpred.Print(LcOutPath+"_Npred.png")
            CanvNpred.Print(LcOutPath+"_Npred.eps")
            CanvNpred.Print(LcOutPath+"_Npred.C")

            gTHTS,TgrTS = plotting.PlotTS(Time,TimeErr,TS)
            CanvTS = _GetCanvas()
            gTHTS.Draw()
            TgrTS.Draw('zP')
            CanvTS.Print(LcOutPath+'_TS.png')
            CanvTS.Print(LcOutPath+'_TS.eps')
            CanvTS.Print(LcOutPath+'_TS.C')


#    Plot the LC itself. This function return a TH2F for a nice plot
#    a TGraph and a list of TArrow for the ULs
        if folded:
            phase = np.linspace(0,1,self.Nbin+1)
            Time = (phase[1:]+phase[:-1])/2.
            TimeErr = (phase[1:]-phase[:-1])/2.
            gTHLC,TgrLC,ArrowLC = plotting.PlotFoldedLC(Time,TimeErr,Flux,FluxErr)
            gTHIndex,TgrIndex,ArrowIndex = plotting.PlotFoldedLC(Time,TimeErr,Index,IndexErr)
            if CutoffName is not None:
                gTHCutoff,TgrCutoff,ArrowCutoff = plotting.PlotFoldedLC(Time,TimeErr,Cutoff,CutoffErr)
        else :
            gTHLC,TgrLC,ArrowLC = plotting.PlotLC(Time,TimeErr,Flux,FluxErr)
            gTHIndex,TgrIndex,ArrowIndex = plotting.PlotLC(Time,TimeErr,Index,IndexErr)
            if CutoffName is not None:
                gTHCutoff,TgrCutoff,ArrowCutoff = plotting.PlotFoldedLC(Time,TimeErr,Cutoff,CutoffErr)

        ### plot and save the flux LC
        CanvLC = ROOT.TCanvas()
        gTHLC.Draw()
        TgrLC.Draw('zP')

        #plot the ul as arrow
        for i in xrange(len(ArrowLC)):
            ArrowLC[i].Draw()

        # compute Fvar and probability of being cst

        self.info("Flux vs Time: infos")
        self.FitWithCst(TgrLC)
        self.Fvar(Flux,FluxErr)

        #Save the canvas in the LightCurve subfolder
        CanvLC.Print(LcOutPath+'_LC.png')
        CanvLC.Print(LcOutPath+'_LC.eps')
        CanvLC.Print(LcOutPath+'_LC.C')

        ### plot and save the Index LC
        CanvIndex = ROOT.TCanvas()
        gTHIndex.Draw()
        TgrIndex.Draw('zP')

        #plot the ul as arrow
        for i in xrange(len(ArrowIndex)):
            ArrowIndex[i].Draw()

        #Save the canvas in the LightCurve subfolder
        if self.config["LightCurve"]["SpectralIndex"] == 0 :
            self.info("Index vs Time")
            self.FitWithCst(TgrIndex)
            CanvIndex.Print(LcOutPath+'_Index.png')
            CanvIndex.Print(LcOutPath+'_Index.eps')
            CanvIndex.Print(LcOutPath+'_Index.C')

        if len(Cutoff) > 0:
            ### plot and save the Cutoff LC
            CanvCutoff = ROOT.TCanvas()
            gTHCutoff.Draw()
            TgrCutoff.Draw('zP')

            #plot the ul as arrow
            for i in xrange(len(ArrowCutoff)):
                ArrowCutoff[i].Draw()

            print "Cutoff vs Time: infos"
            self.FitWithCst(TgrCutoff)
            CanvCutoff.Print(LcOutPath+'_Cutoff.png')
            CanvCutoff.Print(LcOutPath+'_Cutoff.eps')
            CanvCutoff.Print(LcOutPath+'_Cutoff.C')

        #Dump into ascii
        lcfilename = LcOutPath+"_results.dat"
        self.info("Write to Ascii file : "+lcfilename)
        WriteToAscii(Time,TimeErr,Flux,FluxErr,Index,IndexErr,Cutoff,CutoffErr,TS,Npred,lcfilename)

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
        try :
            fvar=sqrt(intvar)/moy
            err_fvar = sqrt( ( 1./sqrt(2*len(Flux))*expvar/moy**2/fvar)**2 + (sqrt(expvar/len(Flux))*1./moy)**2)
            self.info("Calculation of Fvar (Vaughan et al. 2003)")
            print "\tFvar = ",fvar," +/- ",err_fvar
        except :
            print  "\tFvar is negative, Fvar**2 = %2.2e +/- %2.2e"%(intvar/(moy*moy), ((1./sqrt(2*len(Flux))*expvar/moy**2)**2/(intvar/(moy*moy)) + (sqrt(expvar/len(Flux))*1./moy)**2))
        print 

    def FitWithCst(self,tgrFlux):
        """Fit the LC with a constant function an
           print the chi2 and proba"""
        func = ROOT.TF1('func','pol0')
        func.SetLineColor(15)
        func.SetLineStyle(3)
        tgrFlux.Fit('func','Q')
        self.info("Fit with a constant function")
        print '\tChi2 = ',func.GetChisquare()," NDF = ",func.GetNDF()
        print '\tprobability of being cst = ',func.GetProb()
        print
        del func

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
                self._errorReading("fail reading the config file ",i)
#                print "WARNING : fail reading the config file : ",CurConfig
#                print "Job Number : ",i
#                print "Please have a look at this job log file"
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
            self.info("Spectral index frozen to a value of 2")
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
        self.info("Variability index calculation") 
        print "\t TSvar = ",2*(sum(LogL1)-sum(LogL0))
        print "\t NDF = ",len(LogL0)-1
        print "\t Chi2 prob = ",ROOT.TMath.Prob(2*(sum(LogL1)-sum(LogL0)),len(LogL0)-1)
        print 
def _GetCanvas():
    Canv = ROOT.TCanvas()
    Canv.SetGridx()
    Canv.SetGridy()
    return Canv

def WriteToAscii(Time, TimeErr, Flux, FluxErr, Index, IndexErr, Cutoff, CutoffErr, TS, Npred, filename):
    """Write the results of the LC in a Ascii file"""
    flc = open(filename, 'w')
    if len(Cutoff) == 0:
        flc.write('# Time (MET) Delta_Time Flux(ph cm-2 s-1) '
                  'Delta_Flux Index Delta_Index TS Npred\n')
        for i in xrange(len(Time)):
            flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                      str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                      str(Index[i]) + "\t" + str(IndexErr[i]) + "\t" +
                      str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    else:
        flc.write('# Time (MET) Delta_Time Flux(ph cm-2 s-1) '
                  'Delta_Flux Index Delta_Index Cutoff Delta_Cutoff TS Npred\n')
        for i in xrange(len(Time)):
            flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                      str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                      str(Index[i]) + "\t" + str(IndexErr[i]) + "\t" +
                      str(Cutoff[i]) + "\t" + str(CutoffErr[i]) + "\t" +
                      str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    flc.close()

