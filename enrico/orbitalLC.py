import os
from os.path import join
import numpy as np
import ROOT
import utils
import root_style
import plotting
import environ
from math import sqrt
from enrico.config import get_config
from submit import call
from enrico.RunGTlike import run,GenAnalysisObjects
from enrico.lightcurve import LightCurve, WriteToAscii


class OrbitalLC(LightCurve):
    """Class to calculate orbital light curves and variability indexes."""
    def __init__(self, config):
        ROOT.gROOT.SetBatch(ROOT.kTRUE) #Batch mode
        self.config = config

        #Read the config
        self.srcname = self.config['target']['name'] #src name
        self.Tag     = self.config['file']['tag']
        self.Nbin    = self.config['LightCurve']['NLCbin']
        # One point of the LC will be computed as a spectrum plot.
        # enrico_sed will be used
        # Do fits files will be generated
        self.config['Spectrum']['FitsGeneration'] = self.config['LightCurve']['FitsGeneration']
        #Froze the Spectral index at a value of self.config['LightCurve']['SpectralIndex']
        self.config['Spectrum']['FrozenSpectralIndex'] = self.config['LightCurve']['SpectralIndex']
        #TS limit. Compute an UL if the TS is below TSLightCurve
        self.config['UpperLimit']['TSlimit'] = self.config['LightCurve']['TSLightCurve']

        self.folder = self.config['out']
        #All files will be stored in a subfolder name OrbitalLC + NLCbin
        #Create a subfolder name OrbitalLC
        self.LCfolder =  self.folder+"/OrbitalLC_"+str(self.config['LightCurve']['NLCbin'])+"bins/"
        os.system("mkdir -p "+self.LCfolder)
        self.config['out'] = self.LCfolder

        #No plot, no bin in energy, Normal UL
        self.config['Spectrum']['ResultPlots'] = 'no'
        self.config['Ebin']['NumEnergyBins'] = 0
        self.config['UpperLimit']['envelope'] = 'no'
        #No submition. Submission will be directly handle by this soft
        self.config['Spectrum']['Submit'] = 'no'
        self.config['verbose'] ='no' #Be quiet

        self.configfile = []#All the config file in the disk are stored in a list

    def PrepareLC(self,write = 'no'):
        """Simple function to prepare the LC generation : generate and write the config files"""

        for i in xrange(self.Nbin):
            self.config['OrbitalLC']['phasemin'] = i / float(self.Nbin)
            self.config['OrbitalLC']['phasemax'] = (i + 1) / float(self.Nbin)
            self.config['file']['tag'] = self.Tag + '_LC_' + str(i)
            filename = (self.config['out'] + "LC_{0}_{1:03.0f}_{2:03.0f}.conf".format(
                    i,
                    self.config['OrbitalLC']['phasemin']*1000,
                    self.config['OrbitalLC']['phasemax']*1000,))#Name of the config file
            if write == 'yes':
                self.config.write(open(filename, 'w'))

            self.configfile.append(filename)

    def PlotLC(self):
        '''Plot a lightcurve which have been generated previously'''
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

        self.PrepareLC()#Get the config file

        for i in xrange(self.Nbin):
            CurConfig = get_config(self.configfile[i])
            #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
            try :
                ResultDic = utils.ReadResult(CurConfig)
            except :
                self._errorReading("fail reading config file",i)
                continue

            #Update the time and time error array
            Time.append((ResultDic.get("phasemax")+ResultDic.get("phasemin"))/2.)
            TimeErr.append((ResultDic.get("phasemax")-ResultDic.get("phasemin"))/2.)
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
        # Flux values must be multiplied by Nbin because phase selection through
        # gtpphase is not taken into account by gtmktime
        Flux = np.array(Flux)*self.Nbin
        FluxErr = np.array(FluxErr)*self.Nbin
        FluxForNpred = np.array(FluxForNpred)*self.Nbin
        FluxErrForNpred = np.array(FluxErrForNpred)*self.Nbin

        fittedFunc = self.CheckNpred(Npred,FluxForNpred,FluxErrForNpred,Npred_detected_indices)#check the errors calculation

        #Plots the diagnostic plots is asked
        # Plots are : Npred vs flux
        #             TS vs Time
        if self.config['LightCurve']['DiagnosticPlots'] == 'yes':
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

            gTHTS,TgrTS = plotting.PlotOrbTS(Time,TimeErr,TS)
            CanvTS = _GetCanvas()
            gTHTS.Draw()
            TgrTS.Draw('zP')
            CanvTS.Print(LcOutPath+'_TS.eps')
            CanvTS.Print(LcOutPath+'_TS.C')

#    Plot the LC itself. This function return a TH2F for a nice plot
#    a TGraph and a list of TArrow for the ULs
        gTHLC,TgrLC,ArrowLC = plotting.PlotOrbLC(Time,TimeErr,Flux,FluxErr)
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


