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
from enrico.config import get_config
from enrico.constants import LightcurvePath,FoldedLCPath
from enrico.submit import call
from enrico.RunGTlike import run, GenAnalysisObjects
from enrico import Loggin
from enrico.plotting import plot_errorbar_withuls
import astropy.io.fits as pyfits

pol0 = lambda x,p1: p1
pol1 = lambda x,p1,p2: p1+p2*x


class LightCurve(Loggin.Message):
    """Class to calculate light curves and variability indexes."""
    def __init__(self, config, parent_filename=""):
        super(LightCurve,self).__init__()
        Loggin.Message.__init__(self)
        self.parent_filename = os.path.abspath(parent_filename)
        self.config        = get_config(config)
        self.generalconfig = get_config(config)
        print((self.generalconfig))
        #Read the config
        self.srcname = self.config['target']['name'] #src name
        self.Tag = self.config['file']['tag']
        self.tmin = self.config['time']['tmin']
        self.tmax = self.config['time']['tmax']

        self.submit = self.config['Submit']
        # One point of the LC will be computed as a spectrum plot.
        # enrico_sed will be used
        # Do fits files will be generated
        #self.config['target']['spectrum'] = 'PowerLaw' # simplify the spectrum
        self.config['Spectrum']['FitsGeneration'] = self.config['LightCurve']['FitsGeneration']
        #Freeze the Spectral index at a value of self.config['LightCurve']['SpectralIndex'] (no effect if 0)
        self.config['Spectrum']['FrozenSpectralIndex'] = self.config['LightCurve']['SpectralIndex']
        if (self.config['LightCurve']['SpectralIndex'] != 0):
            self.config['UpperLimit']['SpectralIndex'] = self.config['LightCurve']['SpectralIndex']
        #TS limit. Compute an UL if the TS is below TSLightCurve
        self.config['UpperLimit']['TSlimit'] = self.config['LightCurve']['TSLightCurve']

        self.folder = self.config['out']
        # Do not create plots
        self.config['Spectrum']['ResultPlots'] = 'no' # no
        self.config['Spectrum']['ResultParentPlots'] = 'no' # no
        self.config['Ebin']['NumEnergyBins'] = 0
        self.config['energy']['decorrelation_energy'] = 'yes' # no
        self.config['UpperLimit']['envelope'] = 'no'
        # Submission will be directly handle by this soft
        self.config['Submit'] = 'no'
        # self.config['verbose'] ='no' #Be quiet

        # Speed-up the analysis by reusing the evt file from the main analysis
        self._RecycleEvtCoarse()

        self.configfile = [] #All the config file in the disk are stored in a list

    def _RecycleEvtCoarse(self):
        ''' Try to guess if there's an EvtCoarse file with the events extracted, reuse it '''
        import os.path
        evtcoarsefile = str("%s/%s_%s_EvtCoarse.fits"%(self.folder,self.srcname,self.Tag))
        if os.path.isfile(evtcoarsefile):
            print(("reusing %s as event file to speed-up the analysis" %evtcoarsefile))
            self.config['file']['event'] = evtcoarsefile

    def _MakeTimeBins(self):
        self.time_array = np.zeros(0)
        self.Nbin = 0
        self.gtifile = []
        if self.config['time']['file'] != '': 
            if ".fit" not in self.config['time']['file']:
                # Assume it is a text file
                print(("use "+self.config['time']['file']))
                self.gtifile.append(self.config['time']['file'])
                times = np.genfromtxt(self.gtifile[0],dtype="float",unpack=True)
                self.Nbin = int(times.size/2)
                self.time_array=np.reshape(times,times.size,'F')
                    
                if self.config['time']['type']=='MJD':
                     self.time_array = utils.MJD_to_met(self.time_array)
                elif self.config['time']['type']=='JD':
                     self.time_array = utils.JD_to_met(self.time_array)
            else:
                # Assume it is a catalog.fits file
                # get from the header the BEGIN and END time 
                with pyfits.open(self.config['time']['file']) as catfile:
                    self.tmin = catfile[1].header['TSTART']
                    self.tmax = catfile[1].header['TSTOP']
                    self.Nbin = self.config['LightCurve']['NLCbin']
                    self.time_array = np.zeros(self.Nbin*2)
                    t = np.arange(self.tmin,self.tmax+1e-5,\
                        (self.tmax - self.tmin) / self.Nbin)
                    for i in range(self.Nbin):
                        self.time_array[2*i] = t[i]
                        self.time_array[2*i+1]= t[i+1]
        
        else:
            self.Nbin = int(self.config['LightCurve']['NLCbin'])
            self.time_array = np.zeros(self.Nbin*2)
#            self.dt = (self.tmax - self.tmin) / self.Nbin
            t = np.arange(self.tmin,self.tmax+0.000001,(self.tmax - self.tmin) / self.Nbin)
            for i in range(self.Nbin):
                self.time_array[2*i] = t[i]
                self.time_array[2*i+1]= t[i+1]

        self.info("Running LC with "+str(self.Nbin)+" bins")
        for i in range(self.Nbin):
            print(("Bin ",i," Start=",self.time_array[2*i]," Stop=",self.time_array[2*i+1]))
        print()

    def _errorReading(self,message,i):
        self.warning(message+" : "+self.configfile[i])
        print("Job Number : "+str(i))
        self.warning("Please have a look at this job log file")


    def _ManageFolder(self,path):
        """   All files will be stored in a subfolder name path + NLCbin
        Create a subfolder"""
        self.LCfolder =  self.folder+"/LightCurve_"+str(self.Nbin)+"bins/"
        utils.mkdir_p(self.LCfolder)
        self.config['out'] = self.LCfolder

    def PrepareLC(self,write = 'no'):
        """Simple function to prepare the LC generation : generate and write the config files"""
        for i in range(self.Nbin):
            self.config['time']['tmin'] = self.time_array[2*i]
            self.config['time']['tmax'] = self.time_array[2*i+1]
            self.config['file']['tag'] = self.Tag + '_LC_' + str(i)
            self.config['target']['spectrum'] = 'PowerLaw' # simplify the spectrum
            filename = (self.config['out'] + "Config_" + str(i) + "_" +
                    str(self.config['time']['tmin']) + "_" +
                    str(self.config['time']['tmax']))#Name of the config file
            xmlfilename = (self.config['out'] + "" + str(i) + "_" +
                    str(self.config['time']['tmin']) + "_" +
                    str(self.config['time']['tmax'])) + ".xml" #Name of the xml file
            self.config['file']['xml'] = xmlfilename
            # Do not produce spectral plots

            if len(self.gtifile)==1:
                self.config['time']['file']=self.gtifile[0]
            elif len(self.gtifile)>1:
                print(('Time selection file for bin {0} = {1}'.format(i,self.gtifile[i])))
                self.config['time']['file']=self.gtifile[i]

            if write == 'yes':
                self.config.write(open(filename, 'wb'))

            self.configfile.append(filename)

    def _MakeLC(self,Path=LightcurvePath) :
        #import gc
        import os
        #gc.enable()
        '''Main function of the Lightcurve script. Read the config file and run the gtlike analysis'''
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')

        self.PrepareLC(self.config['LightCurve']['MakeConfFile'])#Get the config file

        for i in range(self.Nbin):
            #gc.collect()
            cmd = str("enrico_sed %s && enrico_plot_lc %s" %(self.configfile[i], self.parent_filename))
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
        self.info("Reading files produced by enrico")
        LcOutPath = self.LCfolder + self.config['target']['name']

        #Result are stored into list. This allow to get rid of the bin which failled
        Time = []
        TimeErr = []
        Flux = []
        FluxErr = []
        # FluxErrChi2 = []
        Index = []
        IndexErr = []
        Cutoff = []
        CutoffErr = []
        FluxForNpred = []
        # FluxErrForNpred = []
        Npred = []
        Npred_detected_indices = []
        TS = []
        uplim = []

        # Find name used for index parameter
        if ((self.config['target']['spectrum'] == 'PowerLaw' or
             self.config['target']['spectrum'] == 'PowerLaw2') and
             self.config['target']['redshift'] == 0):
                IndexName = 'Index'
                CutoffName = None
        elif (self.config['target']['spectrum'] == 'PLExpCutoff' or
                self.config['target']['spectrum'] == 'PLSuperExpCutoff'):
            IndexName = 'Index1'
            CutoffName = 'Cutoff'
            CutoffErrName = 'dCutoff'
        else:
            IndexName = 'alpha'
            CutoffName = None
        IndexErrName = 'd' + IndexName

        Nfail = 0
        for i in range(self.Nbin):
            CurConfig = get_config(self.configfile[i])
            #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
            try :
                ResultDic = utils.ReadResult(CurConfig)
                if ResultDic == {}:
                    raise(ValueError)
            except :
                self._errorReading("Fail reading config file",i)
                Nfail+=1
                continue

            #Update the time and time error array
            Time.append((ResultDic.get("tmax")+ResultDic.get("tmin"))/2.)
            TimeErr.append((ResultDic.get("tmax")-ResultDic.get("tmin"))/2.)
            #Check is an ul have been computed. The error is set to zero for the TGraph.
            if 'Ulvalue' in ResultDic :
                uplim.append(1)
                Flux.append(ResultDic.get("Ulvalue"))
                # FluxErr.append(0)
                # FluxErrChi2.append(ResultDic.get("dFlux"))
                # Index.append(ResultDic.get(IndexName))
                # IndexErr.append(0)
            else :
                uplim.append(0)
                Flux.append(ResultDic.get("Flux"))
            FluxErr.append(ResultDic.get("dFlux"))
            # FluxErrChi2.append(ResultDic.get("dFlux"))
            Index.append(ResultDic.get(IndexName))
            IndexErr.append(ResultDic.get(IndexErrName))
                # if CutoffName is not None:
                    # Cutoff.append(ResultDic.get(CutoffName))
                    # CutoffErr.append(ResultDic.get(CutoffErrName))
            # FluxErrForNpred.append(ResultDic.get("dFlux"))
            FluxForNpred.append(ResultDic.get("Flux"))
            #Get the Npred and TS values
            Npred.append(ResultDic.get("Npred"))
            TS.append(ResultDic.get("TS"))
            if (CurConfig['LightCurve']['TSLightCurve']<float(ResultDic.get("TS"))):
                Npred_detected_indices.append(i-Nfail)

        # #change the list into np array
        # TS = np.array(TS)
        Npred = np.asarray(Npred)
        Npred_detected = np.asarray(Npred[Npred_detected_indices])
        Time = np.asarray(Time)
        TimeErr = np.asarray(TimeErr)
        Flux = np.asarray(Flux)
        FluxErr = np.asarray(FluxErr)
        # Index = np.array(Index)
        # IndexErr = np.array(IndexErr)
        # Cutoff = np.array(Cutoff)
        # CutoffErr = np.array(CutoffErr)
        FluxForNpred = np.asarray(FluxForNpred)
        # FluxErrForNpred = np.array(FluxErrForNpred)
        uplim = np.asarray(uplim,dtype=bool)
        #Plots the diagnostic plots is asked
        # Plots are : Npred vs flux
        #             TS vs Time
        if self.config['LightCurve']['DiagnosticPlots'] == 'yes' and len(Npred)>0:
            #plot Npred vs flux
            plt.figure()
            NdN = np.asarray(Npred) /np.sqrt(Npred)
            FdF = np.asarray(FluxForNpred) / (np.asarray(FluxErr) + 1e-20)
            plt.errorbar(NdN, FdF,fmt='+',color='black')

            if len(Npred_detected)>2:
                NdN = np.asarray(Npred_detected) /np.sqrt(Npred_detected)
                FdF = np.asarray(FluxForNpred[Npred_detected_indices]) / (np.asarray(FluxErr[Npred_detected_indices]) + 1e-20)
                plt.errorbar(NdN, FdF,fmt='+',color='red')

                popt,_ = scipy.optimize.curve_fit(pol1, NdN, FdF, p0=[0,1])#, sigma=dydata)


                for i in range(len(FluxForNpred)):
                    if FluxForNpred[i]/FluxErr[i]>2*pol1(sqrt(Npred[i]),popt[0],popt[1]):
                        self._errorReading("problem in errors calculation for",i)
                        print(("Flux +/- error = ",FluxForNpred[i]," +/- ",FluxErr[i]))
                        print(("V(Npred) = ",sqrt(Npred[i])))
                        print()

                plt.plot(np.array([0,max(NdN)]),pol1(np.array([0,max(NdN)]),popt[0],popt[1]),'--',color='black')
                plt.xlabel(r"${\rm Npred/\sqrt{Npred}}$")
                plt.ylabel(r"${\rm Flux/\Delta Flux}$")
                plt.savefig(LcOutPath+"_Npred.png", dpi=150, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)
            else :
                print("No Npred Plot produced")

            #plot TS vs Time
            plt.figure()
            plt.xlabel(r"Time (s)")
            plt.ylabel(r"Test Statistic")
            plt.errorbar(x=Time,y=TS,xerr=TimeErr,fmt='+',color='black',ls='None')
            plt.ylim(ymin=min(TS)*0.8,ymax=max(TS)*1.2)
            plt.xlim(xmin=max(plt.xlim()[0],1.02*min(Time)-0.02*max(Time)),xmax=min(plt.xlim()[1],1.02*max(Time)-0.02*min(Time)))

            # Move the offset to the axis label
            ax = plt.gca()
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            offset_factor = int(np.mean(np.log10(np.abs(ax.get_ylim()))))
            if (offset_factor != 0):
                ax.set_yticklabels([float(round(k,5)) for k in ax.get_yticks()*10**(-offset_factor)])
                ax.yaxis.set_label_text(ax.yaxis.get_label_text() + r" [${\times 10^{%d}}$]" %offset_factor)

            # Secondary axis with MJD
            mjdaxis = ax.twiny()
            mjdaxis.set_xlim([utils.met_to_MJD(k) for k in ax.get_xlim()])
            mjdaxis.set_xlabel(r"Time (MJD)")
            mjdaxis.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
            plt.setp( mjdaxis.xaxis.get_majorticklabels(), rotation=15 )
            plt.tight_layout()

            plt.savefig(LcOutPath+"_TS.png", dpi=150, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)


#    Plot the LC itself. This function return a TH2F for a nice plot
#    a TGraph and a list of TArrow for the ULs
        # if folded:
        #     phase = np.linspace(0,1,self.Nbin+1)
        #     Time = (phase[1:]+phase[:-1])/2.
        #     TimeErr = (phase[1:]-phase[:-1])/2.
        #     gTHLC,TgrLC,ArrowLC = plotting.PlotFoldedLC(Time,TimeErr,Flux,FluxErr)
        #     gTHIndex,TgrIndex,ArrowIndex = plotting.PlotFoldedLC(Time,TimeErr,Index,IndexErr)
        #     if CutoffName is not None:
        #         gTHCutoff,TgrCutoff,ArrowCutoff = plotting.PlotFoldedLC(Time,TimeErr,Cutoff,CutoffErr)
        # else :
        #     gTHLC,TgrLC,ArrowLC = plotting.PlotLC(Time,TimeErr,Flux,FluxErr)
        #     gTHIndex,TgrIndex,ArrowIndex = plotting.PlotLC(Time,TimeErr,Index,IndexErr)
        #     if CutoffName is not None:
        #         gTHCutoff,TgrCutoff,ArrowCutoff = plotting.PlotFoldedLC(Time,TimeErr,Cutoff,CutoffErr)

        # xmin = min(Time) - max(TimeErr) * 10
        # xmax = max(Time) + max(TimeErr) * 10
        # ymin = min(Flux) - max(FluxErr) * 1.3
        # ymax = max(Flux) + max(FluxErr) * 1.3
        plt.figure()
        plt.xlabel(r"Time (s)")
        plt.ylabel(r"${\rm Flux\ (photon\ cm^{-2}\ s^{-1})}$")
        # plt.ylim(ymin=ymin,ymax=ymax)
        # plt.xlim(xmin=xmin,xmax=xmax)
        #plt.errorbar(Time,Flux,xerr=TimeErr,yerr=FluxErr,i
        #             fmt='o',color='black',ls='None',uplims=uplim)
        plot_errorbar_withuls(Time,TimeErr,TimeErr,Flux,FluxErr,FluxErr,
                              uplim,bblocks=True)
        
        try:
            plt.ylim(ymin=max(plt.ylim()[0],np.percentile(Flux[~uplim],1)*0.1),
                     ymax=min(plt.ylim()[1],np.percentile(Flux[~uplim],99)*2.0))
            plt.xlim(xmin=max(plt.xlim()[0],1.02*min(Time)-0.02*max(Time)),
                     xmax=min(plt.xlim()[1],1.02*max(Time)-0.02*min(Time)))
        except IndexError:
            pass

        # Move the offset to the axis label
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        offset_factor = int(np.mean(np.log10(np.abs(ax.get_ylim()))))
        if (offset_factor != 0):
            ax.set_yticklabels([float(round(k,5)) \
              for k in ax.get_yticks()*10**(-offset_factor)])
            ax.yaxis.set_label_text(ax.yaxis.get_label_text() +\
              r" [${\times 10^{%d}}$]" %offset_factor)

        # Secondary axis with MJD
        mjdaxis = ax.twiny()
        mjdaxis.set_xlim([utils.met_to_MJD(k) for k in ax.get_xlim()])
        mjdaxis.set_xlabel(r"Time (MJD)")
        mjdaxis.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
        plt.setp( mjdaxis.xaxis.get_majorticklabels(), rotation=15 )
        plt.tight_layout()

        plt.savefig(LcOutPath+"_LC.png", dpi=150, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)

        if self.config["LightCurve"]["SpectralIndex"] == 0 :
            plt.figure()
            plt.xlabel(r"Time (s)")
            plt.ylabel(r"${\rm Index}$")
            Index = np.asarray(Index)
            IndexErr = np.asarray(IndexErr)
            uplimIndex = uplim #+ Index<0.55
            plot_errorbar_withuls(Time[~uplimIndex],
                                  TimeErr[~uplimIndex],
                                  TimeErr[~uplimIndex],
                                  Index[~uplimIndex],
                                  IndexErr[~uplimIndex],
                                  IndexErr[~uplimIndex],
                                  uplimIndex[~uplimIndex],
                                  bblocks=True)

            plt.ylim(ymin=max(plt.ylim()[0],np.percentile(Index[~uplimIndex],1)*0.1),
                     ymax=min(plt.ylim()[1],np.percentile(Index[~uplimIndex],99)*2.0))
            plt.xlim(xmin=max(plt.xlim()[0],1.02*min(Time)-0.02*max(Time)),
                     xmax=min(plt.xlim()[1],1.02*max(Time)-0.02*min(Time)))

            # Move the offset to the axis label
            ax = plt.gca()
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            offset_factor = int(np.mean(np.log10(np.abs(ax.get_ylim()))))
            if (offset_factor != 0):
                ax.set_yticklabels([float(round(k,5)) \
                  for k in ax.get_yticks()*10**(-offset_factor)])
                ax.yaxis.set_label_text(ax.yaxis.get_label_text() +\
                   r" [${\times 10^{%d}}$]" %offset_factor)

            # Secondary axis with MJD
            mjdaxis = ax.twiny()
            mjdaxis.set_xlim([utils.met_to_MJD(k) for k in ax.get_xlim()])
            mjdaxis.set_xlabel(r"Time (MJD)")
            mjdaxis.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
            plt.setp( mjdaxis.xaxis.get_majorticklabels(), rotation=15 )
            plt.tight_layout()
            plt.savefig(LcOutPath+"_Index.png", dpi=150,
                facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)


        # compute Fvar and probability of being cst

        self.info("Flux vs Time: infos")
        self.FitWithCst(Time,Flux,FluxErr)
        self.Fvar(Flux,FluxErr)

        # ### plot and save the Index LC
        # CanvIndex = ROOT.TCanvas()
        # gTHIndex.Draw()
        # TgrIndex.Draw('zP')

        # #plot the ul as arrow
        # for i in xrange(len(ArrowIndex)):
        #     ArrowIndex[i].Draw()

        # #Save the canvas in the LightCurve subfolder
        # if self.config["LightCurve"]["SpectralIndex"] == 0 :
        #     self.info("Index vs Time")
        #     self.FitWithCst(Time,Index,IndexErr)
        #     CanvIndex.Print(LcOutPath+'_Index.png')
        #     CanvIndex.Print(LcOutPath+'_Index.eps')
        #     CanvIndex.Print(LcOutPath+'_Index.C')


        #Dump into ascii
        lcfilename = LcOutPath+"_results.dat"
        self.info("Write to Ascii file : "+lcfilename)
        WriteToAscii(Time,TimeErr,Flux,FluxErr,Index,IndexErr,
                     Cutoff,CutoffErr,TS,Npred,lcfilename)

        if self.config["LightCurve"]['ComputeVarIndex'] == 'yes':
             self.VariabilityIndex()

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
            print(("\tFvar = ",fvar," +/- ",err_fvar))
        except :
            print(("\tFvar is negative, Fvar**2 = %2.2e +/- %2.2e"%(intvar/(moy*moy), ((1./sqrt(2*len(Flux))*expvar/moy**2)**2/(intvar/(moy*moy)) + (sqrt(expvar/len(Flux))*1./moy)**2))))
        print()

    def FitWithCst(self,x,y,dy):
        """Fit the LC with a constant function an
           print the chi2 and proba"""
        res,_ = scipy.optimize.curve_fit(pol0,x,y,p0=[np.mean(y)],sigma=dy)
        #print res
        #print y
        #print dy
        cost = np.sum(((pol0(x,res[0])-y)/dy)**2)
        self.info("Fit with a constant function")
        print(('\tChi2 = ',cost," NDF = ",len(y)-1))
        print(('\tprobability of being cst = ',1 - chi2.cdf(cost,len(y)-1)))
        print()

    def VariabilityIndex(self):
        """Compute the variability index as in the 2FGL catalogue. (see Nolan et al, 2012)"""
        LcOutPath = self.LCfolder + self.config['target']['name']

        utils._log('Computing Variability index ')

        self.config['Spectrum']['FitsGeneration'] = 'no'

        try :
            ResultDicDC = utils.ReadResult(self.generalconfig)
        except :
            self.warning("No results file found; please run enrico_sed first.")
            return

        LogL1 = []
        LogL0 = []
        Time = []
        for i in range(self.Nbin):
            CurConfig = get_config(self.configfile[i])
            #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
            try :
                ResultDic = utils.ReadResult(CurConfig)
            except :
                self._errorReading("Fail reading the config file ",i)
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
            parameters = dict()
            parameters['Index']  = -2.
            parameters['alpha']  = +2.
            parameters['Index1'] = -2.
            parameters['beta']   = 0
            parameters['Index2'] = 2.
            parameters['Cutoff'] = 100000. # set the cutoff to be high

            for key in list(parameters.keys()):
                try:
                    utils.FreezeParams(Fit, self.srcname, key, parameters[key])
                except:
                    continue

            LogL1.append(-Fit.fit(0,optimizer=CurConfig['fitting']['optimizer']))

            for key in ["norm","Prefactor","Integral"]:
                try:
                    utils.FreezeParams(Fit,self.srcname,key, utils.fluxNorm(ResultsDicDC[key]))
                except:
                    continue

            LogL0.append(-Fit.fit(0,optimizer=CurConfig['fitting']['optimizer']))

            del Fit #Clean memory


        plt.figure()
        plt.xlabel("Time")
        plt.ylabel("Log(Like) Variability")
        plt.errorbar(Time,LogL0,fmt='o',color='black',ls='None')
        plt.xlim(xmin=max(plt.xlim()[0],1.02*min(Time)-0.02*max(Time)),
                 xmax=min(plt.xlim()[1],1.02*max(Time)-0.02*min(Time)))

        # Move the offset to the axis label
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        offset_factor = int(np.mean(np.log10(np.abs(ax.get_ylim()))))
        if (offset_factor != 0):
            ax.set_yticklabels([float(round(k,5)) \
              for k in ax.get_yticks()*10**(-offset_factor)])
            ax.yaxis.set_label_text(ax.yaxis.get_label_text() +\
               r" [${\times 10^{%d}}$]" %offset_factor)

        # Secondary axis with MJD
        mjdaxis = ax.twiny()
        mjdaxis.set_xlim([utils.met_to_MJD(k) for k in ax.get_xlim()])
        mjdaxis.set_xlabel(r"Time (MJD)")
        mjdaxis.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
        plt.setp( mjdaxis.xaxis.get_majorticklabels(), rotation=15 )
        plt.tight_layout()

        plt.savefig(LcOutPath+"_VarIndex.png", dpi=150, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)

        self.info("Variability index calculation")
        print(("\t TSvar = ",2*(sum(LogL1)-sum(LogL0))))
        print(("\t NDF = ",len(LogL0)-1))
        print(("\t Chi2 prob = ",1 - chi2.cdf(2*(sum(LogL1)-sum(LogL0)),len(LogL0)-1)))
        print()

def WriteToAscii(Time, TimeErr, Flux, FluxErr, Index, IndexErr, Cutoff, CutoffErr, TS, Npred, filename):
    """Write the results of the LC in a Ascii file"""
    flc = open(filename, 'w')
    if len(Cutoff) == 0:
        flc.write('# Time (MET) Delta_Time Flux(ph cm-2 s-1) '
                  'Delta_Flux Index Delta_Index TS Npred\n')
        for i in range(len(Time)):
            flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                      str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                      str(Index[i]) + "\t" + str(IndexErr[i]) + "\t" +
                      str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    else:
        flc.write('# Time (MET) Delta_Time Flux(ph cm-2 s-1) '
                  'Delta_Flux Index Delta_Index Cutoff Delta_Cutoff TS Npred\n')
        for i in range(len(Time)):
            flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                      str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                      str(Index[i]) + "\t" + str(IndexErr[i]) + "\t" +
                      str(Cutoff[i]) + "\t" + str(CutoffErr[i]) + "\t" +
                      str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    flc.close()
