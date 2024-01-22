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
from enrico.plotting import plot_bayesianblocks

pol0 = lambda x,p1: p1
pol1 = lambda x,p1,p2: p1+p2*x


class BayesianBlocks(lightcurve.LightCurve):
    """Class to calculate light curves and variability indexes."""
    def __init__(self, config, parent_filename=""):
        super(BayesianBlocks, self).__init__(config, parent_filename)

        self.LCfolder =  self.folder+"/BayesianBlocks/"
        utils.mkdir_p(self.LCfolder)

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
        expfile = str("%s/AppertureLightCurve/%s_%s_applc.fits"%(self.folder,self.srcname,self.Tag))
        apfile = str("%s/AppertureLightCurve/TimeExposureCount.txt"%(self.folder))
        if not os.path.isfile(evtfile) or not os.path.isfile(expfile) or not os.path.isfile(apfile):
            raise Exception('The apperture photometry events list doesn\'t exist\nPlease run enrico_applc first')

    def _MakeTimeBins(self):
        from astropy.stats import bayesian_blocks
        from astropy.table import Table

        evtfile = str("%s/AppertureLightCurve/%s_%s_MkTime.fits"%(self.folder,self.srcname,self.Tag))
        evtlist = Table.read(evtfile, hdu='EVENTS')['TIME'].data
        expfile = str("%s/AppertureLightCurve/%s_%s_applc.fits"%(self.folder,self.srcname,self.Tag))
        expbins = Table.read(expfile, hdu='RATE')

        meanRate = float(len(evtlist))/float(self.tmax-self.tmin)
        print(("Mean photon rate %s s^-1" %meanRate))
        print(("Mean photon rate %s day^-1" %(meanRate*3600.*24)))

        #Sort table in function of time just to be sure
        evtlist.sort()
        evtlistExpCorrected = np.empty_like(evtlist)
        expbins[expbins.argsort('TIME')]

        # Calculate relative exposure time and time correction associated for each exposure bins
        j=0
        surfaceFermi = 10000 # in cm^2
        timeCorrection = np.zeros((len(expbins)+1,2))
        exposure = np.zeros(len(expbins))
        timeCorrection[j, 0] = expbins['TIME'][j]-0.5*expbins['TIMEDEL'][j]
        timeCorrection[j, 1] = 0.
        exposure[j] = expbins['EXPOSURE'][j]/(expbins['TIMEDEL'][j]*surfaceFermi)
        for j in range(1, len(expbins)):
            exposure[j] = expbins['EXPOSURE'][j]/(expbins['TIMEDEL'][j]*surfaceFermi)
            timeCorrection[j, 0] = expbins['TIME'][j]-0.5*expbins['TIMEDEL'][j]
            timeCorrection[j, 1] = timeCorrection[j-1, 1]+exposure[j-1]*expbins['TIMEDEL'][j-1]
        timeCorrection[j+1, 0] = expbins['TIME'][j]+0.5*expbins['TIMEDEL'][j]
        timeCorrection[j+1, 1] = timeCorrection[j, 1]+exposure[j]*expbins['TIMEDEL'][j]

        #Apply exposure time correction
        evtlistcorrected = np.interp(evtlist, timeCorrection[:, 0], timeCorrection[:, 1])
        meanRateCorrected = float(len(evtlistcorrected))/float(timeCorrection[-1, 1]-timeCorrection[0, 1])
        print(("Mean photon rate exposure corrected %s s^-1" %meanRateCorrected))
        print(("Mean photon rate exposure corrected %s day^-1" %(meanRateCorrected*3600.*24)))

        #Calculate bayesian block
        edgesCorrected = bayesian_blocks(evtlistcorrected, fitness='events', p0=self.p0)
        edgesCorrected[0] = timeCorrection[0, 1]
        edgesCorrected[-1] = timeCorrection[-1, 1]

        #Calculate bin event for apperture photometry
        count, tmp = np.histogram(evtlistcorrected, bins=edgesCorrected)
        errcount = np.sqrt(count)

        #Correct edges from exposure
        edges = np.interp(edgesCorrected, timeCorrection[:, 1], timeCorrection[:, 0])
        edges[0] = self.tmin
        edges[-1] = self.tmax

        #Calculate apperture phtometry flux
        flux = np.array(count/(edgesCorrected[1:]-edgesCorrected[:-1]))
        errflux = np.array(errcount/(edgesCorrected[1:]-edgesCorrected[:-1]))

        self.Nbin = len(edges)-1
        self.time_array = np.zeros(self.Nbin*2)
        self.gtifile = []
        for i in range(self.Nbin):
            self.time_array[2*i] = edges[i]
            self.time_array[2*i+1]= edges[i+1]

        self.info("Running LC with "+str(self.Nbin)+" bins")
        for i in range(self.Nbin):
            print(("Bin ",i," Start=",self.time_array[2*i]," Stop=",self.time_array[2*i+1], 'Apperture Photometry=', flux[i], '+/-', errflux[i], 'ph.s^-1'))

        #Dump into ascii
        bbfile = str("%s/BayesianBlocks/%s_bb.dat"%(self.folder,self.srcname))
        np.savetxt(bbfile, np.transpose(np.array([np.array(edges[:-1]), np.array(edges[1:]), np.array(edgesCorrected[1:]-edgesCorrected[:-1]), np.array(count)])),
                   header='tstart tend dt_exposure_corrected count')

        #Load apperture flux point
        time_pt, dTime_pt, flux_pt, errflux_pt = self.readApperturePhotometryPoint()
        
        plt.figure()
        plt.xlabel(r"Time (s)")
        plt.ylabel(r"${\rm Flux\ (photon\ cm^{-2}\ s^{-1})}$")
        
        plot_bayesianblocks(np.array(edges[:-1]), np.array(edges[1:]),
                            flux/surfaceFermi, errflux/surfaceFermi, errflux/surfaceFermi,
                            np.zeros(flux.shape).astype(np.bool))
        plt.errorbar(time_pt, flux_pt/surfaceFermi, yerr=errflux_pt/surfaceFermi, xerr=dTime_pt/2., color='k', ls='None')

        plt.ylim(ymin=max(plt.ylim()[0],np.percentile(flux/surfaceFermi,1)*0.1),
                 ymax=min(plt.ylim()[1],np.percentile(flux/surfaceFermi,99)*2.0))
        plt.xlim(xmin=max(plt.xlim()[0],1.02*min(np.array(edges[:-1]))-0.02*max(np.array(edges[1:]))),
                 xmax=min(plt.xlim()[1],1.02*max(np.array(edges[1:]))-0.02*min(np.array(edges[:-1]))))
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
        plt.setp(mjdaxis.xaxis.get_majorticklabels(), rotation=15)
        plt.tight_layout()

        LcOutPath = self.LCfolder + self.config['target']['name']
        plt.savefig(LcOutPath+"_AP.png", dpi=150, facecolor='w', edgecolor='w',
                    orientation='portrait', format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1)

    def _ManageFolder(self,path):
        """   All files will be stored in a subfolder name path + NLCbin
        Create a subfolder"""
        self.config['out'] = self.LCfolder

    def _MakeLC(self,Path=LightcurvePath) :
        #import gc
        import os
        #gc.enable()
        '''Main function of the Lightcurve script. Read the config file and run the gtlike analysis'''
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')

        self.PrepareLC(self.config['BayesianBlocks']['MakeConfFile'])#Get the config file

        for i in range(self.Nbin):
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
            if (CurConfig['BayesianBlocks']['TSLightCurve']<float(ResultDic.get("TS"))):
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
        if self.config['BayesianBlocks']['DiagnosticPlots'] == 'yes' and len(Npred)>0:
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
                    orientation='portrait', format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1)
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
                    orientation='portrait', format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1)


        if len(Time) > 0:
            plt.figure()
            plt.xlabel(r"Time (s)")
            plt.ylabel(r"${\rm Flux\ (photon\ cm^{-2}\ s^{-1})}$")
            plot_bayesianblocks(Time-TimeErr, Time+TimeErr, Flux, FluxErr, FluxErr, uplim,LcOutPath=LcOutPath)
            plt.ylim(ymin=max(plt.ylim()[0],np.percentile(Flux[~uplim],1)*0.1),
                     ymax=min(plt.ylim()[1],np.percentile(Flux[~uplim],99)*2.0))
            plt.xlim(xmin=max(plt.xlim()[0],1.02*min(Time-TimeErr)-0.02*max(Time+TimeErr)),
                     xmax=min(plt.xlim()[1],1.02*max(Time+TimeErr)-0.02*min(Time-TimeErr)))
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
                    orientation='portrait', format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1)
        else:
            print("[BayesianBlocks] Warning : No valid data")



        if self.config["BayesianBlocks"]["SpectralIndex"] == 0 :
            if len(Time[~uplimIndex]) > 0:
                plt.figure()
                plt.xlabel(r"Time (s)")
                plt.ylabel(r"${\rm Index}$")
                Index = np.asarray(Index)
                IndexErr = np.asarray(IndexErr)
                uplimIndex = uplim #+ Index<0.55
                plot_bayesianblocks(Time[~uplimIndex]-TimeErr[~uplimIndex],
                                    Time[~uplimIndex]+TimeErr[~uplimIndex],
                                    Index[~uplimIndex],
                                    IndexErr[~uplimIndex],
                                    IndexErr[~uplimIndex],
                                    uplimIndex[~uplimIndex])

                plt.ylim(ymin=max(plt.ylim()[0],np.percentile(Index[~uplimIndex],1)*0.1),
                         ymax=min(plt.ylim()[1],np.percentile(Index[~uplimIndex],99)*2.0))
                plt.xlim(xmin=max(plt.xlim()[0],1.02*min(Time-TimeErr)-0.02*max(Time+TimeErr)),
                     	 xmax=min(plt.xlim()[1],1.02*max(Time+TimeErr)-0.02*min(Time-TimeErr)))

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
                    orientation='portrait', format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1)
            else:
               print("[BayesianBlocks] Warning : No valid data")



        #Dump into ascii
        lcfilename = LcOutPath+"_results.dat"
        self.info("Write to Ascii file : "+lcfilename)
        lightcurve.WriteToAscii(Time,TimeErr,Flux,FluxErr,Index,IndexErr,
                     Cutoff,CutoffErr,TS,Npred,lcfilename)

    def readApperturePhotometryPoint(self):
        apfile = str("%s/AppertureLightCurve/TimeExposureCount.txt"%(self.folder))
        time, dTime, exposure, counts = np.loadtxt(apfile, skiprows=1, unpack=True)
        time, dTime, exposure, counts = resampleCount(time, dTime, exposure, counts)
        errcounts = np.sqrt(counts)
        surfaceFermi = 10000 # in cm^2
        flux = counts*surfaceFermi/(exposure)
        errflux = errcounts*surfaceFermi/(exposure)
        time = utils.MJD_to_met(time)
        dTime = utils.MJD_to_met(dTime)-utils.MJD_to_met(0.)
        return time, dTime, flux, errflux

#Functions in order to resample counts on apperture photometry to achieve a given error
def resampleCount(time, dTime, exposure, counts, errorObj=0.1):
    i=1        
    while i < len(time):
        if np.sqrt(float(counts[i-1]))/float(counts[i-1]) > errorObj:
            counts[i-1] += counts[i]
            time[i-1] -= dTime[i-1]/2.
            dTime[i-1] += dTime[i]
            time[i-1] += dTime[i-1]/2.
            exposure[i-1] += exposure[i]
            counts = np.delete(counts, (i), axis=0)
            time = np.delete(time, (i), axis=0)
            dTime = np.delete(dTime, (i), axis=0)
            exposure = np.delete(exposure, (i), axis=0)
        else:
            i += 1
    i = len(time)-1
    if len(time) > 1 and np.sqrt(float(counts[i]))/float(counts[i]) > errorObj:
        counts[i-1] += counts[i]
        time[i-1] -= dTime[i-1]/2.
        dTime[i-1] += dTime[i]
        time[i-1] += dTime[i-1]/2.
        exposure[i-1] += exposure[i]
        counts = np.delete(counts, (i), axis=0)
        time = np.delete(time, (i), axis=0)
        dTime = np.delete(dTime, (i), axis=0)
        exposure = np.delete(exposure, (i), axis=0)
    return time, dTime, exposure, counts
            
