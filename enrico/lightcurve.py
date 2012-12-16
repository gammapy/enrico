import os
from os.path import join
import numpy as np
import ROOT
import utils
import root_style
import plotting
import environ
from config import get_config
from submit import call
from enrico.RunGTlike import run

def PrepareLC(infile,write = 'no'):
    """Simple function to prepare the LC generation : generate and write the config files"""
    config = get_config(infile)

    #Read the config
    Tag = config['file']['tag']
    Nbin = config['LightCurve']['NLCbin']
    tmin = config['time']['tmin']
    tmax = config['time']['tmax']

    # One point of the LC will be computed as a spectrum plot. enrico_fit will be used
    # Do fits files will be generated
    config['Spectrum']['FitsGeneration'] = config['LightCurve']['FitsGeneration']

    #Froze the Spectral index at a value of 2
    config['Spectrum']['FrozenSpectralIndex'] = 2


    #TS limit. Compute an UL if the TS is below TSLightCurve
    config['UpperLimit']['TSlimit'] = config['LightCurve']['TSLightCurve']

    #All files will be stored in a subfolder name LightCurve + NLCbin
    config['out'] +='/LightCurve_'+str(config['LightCurve']['NLCbin'])+'bins'

    config['verbose'] ='no' #Be quiet

    #No plot, no bin in energy, Normal UL
    config['Spectrum']['ResultPlots'] = 'no'
    config['Ebin']['NumEnergyBins'] = 0
    config['UpperLimit']['envelope'] = 'no'
    #No submition
    config['Spectrum']['Submit'] = 'no'

    AllConfigFile = []#All the config file in the disk are stored in a list
    for i in xrange(Nbin):
        dt = (tmax - tmin) / Nbin
        config['time']['tmin'] = tmin + i * dt
        config['time']['tmax'] = tmin + (i + 1) * dt
        config['file']['tag'] = Tag + '_LC_' + str(i)
        #Name of the config file
        filename = (config['out'] + "/Config_" +
                    str(config['time']['tmin']) + "_" +
                    str(config['time']['tmax']))

        if write == 'yes':
            config.write(open(filename, 'w'))

        AllConfigFile.append(filename)
    #Return the list of config
    return AllConfigFile


def WriteToAscii(Time, TimeErr, Flux, FluxErr, TS, Npred, filename):
    """Write the results of the LC in a Ascii file"""
    flc = open(filename, 'w')
    flc.write('Time (MET) Delta_Time Flux(ph cm-2 s-1) '
              'Delta_Flux TS Npred\n')
    for i in xrange(len(Time)):
        flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                  str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                  str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    flc.close()


def MakeLC(infile) :
    '''Main function of the Lightcurve script. Read the config file and run the gtlike analysis'''
    ROOT.gROOT.SetBatch(ROOT.kTRUE) #Batch mode

    enricodir = environ.DIRS.get('ENRICO_DIR')
    fermidir = environ.DIRS.get('FERMI_DIR')
    config = get_config(infile)

    folder = config['out']
    #Create a subfolder name LightCurve
    LCoutfolder = folder+"/LightCurve_"+str(config['LightCurve']['NLCbin'])+"bins/"
    os.system("mkdir -p "+LCoutfolder)

    AllConfigFile = PrepareLC(infile, config['LightCurve']['MakeConfFile'])#Get the config file

    for i in xrange(config['LightCurve']['NLCbin']):
        if config['LightCurve']['Submit'] == 'yes':
            cmd = "enrico_sed "+AllConfigFile[i]
            scriptname=LCoutfolder+"/LC_Script_"+str(i)+".sh"
            JobLog = LCoutfolder+"LC_Job_"+str(i)+".log"
            JobName = (config['target']['name'] + "_" +
                   config['analysis']['likelihood'] +
                   "_LC_" + config['file']['tag'])+"_"+str(i)+".log"

            call(cmd,enricodir,fermidir,scriptname,JobLog,JobName)#Submit the job
        else :
            run(AllConfigFile[i])#run in command line

def PlotLC(infile):
    '''Plot a lightcurve which have been generated previously'''
    root_style.RootStyle()#Nice plot style
    config = get_config(infile) # read config file

    folder = config['out']
    print "Reading files produced by enrico"
    Nbin = config['LightCurve']['NLCbin']

    LcOutPath = folder+"/LightCurve_"+str(config['LightCurve']['NLCbin'])+"bins/"+ config['target']['name']

    #Result are stored into list. This allow to get rid of the bin which failled
    Time = []
    TimeErr = []
    Flux = []
    FluxErr = []
    Npred = []
    TS = []

    AllConfigFile = PrepareLC(infile)#Get the config file

    for i in xrange(Nbin):
        CurConfig = get_config(AllConfigFile[i])
     #Read the result. If it fails, it means that the bins has not bin computed. A warning message is printed
        try :
            ResultDic = utils.ReadResult(CurConfig)
        except :
            print "WARNING : fail reading the config file : ",AllConfigFile[i]
            print "Job Number : ",i
            print "Please have a look at this job log file"
            continue

        #Update the time and time error array
        Time.append((ResultDic.get("tmax")+ResultDic.get("tmin"))/2.)
        TimeErr.append((ResultDic.get("tmax")-ResultDic.get("tmin"))/2.)
        #Check is an ul have been computed. The error is set to zero for the TGraph.
        if ResultDic.has_key('Ulvalue') :
            Flux.append(ResultDic.get("Ulvalue"))
            print ResultDic.get("Ulvalue")
            FluxErr.append(0)
        else :
            Flux.append(ResultDic.get("Flux"))
            FluxErr.append(ResultDic.get("dFlux"))
        #Get the Npred and TS values
        Npred.append(ResultDic.get("Npred"))
        TS.append(ResultDic.get("TS"))

    #change the list into np array
    TS = np.array(TS)
    Npred = np.array(Npred)
    Time = np.array(Time)
    TimeErr = np.array(TimeErr)
    Flux = np.array(Flux)
    FluxErr = np.array(FluxErr)

    #Plots the diagnostic plots is asked
    # Plots are : Npred vs flux
    #             TS vs Time
    if config['LightCurve']['DiagnosticPlots'] == 'yes':
        gTHNpred,TgrNpred = plotting.PlotNpred(Npred,Flux,FluxErr)
        CanvNpred = ROOT.TCanvas()
        CanvNpred.SetGridx()
        CanvNpred.SetGridy()
        gTHNpred.Draw()
        TgrNpred.Draw('zP')
        CanvNpred.Print(LcOutPath+"_Npred.eps")
        CanvNpred.Print(LcOutPath+"_Npred.C")

        gTHTS,TgrTS = plotting.PlotTS(Time,TimeErr,TS)
        CanvTS = ROOT.TCanvas()
        gTHTS.Draw()
        TgrTS.Draw('zP')
        CanvTS.Print(LcOutPath+'_TS.eps')
        CanvTS.Print(LcOutPath+'_TS.C')

#    Plot the LC itself. This function return a TH2F for a nice plot
#    a TGraph and a list of TArrow for the ULs
    gTHLC,TgrLC,ArrowLC = plotting.PlotLC(Time,TimeErr,Flux,FluxErr)
    CanvLC = ROOT.TCanvas()
    gTHLC.Draw()
    TgrLC.Draw('zP')

    #plot the ul as arrow
    for i in xrange(len(ArrowLC)):
        ArrowLC[i].Draw()

    #Save the canvas in the LightCurve subfolder
    CanvLC.Print(LcOutPath+'_LC.eps')
    CanvLC.Print(LcOutPath+'_LC.C')

    #Dump into ascii
    lcfilename = LcOutPath+"_results.dat"
    print "Write to Ascii file : ",lcfilename
    WriteToAscii(Time,TimeErr,Flux,FluxErr,TS,Npred,lcfilename)
