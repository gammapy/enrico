#!/usr/bin/env python
import os
import numpy as np
import astropy.io.fits as fits
from enrico.constants import  DAY_IN_SECOND, AppLCPath #met_ref, mdj_ref,
from enrico.gtfunction import Observation
from enrico.config import get_config
from enrico import environ
from enrico import utils
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt

def AppLC(infile):
    '''Main function of the apperture photometrie Lightcurve script. Read the config file and run the analysis'''

    enricodir = environ.DIRS.get('ENRICO_DIR')
    fermidir = environ.DIRS.get('FERMI_DIR')
    config = get_config(infile)

    folder = config['out']
    #Create a subfolder name LightCurve
    LCoutfolder = folder+"/"+AppLCPath
    utils.mkdir_p(LCoutfolder)

    #Change the ROI to the desired radius in degree, legacy 1 deg.
    try: config['space']['rad'] = config['AppLC']['rad']
    except NameError: config['space']['rad'] = 1

    #Change the minimum energy, legacy std low energy cut.
    try:
        if config['AppLC']['emin'] != -1:
            config['energy']['emin'] = config['AppLC']['emin']
    except NameError: pass
    #Change the maximum energy, legacy std high energy cut.
    try:
        if config['AppLC']['emax'] != -1:
            config['energy']['emax'] = config['AppLC']['emax']
    except NameError: pass

    Nbins = config['AppLC']['NLCbin']#Number of bins
    #Get The time bin
    dt = (config['time']['tmax']-config['time']['tmin'])/Nbins #sec

    Obs = Observation(LCoutfolder, config, tag="")
    if config['AppLC']["FitsGeneration"] == "yes":
        _log('gtselect', 'Select data from library')#run gtselect
        Obs.FirstCut()
        Obs.SelectEvents()
        _log('gtmktime', 'Update the GTI and cut data based on ROI')#run gtdiffresp
        Obs.MkTime()

        #Binning from data or using a fix bin size
        if config['AppLC']['binsFromData'] == "no":
            _log('gtbin', 'bin the data into a light-curve using fixe time bin')#run gtbin
            print(("Use a dt of %2.2e seconds"%(dt)))
            Obs.GtLCbin(dt = dt)
        else:
            spfile=fits.open(Obs.eventfile)
            diff = spfile[1].data.field(9)[1:-1]-spfile[1].data.field(9)[:-2]
            dt = np.min(diff)/2.  ##Compute the delta T as being the min delta t between 2 events divided by 2
            timefile = LCoutfolder+"/Timebin.txt"
            MakeTimebinFile(Obs,timefile)
            _log('gtbindef', 'define de bins')#run gtbindef
            Obs.GtBinDef(timefile)
            _log('gtbin', 'bin the data into a light-curve using bins based on data')#run gtbin
            Obs.GtLCbin(dt = 0)

        _log('gtexposure', 'compute the exposure')#run gtexposure
        Obs.GtExposure()

    #Get Some usefull value here. This allow PlotAppLC to be call independently
    Nbins = config['AppLC']['NLCbin']#Number of bins
    #Plot the results and dump into ascii files
    PlotAppLC(Nbins,LCoutfolder,Obs.lcfile)

def MakeTimebinFile(Obs,timefile):
    spfile = fits.open(Obs.eventfile)
    spfile[1].data.sort(order='TIME')
    Time = spfile[1].data.field(9)[:-1]

    bounds = ((Time[1:] + Time[:-1])/2.).tolist() # compute the edge of the cells
    bounds.insert(0,Obs.t1)
    bounds.append(Obs.t2)

    print((len(Time)))
    fil = open(timefile,"w")
    for i in range(len(bounds)-1):
        fil.write(str(bounds[i])+'\t'+str(bounds[i+1])+'\n')
    fil.close()

def PlotAppLC(Nbins,LCoutfolder,FITSfile):

    spfile=fits.open(FITSfile)
    spfile[1].data.sort(order='TIME')

    Time =  utils.met_to_MJD(spfile[1].data.field(0)[:-1])#mdj_ref+(spfile[1].data.field(0)[:-1]-met_ref)/DAY_IN_SECOND
    dTime = (spfile[1].data.field(1)[:-1])/DAY_IN_SECOND
    Counts = (spfile[1].data.field(2)[:-1])
    Exposure = (spfile[1].data.field(4)[:-1])

    count_histo,count_edges = np.histogram(Time,Nbins,weights=Counts)
    expo_histo,expo_edges = np.histogram(Time,Nbins,weights=Exposure)
    flux_histo,flux_edges = np.histogram(Time,Nbins,weights=Counts/Exposure)
    dCounts = np.sqrt(Counts)
    dflux_histo,_ = np.histogram(Time,Nbins,weights=dCounts/Exposure)

    ######################################################################################
    #Save event time exposure and count
    file_evt = open(LCoutfolder+'/TimeExposureCount.txt',"w")

    #Write into file
    file_evt.write("Time\tdTime\tExposure\tCounts\n")
    for i in range(len(Time)):
       file_evt.write(str(Time[i])+"\t"+str(dTime[i])+"\t"+str(Exposure[i])+"\t"+str(Counts[i])+"\n")
    file_evt.close()
    ######################################################################################

    file_lc = open(LCoutfolder+'/AppLC.txt',"w")

    for i in range(len(Time)):
        file_lc.write(str(Time[i])+"\t"+str(Counts[i]/Exposure[i])+"\t"+str(dCounts[i]/Exposure[i])+"\n")
    file_lc.close()

    plt.figure()
    print(count_edges)
    print(count_histo)
    plt.errorbar((count_edges[:1]+count_edges[1:])/2.,count_histo,fmt="o")
    plt.xlabel("Time (MJD)")
   #Save the canvas in the Apperture LightCurve subfolder
    plt.savefig(LCoutfolder+'/Counts.png', dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)


    plt.figure()
    plt.errorbar((expo_edges[:1]+expo_edges[1:])/2.,expo_histo,fmt="o")
    plt.xlabel("Time (MJD)")
   #Save the canvas in the Apperture LightCurve subfolder
    plt.savefig(LCoutfolder+'/Exposure.png', dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

    plt.figure()
    plt.errorbar((flux_edges[:1]+flux_edges[1:])/2.,flux_histo,yerr=dflux_histo,fmt="o")
    plt.ylabel("Flux")
    plt.xlabel("Time (MJD)")
   #Save the canvas in the Apperture LightCurve subfolder
    plt.savefig(LCoutfolder+'/AppLC.png', dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

def _log(task='', description=''):
    print()
    print(("\033[34m"+'# ' + '*' * 60))
    if task:
        task = '%10s --- ' % task
    print(("\033[34m"+'# *** %s%s' %
        ( task, description)))
    print(("\033[34m"+'# ' + '*' * 60+"\033[0m"))


if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
    except:
        from enrico import Loggin
        mes = Loggin.Message()
        mes.error('Config file not found.')
    AppLC(infile)
