#!/usr/bin/env python
import os
import utils
from submit import call
from config import get_config
import environ
import numpy as np
import ROOT
import pyfits
from gtfunction import Observation

def _log(self, task='', description=''):
    print
    print('# ' + '*' * 60)
    if task:
        task = '%10s --- ' % task
    print('# *** %s%s' %
        ( task, description))
    print '# ' + '*' * 60

def AppLC(infile):
    '''Main function of the apperture photometrie Lightcurve script. Read the config file and run the analysis'''
    ROOT.gROOT.SetBatch(ROOT.kTRUE) #Batch mode

    enricodir = environ.DIRS.get('ENRICO_DIR')
    fermidir = environ.DIRS.get('FERMI_DIR')
    config = get_config(infile)

    folder = config['out']
    #Create a subfolder name LightCurve
    LCoutfolder = folder+"/AppertureLightCurve"
    os.system("mkdir -p "+LCoutfolder)

    #Change the ROI to 1 degree
    config['space']['rad'] = 1

    Nbins = config['AppLC']['NLCbin']#Number of bins
    #Get The time bin
    dt = (config['time']['tmax']-config['time']['tmin'])/Nbins #sec

    Obs = Observation(LCoutfolder, config, config['analysis']['convtype'], tag="")
    if config['AppLC']["FitsGeneration"] == "yes":
        _log('gtselect', 'Select data from library')#run gtselect
        Obs.FirstCut()
        _log('gtmktime', 'Update the GTI and cut data based on ROI')#run gtdiffresp
        Obs.MkTime()

        #Binning from data or using a fix bin size
        if config['AppLC']['binsFromData'] == "no":
            _log('gtbin', 'bin the data into a light-curve using fixe time bin')#run gtbin
            print "Use a dt of %2.2e seconds"%(dt)
            Obs.GtLCbin(dt = dt)
        else:
            spfile=pyfits.open(Obs.eventfile)
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

    spfile = pyfits.open(Obs.eventfile)
    Time = spfile[1].data.field(9)[:-1]

    bounds = ((Time[1:] + Time[:-1])/2.).tolist() # compute the edge of the cells
    bounds.insert(0,Obs.t1)
    bounds.append(Obs.t2)

    print len(Time)
    fil = open(timefile,"w")
    for i in xrange(len(bounds)-1):
        fil.write(str(bounds[i])+'\t'+str(bounds[i+1])+'\n')
    fil.close()

def PlotAppLC(Nbins,LCoutfolder,FITSfile):

    ROOT.gStyle.SetOptStat(0)
    # time MET in MJD
    met_ref = 240106987.-23*60.-6
    mdj_ref= 54689.

    spfile=pyfits.open(FITSfile)

    Time = mdj_ref+(spfile[1].data.field(0)[:-1]-met_ref)/3600./24
    dTime = (spfile[1].data.field(1)[:-1])/3600./24
    Counts = (spfile[1].data.field(2)[:-1])
    Exposure = (spfile[1].data.field(4)[:-1])

    count_histo =ROOT.TH1F("count","count",Nbins,Time[0],Time[-1])
    expo_histo =ROOT.TH1F("exposure","exposure",Nbins,Time[0],Time[-1])
    for i in xrange(len(Time)):
      if Counts[i]>0 and Exposure[i] :
          count_histo.Fill(Time[i],Counts[i])
          expo_histo.Fill(Time[i],Exposure[i])

    flux_histo =ROOT.TH1F()
    count_histo.Copy(flux_histo) #correct for exposure
    flux_histo.Divide(expo_histo)

    for i in xrange(Nbins):#Correct error for exposure
        if expo_histo.GetBinContent(i+1)>0:
            flux_histo.SetBinError(i+1,count_histo.GetBinError(i+1)/expo_histo.GetBinContent(i+1))

    ######################################################################################
    #Save event time exposure and count
    file_evt = open(LCoutfolder+'/TimeExposureCount.txt',"w")

    #Write into file
    file_evt.write("Time\tdTime\tExposure\tCounts\n")
    for i in xrange(len(Time)):
       file_evt.write(str(Time[i])+"\t"+str(dTime[i])+"\t"+str(Exposure[i])+"\t"+str(Counts[i])+"\n")
    file_evt.close()
    ######################################################################################

    file_lc = open(LCoutfolder+'/AppLC.txt',"w")
    for i in xrange(Nbins):
        file_lc.write(str(flux_histo.GetBinCenter(i))+"\t"+str(flux_histo.GetBinContent(i))+"\t"+str(expo_histo.GetBinContent(i))+"\n")
    file_lc.close()

    CanvCount = ROOT.TCanvas()
    count_histo.SetMarkerStyle(20)
    count_histo.SetLineColor(1)
    count_histo.SetMarkerColor(1)
    count_histo.SetXTitle("Time (MJD)")
    count_histo.Draw("ep")

   #Save the canvas in the Apperture LightCurve subfolder
    CanvCount.Print(LCoutfolder+'/Counts.eps')
    CanvCount.Print(LCoutfolder+'/Counts.C')

    CanvExposure = ROOT.TCanvas()
    CanvExposure.SetGridx()
    CanvExposure.SetGridy()
    expo_histo.SetMarkerStyle(20)
    expo_histo.SetLineColor(1)
    expo_histo.SetMarkerColor(1)
    expo_histo.SetXTitle("Time (MJD)")
    expo_histo.Draw("ep")

   #Save the canvas in the Apperture LightCurve subfolder
    CanvExposure.Print(LCoutfolder+'/Exposure.eps')
    CanvExposure.Print(LCoutfolder+'/Exposure.C')

    CanvFlux = ROOT.TCanvas()
    flux_histo.SetNameTitle("Apperture_Photometry_Flux","Apperture_Photometry_Flux")
    flux_histo.SetYTitle("Flux")
    flux_histo.SetXTitle("Time (MJD)")
    flux_histo.SetMarkerStyle(20)
    flux_histo.SetLineColor(1)
    flux_histo.SetMarkerColor(1)
    flux_histo.Draw("ep")

   #Save the canvas in the Apperture LightCurve subfolder
    CanvFlux.Print(LCoutfolder+'/AppLC.eps')
    CanvFlux.Print(LCoutfolder+'/AppLC.C')

