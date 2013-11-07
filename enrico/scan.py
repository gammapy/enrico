"""
Make a profile likelihood for each free parameters of the source of interest
begun November 2013
"""
import enrico.constants as cst
import RunGTlike
import ROOT
import numpy,os

def MakeScan(Fit,spectrum,par,bmin,bmax,opt):
    N=30
    Param = numpy.zeros(N)
    loglike = numpy.zeros(N)
    for i in xrange(N):
        Param[i] = bmin + (bmax-bmin)*i/(N-1.)
        spectrum.getParam(par).setValue(Param[i])
        spectrum.getParam(par).setFree(0)
        loglike[i] = Fit.fit(0,covar=False,optimizer=opt)
    
    return Param,loglike

def Scan(config):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    cres = ROOT.TCanvas("Scan")
    config["Spectrum"]["FitsGeneration"] = "no"
    FitRunner,Fit = RunGTlike.GenAnalysisObjects(config)
    spectrum = Fit[FitRunner.obs.srcname].funcs['Spectrum']
    ParName = spectrum.paramNames

    Fit.fit(0,covar=False,optimizer=config['fitting']['optimizer'])

    for par in ParName : #Loop over the parameters and get value, error and scale
      if  spectrum.getParam(par).isFree():
          print "Scan for parameter ",par
          ParValue = spectrum.getParam(par).value()
          ParError = spectrum.getParam(par).error()
          bmin,bmax = spectrum.getParam(par).getBounds()

          bmin = max(bmin,ParValue-10*ParError)
          bmax = min(bmax,ParValue+10*ParError)

          param,loglike = MakeScan(Fit,spectrum,par,bmin,bmax,config['fitting']['optimizer'])

         #restore best fit parameters
          spectrum.getParam(par).setFree(1)
          ParValue = spectrum.getParam(par).setValue(ParValue)

          tgr = ROOT.TGraph(param.size,param,loglike)
          tgr.GetHistogram().SetXTitle("Parameter: "+par)
          tgr.GetHistogram().SetYTitle("Log(Like)")
          tgr.SetLineColor(2)
          tgr.Draw('AL')

          os.system("mkdir -p "+config["out"]+"/"+cst.ScanPath)
          cres.Print(config["out"]+"/"+cst.ScanPath+ "/Scan_"+par+".eps")
          cres.Print(config["out"]+"/"+cst.ScanPath+ "/Scan_"+par+".C")
          cres.Print(config["out"]+"/"+cst.ScanPath+ "/Scan_"+par+".png")

