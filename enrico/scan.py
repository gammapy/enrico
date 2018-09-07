"""
Make a profile likelihood for each free parameters of the source of interest
begun November 2013
"""
import enrico.constants as cst
import RunGTlike
import ROOT
import numpy,os,string,array
from enrico import Loggin
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt

def MakeScan(Fit,spectrum,par,bmin,bmax,opt,N=100):
    Param = numpy.zeros(N)
    loglike = numpy.zeros(N)
    for i in xrange(N):
        Param[i] = bmin + (bmax-bmin)*i/(N-1.)
        spectrum.getParam(par).setValue(Param[i])
#        spectrum.getParam(par).setFree(0)
#        loglike[i] = Fit.fit(0,covar=False,optimizer=opt)
#        Fit.optimize(0) 
        Fit.logLike.syncParams()
        loglike[i] = -Fit.logLike.value()

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

          plt.figure()
          plt.plot(param,loglike,"-r")
          plt.title(par)
          plt.xlabel("Parameter: "+par)
          plt.ylabel("Log(Like)")

          os.system("mkdir -p "+config["out"]+"/"+cst.ScanPath)
          savefile = open(config["out"]+"/"+cst.ScanPath+ "/Scan_"+par+".txt","w")
          for i in xrange(param.size):
             savefile.write(str(param[i])+" "+str(loglike[i])+"\n")
          savefile.close()
          plt.savefig(config["out"]+"/"+cst.ScanPath+ "/Scan_"+par+".png", dpi=150, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

def Contour(config):
    # ROOT.gROOT.SetBatch(ROOT.kTRUE)
#    cres = ROOT.TCanvas("Contour")
    config["Spectrum"]["FitsGeneration"] = "no"
    parname1 = config["Contours"]["parname1"]
    parname2 = config["Contours"]["parname2"]

    FitRunner,Fit = RunGTlike.GenAnalysisObjects(config)
    spectrum = Fit[FitRunner.obs.srcname].funcs['Spectrum']
    
    ParName = spectrum.paramNames

    mes = Loggin.Message()
    mes.info("Computing Contours for "+parname1+" and "+parname2)

    ### Check part !!!!
    findpar2 = findpar1 = False
    for par in ParName : #Loop over the parameters to check
        if par == parname1:
            findpar1 = True
            if  not(spectrum.getParam(par).isFree()):
                mes.error(parname1+" is not a free parameter")
        if par == parname2:
            findpar2 = True
            if  not(spectrum.getParam(par).isFree()):
                mes.error(parname2+" is not a free parameter")

    if not(findpar1):
        mes.error(parname1+" is not a valid parameter")
    if not(findpar2):
        mes.error(parname2+" is not a valid parameter")

    bestloglike = Fit.fit(0,covar=False,optimizer=config['fitting']['optimizer'])
    print spectrum
    print "Min LogLikelihood =",bestloglike

    ## get values
    ParValue1 = spectrum.getParam(parname1).value()
    ParError1 = spectrum.getParam(parname1).error()
    bmin1,bmax1 = spectrum.getParam(parname1).getBounds()

    bmin1 = max(bmin1,ParValue1-20*ParError1)
    bmax1 = min(bmax1,ParValue1+20*ParError1)

    ParValue2 = spectrum.getParam(parname2).value()
    ParError2 = spectrum.getParam(parname2).error()
    bmin2,bmax2 = spectrum.getParam(parname2).getBounds()

    bmin2 = max(bmin2,ParValue2-20*ParError2)
    bmax2 = min(bmax2,ParValue2+20*ParError2)

    N = 100
    param2  = numpy.zeros(N)
    loglike = ROOT.TH2F("loglike","Contours (68%, 95%, 99%)",N,bmin1,bmax1,N,bmin2,bmax2)
    spectrum.getParam(parname2).setFree(0)

    mes.info("Boundaries for "+parname1+" ["+str(bmin1)+","+str(bmax1)+"]")
    mes.info("Boundaries for "+parname2+" ["+str(bmin2)+","+str(bmax2)+"]")

    for i in xrange(N):

      param2[i] = bmin2 + (bmax2-bmin2)*i/(N-1.)

      spectrum.getParam(parname2).setValue(param2[i])

      param1,ll = MakeScan(Fit,spectrum,parname1,bmin1,bmax1,config['fitting']['optimizer'],N)
      
      for j in xrange(N):
         loglike.Fill(param1[j],param2[i],ll[j])

    os.system("mkdir -p "+config["out"]+"/"+cst.ScanPath)
    cres = ROOT.TCanvas("Contours")
    loglike.SetMinimum(bestloglike);
    loglike.SetMaximum(bestloglike+3);
    loglike.SetXTitle(parname1);
    loglike.SetYTitle(parname2);

    loglike.SetStats(000)
    loglike.SetContour(3)
    loglike.SetContourLevel(0,bestloglike+0.5)
    loglike.SetContourLevel(1,bestloglike+4./2.)
    loglike.SetContourLevel(2,bestloglike+6.63/2.)
    loglike.Draw("CONT1");

    tgrres = ROOT.TGraphErrors(2,array.array('f',[ParValue1,ParValue1]),array.array('f',[ParValue2,ParValue2]),array.array('f',[ParError1,ParError1]),array.array('f',[ParError2,ParError2]))
    tgrres.Draw(".pz")
    cres.Print(config["out"]+"/"+cst.ScanPath+ "/Contours_"+parname1+"_"+parname2+".eps")
    cres.Print(config["out"]+"/"+cst.ScanPath+ "/Contours_"+parname1+"_"+parname2+".C")
    cres.Print(config["out"]+"/"+cst.ScanPath+ "/Contours_"+parname1+"_"+parname2+".png")



    mes.success("Scan Performed")
