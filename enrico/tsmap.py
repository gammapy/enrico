#!/usr/bin/env python
import os
import pyfits
import SummedLikelihood
import utils
from submit import call
from config import get_config
import environ
from gtfunction import Observation

def _PixelFile(folder,i,j):
    return folder+'/TSMap/Pixel_'+str(i)+'_'+str(j)

def GetSrc(Fit,ra,dec):
    ind = 0
    for name in  Fit.logLike.srcNames() :
        src = Fit.logLike.getSource(name).clone()
        if str(type(src))=="<type 'NoneType'>" :
            continue
        if src.getType() == 'Point' :
            src.setName("Spurious")
            src.setSpectrum("PowerLaw")
            src.getSrcFuncs()['Position'].getParam('RA').setValue(ra)
            src.getSrcFuncs()['Position'].getParam('DEC').setValue(dec)

            src.getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-5,1e5)
            src.getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-9)
            src.getSrcFuncs()['Spectrum'].getParam('Index').setBounds(-5,0)

            src.getSrcFuncs()['Spectrum'].getParam('Scale').setValue(300)
            src.getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(1e-5,1e5)

            return src


def FitOnePixel(infile,ra,dec,i,j) :
    """@todo: document me"""
    config = get_config(infile)
    outXml = utils._dump_xml(config)

    folder = config['out']

    if config['Spectrum']['SummedLike'] == 'yes':
        # Create two obs instances
        runfitfront, _ = utils.Analysis(folder, config, tag="FRONT", convtyp=0)
        runfitback, _ = utils.Analysis(folder, config, tag="BACK", convtyp=1)
        runfitfront.obs.xmlfile = outXml
        runfitback.obs.xmlfile = outXml
        FitB = runfitback.CreateFit()
        FitF = runfitfront.CreateFit()
        Fit = SummedLikelihood.SummedLikelihood()
        Fit.addComponent(FitB)
        Fit.addComponent(FitF)
        runfit = runfitback

    else:
        convtype = config['analysis']['convtype']
        # Create one obs instanceFit.addSource
        runfit, _ = utils.Analysis(folder, config, tag="", convtyp=convtype)
        runfit.obs.xmlfile = outXml

        Fit = runfit.CreateFit()

    src = GetSrc(Fit,ra,dec)

    if config['TSMap']['RemoveTarget'] :
        Fit.deleteSource(config['target']['name'])

    if config['TSMap']['Re-Fit']:
       Fit.fit(0,optimizer=config['fitting']['optimizer'])

    for par in xrange(Fit.logLike.getNumParams()):
        Fit[par].setFree(0)

    Fit.addSource(src)

    Fit.writeXml(folder+"/TSMap/model_"+str(ra)+"_"+str(dec)+".xml")
    runfit.obs.xmlfile = folder+"/TSMap/model_"+str(ra)+"_"+str(dec)+".xml"

    TSFit = runfit.CreateFit()
    TSFit.fit(0,optimizer=config['fitting']['optimizer'])
    TS = TSFit.Ts("Spurious")
    loglike = TSFit.logLike.value()

    fsave = open(_PixelFile(folder,i,j),'w')
    fsave.write(str(ra)+"\t"+str(dec)+"\t"+str(TS)+"\t"+str(loglike))
    fsave.close()

def runTSMap(infile) :
    config = get_config(infile)
    
    folder = config['out']
    os.system('mkdir -p ' + folder+'/TSMap')
    runfit = Observation(folder, config)

    cmap = pyfits.open(runfit.cmapfile)

    npix_im = min(cmap[0].header['NAXIS1'],cmap[0].header['NAXIS2'])
    npix = min(config['TSMap']['npix'],npix_im)
    RAref = cmap[0].header['CRVAL1']
    DECref = cmap[0].header['CRVAL2']

    binsz = cmap[0].header['CDELT1']

    enricodir = environ.DIRS.get('ENRICO_DIR')
    prog = enricodir+"/enrico/tsmap.py "+os.getcwd()+"/"+infile

    for i in xrange(npix):
        for j in xrange(npix):
             ra = RAref + binsz*(i-npix/2.)
             dec = DECref + binsz*(j-npix/2.)
             print 'Run Pixel evaluation at ',ra,' ',dec

             cmd = prog +" "+ str(ra) +" "+ str(dec) +" "+ str(i) +" "+ str(j)
             prefix = config['out'] + "/TSMap/TSMap_" + str(i) +"_"+ str(j)
             scriptname = prefix + "_Script.sh"
             JobLog = prefix + "_Job.log"
             JobName = (config['target']['name'] + "_TSMap_" + str(i) +"_"+ str(j))
             call(cmd, enricodir, scriptname, JobLog, JobName)


def PlotTSmap(infile) :
    config = get_config(infile)
    
    folder = config['out']
    runfit = Observation(folder, config)

    header = pyfits.getheader(runfit.cmapfile)
    data = pyfits.getdata(runfit.cmapfile)*0

    npix_im = min(header['NAXIS1'],header['NAXIS2'])
    npix = min(config['TSMap']['npix'],npix_im)

    Xref = header['CRPIX1']
    Yref = header['CRPIX2']

    binsz = header['CDELT1']

    import string
    for i in xrange(npix):
        for j in xrange(npix):
            try : 
                lines = open(_PixelFile(folder,i,j),"r").readlines()
                Value = float(string.split(lines[0])[2])
            except :
                print "Cannot find, open or read ",_PixelFile(folder,i,j)
                Value = 0.
            data[Xref+ (i-npix/2.)][Yref+ (j-npix/2.)] = Value

    pyfits.writeto(folder+"/TSMap.fits",data,header)
    print "TS Map saved in "+folder+"/TSMap.fits"

if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)

    if len(sys.argv)==2 :
         print "TO BE IMPLEMENTED"
    else :
        FitOnePixel(infile,float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]))
