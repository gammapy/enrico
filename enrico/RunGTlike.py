#!/usr/bin/env python
import os
import SummedLikelihood
import utils
from submit import call
from config import get_config
import environ


def GenAnalysisObject(config, verbose = 1):
    #check is the summed likelihood method should be used and get the 
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    SummedLike = config['Spectrum']['SummedLike']
    folder = config['out']
    if SummedLike == 'yes':
        # Create two obs instances
        runfitfront = utils.Analysis(folder, config, tag="FRONT", convtyp=0, verbose = verbose)
        runfitback = utils.Analysis(folder, config, tag="BACK", convtyp=1, verbose = verbose)
        FitB = runfitback.CreateFit()
        FitF = runfitfront.CreateFit()
        Fit = SummedLikelihood.SummedLikelihood()
        Fit.addComponent(FitB)
        Fit.addComponent(FitF)
        runfit = runfitback
    else:
        convtype = config['analysis']['convtype']
        # Create one obs instance
        runfit = utils.Analysis(folder, config, tag="", convtyp=convtype, verbose = verbose)
        Fit = runfit.CreateFit()
    return runfit,Fit

def run(infile):
    """Run an entire Fermi analysis (spectrum) by reading a config file"""
    config = get_config(infile)
    folder = config['out']
    os.system('mkdir -p ' + folder)

    runfit,Fit = GenAnalysisObject(config)
    # create all the fit files and run gtlike
    runfit.PerformFit(Fit)

    if config['verbose'] == 'yes' :
        utils.GetFluxes(Fit,runfit.obs.Emin,runfit.obs.Emax) #print the flux of all the sources

    Result = runfit.GetAndPrintResults(Fit)#Get and dump the target specific results
    utils.DumpResult(Result, config)

    #plot the SED and model map if possible and asked
    if config['Spectrum']['ResultPlots'] == 'yes' :
    	os.system("mkdir -p "+config['out'] + '/Spectrum/')
        if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
            runfit.ComputeSED(Fit)
        outXml = utils._dump_xml(config)
        if SummedLike == 'yes': # the possiblity of making the model map is checked inside the function
            runfitback.ModelMap(outXml)
            runfitfront.ModelMap(outXml)
        else:
            runfit.ModelMap(outXml)

    #  Make energy bins by run a *new* analysis
    Nbin = config['Ebin']['NumEnergyBins']
    if int(Nbin) > 0:
        configfiles = utils.PrepareEbin(Fit, runfit)
        ind = 0
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')
        for conf in configfiles:
             pathconf = folder + "/Ebin" + str(Nbin) +"/" + conf
             Newconfig = get_config(pathconf)
             cmd = enricodir+"/enrico/RunGTlike.py "+pathconf
             if Newconfig['Ebin']['Submit'] == 'no' : #run directly
                 os.system(cmd)
             else : #submit a job to a cluster
                 prefix = Newconfig['out'] + "/Ebin" + str(ind) 
                 scriptname = prefix + "_Script.sh"
                 JobLog = prefix + "_Job.log"
                 JobName = (Newconfig['target']['name'] + "_" +
                           Newconfig['analysis']['likelihood'] +
                           "_Ebin_" + str(ind) + "_" + Newconfig['file']['tag'])
                 call(cmd, enricodir, fermidir, scriptname, JobLog, JobName)# submition
             ind+=1

# @todo: Should this be a command line utility in bin?
if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
        config = get_config(infile)
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)
    run(infile)
