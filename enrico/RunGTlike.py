#!/usr/bin/env python
import os
import SummedLikelihood
import utils
from submit import call
from config import get_config
import environ

def run(infile):
    """Run an entire Fermi analysis (spectrum) by reading a config file"""
    config = get_config(infile)
    folder = config['out']
    os.system('mkdir -p ' + folder)

    #check is the summed likelihood method should be used and get the 
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    SummedLike = config['Spectrum']['SummedLike']
    if SummedLike == 'yes':
        # Create two obs instances
        runfitfront, _ = utils.Analysis(folder, config, tag="FRONT", convtyp=0)
        runfitback, _ = utils.Analysis(folder, config, tag="BACK", convtyp=1)
        FitB = runfitback.CreateFit()
        FitF = runfitfront.CreateFit()
        Fit = SummedLikelihood.SummedLikelihood()
        Fit.addComponent(FitB)
        Fit.addComponent(FitF)
        runfit = runfitback
    else:
        convtype = config['analysis']['convtype']
        # Create one obs instance
        runfit, _ = utils.Analysis(folder, config, tag="", convtyp=convtype)
        Fit = runfit.CreateFit()
    # create all the fit files and run gtlike
    runfit.PerformFit(Fit)

    utils.GetFluxes(Fit,runfit.obs.Emin,runfit.obs.Emax) #print the flux of all the sources

    Result = runfit.GetAndPrintResults(Fit)#Get and dump the target specific results
    utils.DumpResult(Result, config)

    #plot the SED and model map if possible and asked
    if config['Spectrum']['ResultPlots'] == 'yes' :
        if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
            runfit.ComputeSED(Fit)
        outXml = utils._dump_xml(config)
        if SummedLike == 'yes': # the possiblity of making the model map is checked inside the function
            runfitback.ModelMap(outXml)
            runfitfront.ModelMap(outXml)
        else:
            runfit.ModelMap(outXml)

    #  Make energy bins by run a *new* analysis
    if int(config['Ebin']['NumEnergyBins']) > 0:
        configfiles = utils.PrepareEbin(Fit, runfit)
        ind = 0
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')
        for conf in configfiles:
             config = get_config(conf)
             cmd = "enrico_sed " + conf
             if config['Ebin']['Submit'] == 'no' : #run directly
                 os.system(cmd)
             else : #submit a job to a cluster
                 prefix = config['out'] + "/Ebin" + str(ind) 
                 scriptname = prefix + "_Script.sh"
                 JobLog = prefix + "_Job.log"
                 JobName = (config['target']['name'] + "_" +
                           config['analysis']['likelihood'] +
                           "_Ebin_" + str(ind) + "_" + config['file']['tag'])
                 call(cmd, enricodir, fermidir, scriptname, JobLog, JobName)# submition
             ind+=1

# @todo: Should this be a command line utility in bin?
if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)
    run(infile)
