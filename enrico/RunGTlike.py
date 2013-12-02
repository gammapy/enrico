#!/usr/bin/env python
import os
from enrico import utils
from enrico import energybin
from enrico.config import get_config
from enrico.gtfunction import Observation
from enrico.fitmaker import FitMaker


def Analysis(folder, config, tag="", convtyp='-1', verbose = 1):
    """ run an analysis"""
    Obs = Observation(folder, config, convtyp, tag=tag)
    if verbose:
        utils._log('SUMMARY: ' + tag)
        Obs.printSum()
    FitRunner = FitMaker(Obs, config)##Class
    if config['Spectrum']['FitsGeneration'] == 'yes':
        FitRunner.GenerateFits() #Generates fits files
    return FitRunner

def GenAnalysisObjects(config, verbose = 1, xmlfile =""):
    #check is the summed likelihood method should be used and get the 
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    SummedLike = config['Spectrum']['SummedLike']
    folder = config['out']
    if SummedLike == 'yes':
        # Create two obs instances
        FitRunnerfront = Analysis(folder, config, tag="FRONT", convtyp=0, verbose = verbose)
        FitRunnerback = Analysis(folder, config, tag="BACK", convtyp=1, verbose = verbose)
        if not(xmlfile ==""):
            FitRunnerfront.obs.xmlfile = xmlfile
            FitRunnerback.obs.xmlfile = xmlfile
        FitB = FitRunnerback.CreateLikeObject()
        FitF = FitRunnerfront.CreateLikeObject()
        import SummedLikelihood
        Fit = SummedLikelihood.SummedLikelihood()
        Fit.addComponent(FitB)
        Fit.addComponent(FitF)
        FitRunner = FitRunnerback
    else:
        convtype = config['analysis']['convtype']
        # Create one obs instance
        FitRunner = Analysis(folder, config, tag="", convtyp=convtype, verbose = verbose)
        if not(xmlfile ==""):
            FitRunner.obs.xmlfile = xmlfile
        Fit = FitRunner.CreateLikeObject()
    return FitRunner,Fit

def run(infile):
    """Run an entire Fermi analysis (spectrum) by reading a config file"""
    config = get_config(infile)
    folder = config['out']
    os.system('mkdir -p ' + folder)

    FitRunner,Fit = GenAnalysisObjects(config)
    # create all the fit files and run gtlike
    FitRunner.PerformFit(Fit)

    Result = FitRunner.GetAndPrintResults(Fit)#Get and dump the target specific results
    if config['verbose'] == 'yes' :
        utils.GetFluxes(Fit,FitRunner.obs.Emin,FitRunner.obs.Emax) #print the flux of all the sources

    utils.DumpResult(Result, config)

    #plot the SED and model map if possible and asked
    if config['Spectrum']['ResultPlots'] == 'yes' :
        from enrico.constants import SpectrumPath
    	os.system("mkdir -p "+config['out'] + '/'+SpectrumPath+'/')
        if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
            FitRunner.ComputeSED(Fit)
        outXml = utils._dump_xml(config)
        if config['Spectrum']['SummedLike'] != 'yes': # the possiblity of making the model map is checked inside the function
            FitRunner.ModelMap(outXml)

    #  Make energy bins by running a *new* analysis
    Nbin = config['Ebin']['NumEnergyBins']
    energybin.RunEbin(folder,Nbin,Fit,FitRunner)

# @todo: Should this be a command line utility in bin?
if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)
    run(infile)
