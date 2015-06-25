#!/usr/bin/env python
import os
from enrico import utils
from enrico import energybin
from enrico.config import get_config
from enrico.gtfunction import Observation
from enrico.fitmaker import FitMaker


def Analysis(folder, config, tag="" ,verbose = 1):
    """ run an analysis"""
    Obs = Observation(folder, config, tag=tag)
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
        sys.exit("not yet working")
        FitRunnerfront = Analysis(folder, config, tag="FRONT", verbose = verbose)
        FitRunnerback = Analysis(folder, config, tag="BACK", verbose = verbose)
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
        # Create one obs instance
        FitRunner = Analysis(folder, config, tag="", verbose = verbose)
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
    del FitRunner

# @todo: Should this be a command line utility in bin?
if __name__ == '__main__':
    import sys
    from enrico import Loggin
    mes = Loggin.Message()
    try:
        infile = sys.argv[1]
    except:
        print('Usage: '+sys.argv[0]+' <config file name>')
        mes.error('Config file not found.')

    run(infile)
