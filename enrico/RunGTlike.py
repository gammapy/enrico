#!/usr/bin/env python
import os,os.path,math
from enrico import utils
from enrico.gtfunction import Observation
from enrico.fitmaker import FitMaker
import Loggin
import SummedLikelihood
from enrico.extern.configobj import ConfigObj
from utils import hasKey, isKey, typeirfs

def Analysis(folder, config, configgeneric=None, tag="", convtyp='-1', verbose = 1):

    mes = Loggin.Message()
    """ run an analysis"""
    Obs = Observation(folder, config, tag=tag)
    if verbose:
        utils._log('SUMMARY: ' + tag)
        Obs.printSum()

    FitRunner = FitMaker(Obs, config)##Class
    if config['Spectrum']['FitsGeneration'] == 'yes':
        FitRunner.FirstSelection(configgeneric) #Generates fits files for the coarse selection
        FitRunner.GenerateFits() #Generates fits files for the rest of the products
    return FitRunner

def GenAnalysisObjects(config, verbose = 1, xmlfile =""):

    mes = Loggin.Message()
    #check is the summed likelihood method should be used and get the
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    folder = config['out']

    Fit = SummedLikelihood.SummedLikelihood()

    EUnBinned = config['ComponentAnalysis']['EUnBinned']
    emintotal = float(config['energy']['emin'])
    emaxtotal = float(config['energy']['emax'])

    evtnum = [config["event"]["evtype"]] #for std analysis
    evtold = evtnum[0] #for std analysis
 
    # Create one obs instance for each component
    if isKey(config['ComponentAnalysis'],'FrontBack') == 'yes':
        evtnum = [1, 2]
    if isKey(config['ComponentAnalysis'],'PSF') == 'yes':
        evtnum = [4,8,16,32]
    if isKey(config['ComponentAnalysis'],'EDISP') == 'yes':
        evtnum = [64,128,256,521]
    oldxml = config['file']['xml']
    for k,evt in enumerate(evtnum):
        config['event']['evtype'] = evt
        config["file"]["xml"] = oldxml.replace(".xml","_"+typeirfs[evt]+".xml").replace("_.xml",".xml")

        if EUnBinned>emintotal and EUnBinned<emaxtotal:
            mes.info("Breaking the analysis in Binned (low energy) and Unbinned (high energies)")
            analysestorun = ["lowE","highE"]

            for k,TYPE in enumerate(analysestorun):
                tag = TYPE
                if typeirfs[evt] != "" : tag += "_"+typeirfs[evt]# handle name of fits file

                # Tune parameters
                if TYPE is "lowE":
                    config['energy']['emin'] = emintotal
                    config['energy']['emax'] = min(config['energy']['emax'],EUnBinned)
                    config['analysis']['likelihood'] = "binned"
                    config['analysis']['ComputeDiffrsp'] = "no"
                elif TYPE is "highE":
                    config['energy']['emin'] = max(config['energy']['emin'],EUnBinned)
                    config['energy']['emax'] = emaxtotal
                    config['analysis']['likelihood'] = "unbinned"
                    config['analysis']['ComputeDiffrsp'] = "yes"

                Analyse = Analysis(folder, config, \
                    configgeneric=config,\
                    tag=TYPE,\
                    verbose=verbose)


                Fit_component = Analyse.CreateLikeObject()
                Fit.addComponent(Fit_component)
            FitRunner = Analyse
            FitRunner.obs.Emin = emintotal
            FitRunner.obs.Emax = emaxtotal

        else:
            Analyse = Analysis(folder, config, \
                configgeneric=config,\
                tag=typeirfs[evt], verbose = verbose)

            # if not(xmlfile ==""): Analyse.obs.xmlfile = xmlfile
            Fit_component = Analyse.CreateLikeObject()
            Fit.addComponent(Fit_component)
    FitRunner = Analyse

    config["event"]["evtype"] = evtold
    FitRunner.config = config

    return FitRunner,Fit

def run(infile):
    from enrico import utils
    from enrico import energybin
    from enrico.config import get_config
    from enrico import Loggin
    mes = Loggin.Message()

    """Run an entire Fermi analysis (spectrum) by reading a config file"""
    config = get_config(infile)
    folder = config['out']
    utils.create_dir(folder)

    FitRunner,Fit = GenAnalysisObjects(config)
    # create all the fit files and run gtlike
    FitRunner.PerformFit(Fit)
    sedresult = None

    #plot the SED and model map if possible and asked
    if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
        if config['Spectrum']['ResultPlots'] == 'yes':
            from enrico.constants import SpectrumPath
            utils.create_dir("%s/%s/" %(config['out'],SpectrumPath))
            sedresult = FitRunner.ComputeSED(Fit,dump=True)
        else:
            sedresult = FitRunner.ComputeSED(Fit,dump=False)
        
        if (config['energy']['decorrelation_energy'] == 'yes'):
            #Update the energy scale to decorrelation energy
            mes.info('Setting the decorrelation energy as new Scale for the spectral parameters')
            spectrum = Fit[FitRunner.obs.srcname].funcs['Spectrum']
            modeltype = spectrum.genericName()
            genericName = Fit.model.srcs[FitRunner.obs.srcname].spectrum().genericName()

            varscale = None
            if genericName=="PowerLaw2":
                varscale = None
            elif genericName in ["PowerLaw", "PLSuperExpCutoff", "EblAtten::PLSuperExpCutoff"]:
                varscale = "Scale"
            elif genericName in ["LogParabola","EblAtten::LogParabola", \
                                 "BrokenPowerLaw", "EblAtten::BrokenPowerLaw"]:
                varscale = "Eb"

            if varscale is not None:
                spectrum.getParam(varscale).setValue(sedresult.decE)
                FitRunner.PerformFit(Fit)
            
    #Get and dump the target specific results
    Result = FitRunner.GetAndPrintResults(Fit)
    utils.DumpResult(Result, config)

    if config['Spectrum']['ResultPlots'] == 'yes' :
        outXml = utils._dump_xml(config)
        # the possibility of making the model map is checked inside the function
        FitRunner.ModelMap(outXml)

    #  Make energy bins by running a *new* analysis
    Nbin = config['Ebin']['NumEnergyBins']

    energybin.RunEbin(folder,Nbin,Fit,FitRunner,sedresult)
    
    del(sedresult)
    del(Result)
    del(FitRunner)

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
