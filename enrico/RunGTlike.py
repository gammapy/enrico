#!/usr/bin/env python

def Analysis(folder, config, configgeneric=None, tag="", convtyp='-1', verbose = 1):
    import os
    from enrico import utils
    from enrico.gtfunction import Observation
    from enrico.fitmaker import FitMaker
    import Loggin

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
    import os
    import os.path
    import math
    import SummedLikelihood
    from enrico.xml_model import XmlMaker
    from enrico.extern.configobj import ConfigObj
    from utils import hasKey, isKey
    import Loggin
    mes = Loggin.Message()
    #check is the summed likelihood method should be used and get the
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    SummedLike = config['Spectrum']['SummedLike']
    folder = config['out']
    
    # If there is no xml file, create it and print a warning
    if (not os.path.isfile(config['file']['xml'])):
        mes.warning("Xml not found, creating one for the given config %s" %config['file']['xml'])
        XmlMaker(config)
    
    Fit = SummedLikelihood.SummedLikelihood()

    if hasKey(config,'ComponentAnalysis') == True:
        # Create one obs instance for each component
        if isKey(config['ComponentAnalysis'],'FrontBack') == 'yes':
            mes.info("Breaking the analysis in Front/Back events")
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            FitRunnerfront = Analysis(folder, config, \
                configgeneric=config,\
                tag="FRONT", verbose = verbose)
            FitRunnerback = Analysis(folder, config, \
                configgeneric=config,\
                tag="BACK", verbose = verbose)
            if not(xmlfile ==""):
                FitRunnerfront.obs.xmlfile = xmlfile
                FitRunnerback.obs.xmlfile = xmlfile
            FitB = FitRunnerback.CreateLikeObject()
            FitF = FitRunnerfront.CreateLikeObject()
            Fit.addComponent(FitB)
            Fit.addComponent(FitF)
            if verbose: print(Fit.components)
            FitRunner = FitRunnerback
        elif isKey(config['ComponentAnalysis'],'PSF') == 'yes':
            mes.info("Breaking the analysis in PSF 0,1,2,3.")
            # Clone the configs
            config_psfs  = [None]*4
            config_xmls  = [None]*4
            FitPSFs      = [None]*4
            AnalysisPSFs = [None]*4
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            for k in xrange(4):
                config_psfs[k] = ConfigObj(config)
                # Tune parameters
                config_psfs[k]['event']['evtype'] = int(2**(k+2))
                oldxml = config_psfs[k]['file']['xml']
                AnalysisPSFs[k] = Analysis(folder, config_psfs[k], \
                    configgeneric=config,\
                    tag="PSF%d"%k,\
                    verbose = verbose)
                if not(xmlfile ==""):
                    AnalysisPSF[k].obs.xmlfile = xmlfile
                FitPSFs[k] = AnalysisPSFs[k].CreateLikeObject()
                Fit.addComponent(FitPSFs[k])
            if verbose: print(Fit.components)
            FitRunner = AnalysisPSFs[0]
        elif isKey(config['ComponentAnalysis'],'EDISP') == 'yes':
            mes.info("Breaking the analysis in EDISP 0,1,2,3.")
            # Clone the configs
            config_edisps  = [None]*4
            config_xmls    = [None]*4
            FitEDISPs      = [None]*4
            AnalysisEDISPs = [None]*4
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            for k in xrange(4):
                config_edisps[k] = ConfigObj(config)
                # Tune parameters
                config_edisps[k]['event']['evtype'] = int(2**(k+6))
                oldxml = config_edisps[k]['file']['xml']
                AnalysisEDISPs[k] = Analysis(folder, config_edisps[k], \
                    configgeneric=config,\
                    tag="EDISP%d"%k,\
                    verbose = verbose)
                if not(xmlfile ==""):
                    AnalysisEDISP[k].obs.xmlfile = xmlfile
                FitEDISPs[k] = AnalysisEDISPs[k].CreateLikeObject()
                Fit.addComponent(FitEDISPs[k])
            if verbose: print(Fit.components)
            FitRunner = AnalysisEDISPs[0]
        elif isKey(config['ComponentAnalysis'],'EUnBinned')>=0:
            mes.info("Breaking the analysis in Binned (low energy) and Unbinned (high energies)")
            # Clone the configs
            config_bin   = [None]*2
            config_xmls  = [None]*2
            FitBIN       = [None]*2
            AnalysisBIN  = [None]*2
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            EUnBinned = config['ComponentAnalysis']['EUnBinned']
            emintotal = float(config['energy']['emin'])
            emaxtotal = float(config['energy']['emax'])
            # Run the following analysis depending on the case
            # this is general, no matter if we are in the total fit or the Ebin #N fit. 

            if EUnBinned<=emintotal:    analysestorun = ["highE"]
            elif EUnBinned>=emaxtotal:  analysestorun = ["lowE"]
            else:                       analysestorun = ["lowE","highE"]
            
            print(analysestorun)

            oldxml = config['file']['xml']
            for k,name in enumerate(analysestorun):
                config_bin[k] = ConfigObj(config)
                # Tune parameters
                if name is "lowE":
                    config_bin[k]['energy']['emin'] = emintotal
                    config_bin[k]['energy']['emax'] = min(config['energy']['emax'],EUnBinned)
                    config_bin[k]['analysis']['likelihood'] = "binned"
                    config_bin[k]['analysis']['ComputeDiffrsp'] = "no"
                elif name is "highE":
                    config_bin[k]['energy']['emin'] = max(config['energy']['emin'],EUnBinned)
                    config_bin[k]['energy']['emax'] = emaxtotal
                    config_bin[k]['analysis']['likelihood'] = "unbinned"
                    config_bin[k]['analysis']['ComputeDiffrsp'] = "yes"
                AnalysisBIN[k] = Analysis(folder, config_bin[k], \
                    configgeneric=config,\
                    tag=name,\
                    verbose=verbose)
                if not(xmlfile ==""):
                    AnalysisBIN[k].obs.xmlfile = xmlfile
                FitBIN[k] = AnalysisBIN[k].CreateLikeObject()
                Fit.addComponent(FitBIN[k])
            if verbose: print(Fit.components)
            FitRunner = AnalysisBIN[0]
            FitRunner.obs.Emin = emintotal
            FitRunner.obs.Emax = emaxtotal
    
    try:
        FitRunner
    except NameError:
        # Create one obs instance
        FitRunner = Analysis(folder, config, tag="", verbose = verbose)
        Fit.addComponent(FitRunner.CreateLikeObject())
        if not(xmlfile ==""):
            FitRunner.obs.xmlfile = xmlfile
    
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

    #plot the SED and model map if possible and asked
    if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
        if config['Spectrum']['ResultPlots'] == 'yes':
            from enrico.constants import SpectrumPath
            utils.create_dir("%s/%s/" %(config['out'],SpectrumPath))
            sedresult = FitRunner.ComputeSED(Fit,dump=True)
        else:
            sedresult = FitRunner.ComputeSED(Fit,dump=False)
    
        # Update the energy scale to decorrelation energy
        mes.info('Setting the decorrelation energy as new Scale for the spectral parameters')
        spectrum = Fit[FitRunner.obs.srcname].funcs['Spectrum']
        modeltype = spectrum.genericName()
        if Fit.model.srcs[FitRunner.obs.srcname].spectrum().genericName()=="PowerLaw":
            varscale = "Scale"
        if Fit.model.srcs[FitRunner.obs.srcname].spectrum().genericName()=="PowerLaw2":
            varscale = None
        elif Fit.model.srcs[FitRunner.obs.srcname].spectrum().genericName()=="PLSuperExpCutoff":
            varscale = "Scale"
        elif Fit.model.srcs[FitRunner.obs.srcname].spectrum().genericName()=="LogParabola":
            varscale = "Eb"
        elif Fit.model.srcs[FitRunner.obs.srcname].spectrum().genericName()=="BrokenPowerLaw":
            varscale = "Eb"
        
        if varscale is not None:
            spectrum.getParam(varscale).setValue(sedresult.decE)
            FitRunner.PerformFit(Fit)

    if config['Spectrum']['ResultPlots'] == 'yes' :
        outXml = utils._dump_xml(config)
        if config['Spectrum']['SummedLike'] != 'yes':
            # the possiblity of making the model map is checked inside the function
            FitRunner.ModelMap(outXml)
    
    #Get and dump the target specific results
    Result = FitRunner.GetAndPrintResults(Fit)
    utils.DumpResult(Result, config)

    #  Make energy bins by running a *new* analysis
    Nbin = config['Ebin']['NumEnergyBins']
    energybin.RunEbin(folder,Nbin,Fit,FitRunner)

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
