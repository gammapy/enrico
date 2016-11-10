#!/usr/bin/env python

def Analysis(folder, config, tag="", convtyp='-1', verbose = 1):
    import os
    from enrico import utils
    from enrico.gtfunction import Observation
    from enrico.fitmaker import FitMaker
    from enrico.xml_model import XmlMaker
    import Loggin

    mes = Loggin.Message()
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
    import os
    import os.path
    import math
    import SummedLikelihood
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
                tag="FRONT", verbose = verbose)
            FitRunnerback = Analysis(folder, config, \
                tag="BACK", verbose = verbose)
            if not(xmlfile ==""):
                FitRunnerfront.obs.xmlfile = xmlfile
                FitRunnerback.obs.xmlfile = xmlfile
            FitB = FitRunnerback.CreateLikeObject()
            FitF = FitRunnerfront.CreateLikeObject()
            Fit.addComponent(FitB)
            Fit.addComponent(FitF)
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
                #config_psfs[k]['file']['tag'] = str("PSF%d" %k)
                oldxml = config_psfs[k]['file']['xml']
                #config_psfs[k]['file']['xml'] = oldxml.replace('model.xml','model_PSF%d.xml'%k)
                AnalysisPSFs[k] = Analysis(folder, config_psfs[k], \
                    tag="PSF%d"%k,\
                    verbose = verbose)
                if not(xmlfile ==""):
                    FitPSFs[k].obs.xmlfile = xmlfile
                FitPSFs[k] = AnalysisPSFs[k].CreateLikeObject()
                Fit.addComponent(FitPSFs[k])
            FitRunner = AnalysisPSFs[0]
        elif isKey(config['ComponentAnalysis'],'EUnBinned') is not -1.:
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
            
            for k,name in enumerate(analysestorun):
                config_bin[k] = ConfigObj(config)
                # Tune parameters
                config_bin[k]['file']['tag'] = str("%s" %name)
                #config_bin[k]['file']['xml'].replace('model','model_%s'%name)
                if name is "lowE":
                    config_bin[k]['energy']['emax'] = min(config_bin[k]['energy']['emax'],EUnBinned)
                    config_bin[k]['analysis']['likelihood'] = "binned"
                elif name is "highE":
                    config_bin[k]['energy']['emin'] = max(config_bin[k]['energy']['emin'],EUnBinned)
                    config_bin[k]['analysis']['likelihood']     = "unbinned"
                    config_bin[k]['analysis']['ComputeDiffrsp'] = "yes"
                oldxml = config_bin[k]['file']['xml']
                #config_bin[k]['file']['xml'] = oldxml.replace('model.xml','model_%s.xml'%name)
                AnalysisBIN[k] = Analysis(folder, config_bin[k], \
                    tag=name,\
                    verbose=verbose)
                if not(xmlfile ==""):
                    FitBINs[k].obs.xmlfile = xmlfile
                AnalysisBIN[k].obs.Emin = emintotal
                AnalysisBIN[k].obs.Emax = emaxtotal
                FitBIN[k] = AnalysisBIN[k].CreateLikeObject()
                Fit.addComponent(FitBIN[k])
            FitRunner = AnalysisBIN[0]
            FitRunner.config = config

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
                config_edisps[k]['file']['tag'] = str("EDISP%d" %k)
                #config_edisps[k]['file']['xml'].replace('model','model_EDISP%d'%k)
                oldxml = config_edisps[k]['file']['xml']
                #config_edisp[k]['file']['xml'] = oldxml.replace('model.xml','model_EDISP%d.xml'%k)
                AnalysisEDISPs[k] = Analysis(folder, config_edisps[k], \
                    tag="EDISP%d"%k,\
                    verbose = verbose)
                if not(xmlfile ==""):
                    FitEDISPs[k].obs.xmlfile = xmlfile
                FitEDISPs[k] = AnalysisEDISPs[k].CreateLikeObject()
                Fit.addComponent(FitEDISPs[k])
            FitRunner = AnalysisEDISPs[0]
    
    try:
        FitRunner
    except NameError:
        # Create one obs instance
        FitRunner = Analysis(folder, config, tag="", verbose = verbose)
        Fit.addComponent(FitRunner.CreateLikeObject())
        if not(xmlfile ==""):
            FitRunner.obs.xmlfile = xmlfile
    
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
        spectrum.getParam("Scale").setValue(sedresult.decE)
        FitRunner.PerformFit(Fit)

    if config['Spectrum']['ResultPlots'] == 'yes' :
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
