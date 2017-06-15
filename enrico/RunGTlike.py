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
        configs  = [None]*4
        Fits     = [None]*4
        Analyses = [None]*4
        if isKey(config['ComponentAnalysis'],'FrontBack') == 'yes':
            from enrico.data import fermievtypes
            mes.info("Breaking the analysis in Front/Back events")
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            oldxml = config['file']['xml']
            for k,TYPE in enumerate(["FRONT", "BACK"]):
                configs[k] = ConfigObj(config)
                configs[k]['event']['evtype'] = fermievtypes[TYPE]
                try:
                    Analyses[k] = Analysis(folder, configs[k], \
                        configgeneric=config,\
                        tag=TYPE, verbose = verbose)
                    if not(xmlfile ==""): Analyses[k].obs.xmlfile = xmlfile
                    Fits[k] = Analyses[k].CreateLikeObject()
                    Fit.addComponent(Fits[k])
                except RuntimeError,e:
                    if 'RuntimeError: gtltcube execution failed' in str(e):
                        mes.warning("Event type %s is empty! Error is %s" %(TYPE,str(e)))
            FitRunner = Analyses[0]

        elif isKey(config['ComponentAnalysis'],'PSF') == 'yes':
            from enrico.data import fermievtypes
            mes.info("Breaking the analysis in PSF 0,1,2,3.")
            # Clone the configs
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            for k,TYPE in enumerate(["PSF0", "PSF1", "PSF2", "PSF3"]):
                configs[k] = ConfigObj(config)
                configs[k]['event']['evtype'] = fermievtypes[TYPE]
                try:
                    Analyses[k] = Analysis(folder, configs[k], \
                        configgeneric=config,\
                        tag=TYPE, verbose = verbose)
                    if not(xmlfile ==""): Analyses[k].obs.xmlfile = xmlfile
                    Fits[k] = Analyses[k].CreateLikeObject()
                    Fit.addComponent(Fits[k])
                except RuntimeError,e:
                    if 'RuntimeError: gtltcube execution failed' in str(e):
                        mes.warning("Event type %s is empty! Error is %s" %(TYPE,str(e)))
            FitRunner = Analyses[0]

        elif isKey(config['ComponentAnalysis'],'EDISP') == 'yes':
            from enrico.data import fermievtypes
            mes.info("Breaking the analysis in EDISP 0,1,2,3.")
            # Clone the configs
            # Set Summed Likelihood to True
            config['Spectrum']['SummedLike'] = 'yes'
            for k,TYPE in enumerate(["EDISP0", "EDISP1", "EDISP2", "EDISP3"]):
                configs[k] = ConfigObj(config)
                configs[k]['event']['evtype'] = fermievtypes[TYPE]
                try:
                    Analyses[k] = Analysis(folder, configs[k], \
                        configgeneric=config,\
                        tag=TYPE, verbose = verbose)
                    if not(xmlfile ==""): Analyses[k].obs.xmlfile = xmlfile
                    Fits[k] = Analyses[k].CreateLikeObject()
                    Fit.addComponent(Fits[k])
                except RuntimeError,e:
                    if 'RuntimeError: gtltcube execution failed' in str(e):
                        mes.warning("Event type %s is empty! Error is %s" %(TYPE,str(e)))
            FitRunner = Analyses[0]

        elif isKey(config['ComponentAnalysis'],'EUnBinned')>=0:
            mes.info("Breaking the analysis in Binned (low energy) and Unbinned (high energies)")
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
            oldxml = config['file']['xml']
            for k,TYPE in enumerate(analysestorun):
                configs[k] = ConfigObj(config)
                # Tune parameters
                if TYPE is "lowE":
                    configs[k]['energy']['emin'] = emintotal
                    configs[k]['energy']['emax'] = min(config['energy']['emax'],EUnBinned)
                    configs[k]['analysis']['likelihood'] = "binned"
                    configs[k]['analysis']['ComputeDiffrsp'] = "no"
                elif TYPE is "highE":
                    configs[k]['energy']['emin'] = max(config['energy']['emin'],EUnBinned)
                    configs[k]['energy']['emax'] = emaxtotal
                    configs[k]['analysis']['likelihood'] = "unbinned"
                    configs[k]['analysis']['ComputeDiffrsp'] = "yes"
                try:
                    Analyses[k] = Analysis(folder, configs[k], \
                        configgeneric=config,\
                        tag=TYPE,\
                        verbose=verbose)
                    if not(xmlfile ==""): Analyses[k].obs.xmlfile = xmlfile
                    Fits[k] = Analyses[k].CreateLikeObject()
                    Fit.addComponent(Fits[k])
                except RuntimeError,e:
                    raise
                    if 'RuntimeError: gtltcube execution failed' in str(e):
                        mes.warning("Event type %s is empty! Error is %s" %(TYPE,str(e)))
            FitRunner = Analyses[0]
            FitRunner.obs.Emin = emintotal
            FitRunner.obs.Emax = emaxtotal

    try:
        FitRunner.config
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

    #Get and dump the target specific results
    Result = FitRunner.GetAndPrintResults(Fit)
    utils.DumpResult(Result, config)

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



    if config['Spectrum']['ResultPlots'] == 'yes' :
        outXml = utils._dump_xml(config)
        if config['Spectrum']['SummedLike'] != 'yes':
            # the possiblity of making the model map is checked inside the function
            FitRunner.ModelMap(outXml)



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
