#!/usr/bin/env python


def Analysis(folder, config, tag="", convtyp='-1', verbose = 1):
    import os
    from enrico import utils
    from enrico.gtfunction import Observation
    from enrico.fitmaker import FitMaker
    from enrico.xml_model import XmlMaker
    import Loggin

    mes = Loggin.Message()

    # If there is no xml file, create it and print a warning
    if (not os.path.isfile(config['file']['xml'])):
        mes.warning("Xml not found, creating one for the given config %s" %config['file']['xml'])
        XmlMaker(config)


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
    from enrico.extern.configobj import ConfigObj
    from utils import hasKey, isKey
    import Loggin
    mes = Loggin.Message()
    #check is the summed likelihood method should be used and get the
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    SummedLike = config['Spectrum']['SummedLike']
    folder = config['out']
    if SummedLike == 'yes':
        import sys
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
    elif hasKey(config,'ComponentAnalysis') == True:
        # Create one obs instance for each component
        if isKey(config['ComponentAnalysis'],'PSF') == 'yes':
            mes.info("Breaking the analysis in PSF 0,1,2,3.")
            # Clone the configs
            config_psfs  = [None]*4
            config_xmls  = [None]*4
            FitPSFs      = [None]*4
            AnalysisPSFs = [None]*4
            import SummedLikelihood
            Fit = SummedLikelihood.SummedLikelihood()
            for k in xrange(4):
                config_psfs[k] = ConfigObj(config)
                # Tune parameters
                config_psfs[k]['event']['evtype'] = int(2**(k+2))
                config_psfs[k]['file']['tag'] = str("PSF%d" %k)
                oldxml = config_psfs[k]['file']['xml']
                config_psfs[k]['file']['xml'] = oldxml.replace('model.xml','model_PSF%d.xml'%k)
                AnalysisPSFs[k] = Analysis(folder, config_psfs[k], \
                    tag=config_psfs[k]['file']['tag'], \
                    verbose = verbose)
                FitPSFs[k] = AnalysisPSFs[k].CreateLikeObject()
                Fit.addComponent(FitPSFs[k])
            FitRunner = AnalysisPSFs[0]
        elif hasKey(config['ComponentAnalysis'],'EUnBinned'):
            mes.info("Breaking the analysis in Binned (low energy) and Unbinned (high energies)")
            # Clone the configs
            config_bin   = [None]*2
            config_xmls  = [None]*2
            FitBIN       = [None]*2
            AnalysisBIN  = [None]*2
            import SummedLikelihood
            Fit = SummedLikelihood.SummedLikelihood()
            for k,name in enumerate(["highE","lowE"]):
                config_bin[k] = ConfigObj(config)
                # Tune parameters
                config_bin[k]['event']['evtype'] = int(2**(k+6))
                config_bin[k]['file']['tag'] = str("%s" %name)
                config_bin[k]['file']['xml'].replace('model','model_%s'%name)
                if name is "lowE":
                    config_bin[k]['analysis']['emax']       = config['ComponentAnalysis']['EUnBinned']
                    config_bin[k]['analysis']['likelihood'] = "binned"
                elif name is "highE":
                    config_bin[k]['analysis']['emin']           = config['ComponentAnalysis']['EUnBinned']
                    config_bin[k]['analysis']['likelihood']     = "unbinned"
                    config_bin[k]['analysis']['ComputeDiffrsp'] = "yes"
                oldxml = config_bin[k]['file']['xml']
                config_bin[k]['file']['xml'] = oldxml.replace('model.xml','model_%s.xml'%name)
                AnalysisBIN[k] = Analysis(folder, config_bin[k], \
                    tag=config_bin[k]['file']['tag'], \
                    verbose=verbose)
                FitBIN[k] = AnalysisBIN[k].CreateLikeObject()
                Fit.addComponent(FitBIN[k])
            FitRunner = AnalysisBIN[0]
        elif hasKey(config['ComponentAnalysis'],'EDISP'):
            mes.info("Breaking the analysis in EDISP 0,1,2,3.")
            # Clone the configs
            config_edisps  = [None]*4
            config_xmls    = [None]*4
            FitEDISPs      = [None]*4
            AnalysisEDISPs = [None]*4
            import SummedLikelihood
            Fit = SummedLikelihood.SummedLikelihood()
            for k in xrange(4):
                config_edisps[k] = ConfigObj(config)
                # Tune parameters
                config_edisps[k]['event']['evtype'] = int(2**(k+6))
                config_edisps[k]['file']['tag'] = str("EDISP%d" %k)
                config_edisps[k]['file']['xml'].replace('model','model_EDISP%d'%k)
                oldxml = config_edisps[k]['file']['xml']
                config_edisp[k]['file']['xml'] = oldxml.replace('model.xml','model_EDISP%d.xml'%k)
                AnalysisEDISPs[k] = Analysis(folder, config_edisps[k], \
                    tag=config_edisps[k]['file']['tag'], \
                    verbose = verbose)
                FitEDISPs[k] = AnalysisEDISPs[k].CreateLikeObject()
                Fit.addComponent(FitEDISPs[k])
            FitRunner = AnalysisEDISPs[0]
    else:
        # Create one obs instance
        FitRunner = Analysis(folder, config, tag="", verbose = verbose)
    if not(xmlfile ==""):
        FitRunner.obs.xmlfile = xmlfile
    Fit = FitRunner.CreateLikeObject()
    return FitRunner,Fit

def run(infile):
    from enrico import utils
    from enrico import energybin
    from enrico.config import get_config

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
    if config['Spectrum']['ResultPlots'] == 'yes' :
        from enrico.constants import SpectrumPath
        utils.create_dir("%s/%s/" %(config['out'],SpectrumPath))
        if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
            FitRunner.ComputeSED(Fit)
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
