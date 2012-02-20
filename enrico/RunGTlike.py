import os
import SummedLikelihood
import Utility
from config import get_config


def run(infile):
    """@todo: document me"""
    config = get_config(infile)
    folder = config['out']
    os.system('mkdir -p ' + folder)
    SummedLike = config['Spectrum']['SummedLike']
    if SummedLike == 'yes':
        # Create two obs instances
        runfitfront, _ = Utility.Analysis(folder, config, tag="FRONT", convtyp=0)
        runfitback, _ = Utility.Analysis(folder, config, tag="BACK", convtyp=1)
        FitB = runfitback.CreateFit()
        FitF = runfitfront.CreateFit()
        Fit = SummedLikelihood.SummedLikelihood()
        Fit.addComponent(FitB)
        Fit.addComponent(FitF)
        runfit = runfitback
    else:
        convtype = config['analysis']['convtype']
        # Create one obs instance
        runfit, _ = Utility.Analysis(folder, config, tag="", convtyp=convtype)
        Fit = runfit.CreateFit()
    Result = runfit.PerformFit(Fit)
    Utility.DumpResult(Result, config)

    if config['Spectrum']['ResultPlots'] == 'yes':
        runfit.PlotSED(Fit)
        outXml = (folder + "/" + runfit.obs.srcname + "_" +
                  config['file']['tag'] + "_out.xml")
        if SummedLike == 'yes':
            runfitback.ModelMap(outXml)
            runfitfront.ModelMap(outXml)
        else:
            runfit.ModelMap(outXml)

    if int(config['Ebin']['NumEnergyBins']) > 0:
        configfiles = Utility.PrepareEbin(Fit, runfit)
        for conf in configfiles:
            os.system('enrico_submit ' + conf)

# @todo: Should this be a command line utility in bin?
if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)
    run(infile)
