#!/usr/bin/env python
import os,glob,os.path,math
from enrico import utils
from enrico.gtfunction import Observation
from enrico.fitmaker import FitMaker
from enrico.plotting import plot_sed_fromconfig
from enrico import Loggin
import SummedLikelihood
from enrico.xml_model import XmlMaker
from enrico.extern.configobj import ConfigObj
from enrico.utils import hasKey, isKey, typeirfs

# Called per component
def Analysis(folder, config, configgeneric=None, tag="", convtyp='-1', verbose = 1):
    mes = Loggin.Message()
    """ run an analysis"""
    # If there are no xml files, create it and print a warning <--- This should be here?
    if len(glob.glob(config['file']['xml']))==0:
        mes.warning("Xml not found, creating one for the given config %s" %config['file']['xml'])
        XmlMaker(config)
    
    Obs = Observation(folder, config, tag=tag)
    if verbose:
        utils._log('SUMMARY: ' + tag)
        Obs.printSum()

    FitRunner = FitMaker(Obs, config)##Class
    if config['Spectrum']['FitsGeneration'] == 'yes':
        FitRunner.FirstSelection(configgeneric) #Generates fits files for the coarse selection
        FitRunner.GenerateFits() #Generates fits files for the rest of the products
    return FitRunner

# Called once
def GenAnalysisObjects(config, verbose = 1, xmlfile =""):
    # Array containing the list of analysis objects (needed to produce the individual residual maps)
    ListOfAnalysisObjects = []
    mes = Loggin.Message()
    #check is the summed likelihood method should be used and get the
    #Analysis objects (observation and (Un)BinnedAnalysis objects)
    folder = config['out']

    # If there are no xml files, create it and print a warning <--- Not sure if this is needed here.
    Fit = SummedLikelihood.SummedLikelihood()
    EUnBinned = config['ComponentAnalysis']['EUnBinned']
    emintotal = float(config['energy']['emin'])
    emaxtotal = float(config['energy']['emax'])

    evtnum = [config["event"]["evtype"]] #for std analysis
    evtold = evtnum[0] #for std analysis
 
    # Create one obs instance for each component. 
    # The first 3 can be combined with splitting in energy. The 4th tries to mimick 4FGL.
    if isKey(config['ComponentAnalysis'],'FrontBack') == 'yes':
        evtnum = [1, 2]
        config['analysis']['likelihood'] = "binned"
    elif isKey(config['ComponentAnalysis'],'PSF') == 'yes':
        evtnum = [4,8,16,32]
        config['analysis']['likelihood'] = "binned"
    elif isKey(config['ComponentAnalysis'],'EDISP') == 'yes':
        evtnum = [64,128,256,521]
        config['analysis']['likelihood'] = "binned"
    elif isKey(config['ComponentAnalysis'],'FGL4') == 'yes':
        # Special case of the PSF component analysis, 
        # where up to 15 components (energy+PSF) are created following 
        # 4FGL prescription.
        from enrico.catalogComponents import evtnum, energybins, nbinsbins, zmaxbins, ringwidths, pixelsizes
        config['analysis']['likelihood'] = "binned"
        oldxml = config['file']['xml']

        bin_i = 0
        roi = 0
        # energybins is a dictionary containing an index and a pair of energies
        for ebin_i in energybins:
            # Restrict the analysis to the specified energy range in all cases.
            if emintotal>=energybins[ebin_i][1]:
                continue
            if emaxtotal<=energybins[ebin_i][0]:
                continue
           
            if (roi==0): roi = 2.*ringwidths[ebin_i]+4.
            zmax    = zmaxbins[ebin_i]
            nbinsE  = nbinsbins[ebin_i]
            energybin = energybins[ebin_i]
            
            for k,evt in enumerate(evtnum):
                pixel_size = pixelsizes[ebin_i][k]
                if pixel_size<0: continue
                tag     = "{0}_En{1}".format(typeirfs[evt],ebin_i)
                # Approximation, in the 4FGL the core radius changes from src to src!
                mes.info("Breaking the analysis in bins ~ 4FGL")
                config['event']['evtype'] = evt
                config["file"]["xml"] = oldxml.replace(".xml","_")+typeirfs[evt]+"_"+\
                                        "En{0}.xml".format(ebin_i)
                config["energy"]["emin"] = max(emintotal,energybin[0])
                config["energy"]["emax"] = min(emaxtotal,energybin[1])
                config["analysis"]["likelihood"] = "binned"
                config["analysis"]["ComputeDiffrsp"] = "no"
                config["analysis"]["enumbins_per_decade"] = \
                    int(1.*nbinsE/math.log10(energybin[1]/energybin[0])+0.5)
                config["space"]["rad"] = roi
                config["analysis"]["zmax"] = zmax
                
                Analyse = Analysis(folder, config, \
                    configgeneric=config,\
                    tag=tag, verbose=verbose)
                ListOfAnalysisObjects.append(Analyse)
                
                if not(xmlfile ==""): Analyse.obs.xmlfile = xmlfile
                mes.info('Creating Likelihood object for component.')
                Fit_component = Analyse.CreateLikeObject()
                mes.info('Adding component to the summed likelihood.')
                Fit.addComponent(Fit_component)
      
		# Exports the weights that internally gtlike uses. These are hyper-cubes (component, energy, position). 
		# This is not fully correct. There should be an unbinned way of doing these steps, so it should be possible to generate weights
		if config['analysis']['likelihood'] == 'binned':
			for Analyse in ListOfAnalysisObjects:
				Analyse.obs.GtEffBkg()
			Analyse.obs.GtAlphaBkg()	
			for Analyse in ListOfAnalysisObjects:
				Analyse.obs.GtWtsMap()
 
        FitRunner = Analyse
        FitRunner.obs.Emin = emintotal
        FitRunner.obs.Emax = emaxtotal
        config["energy"]["emin"] = emintotal
        config["energy"]["emax"] = emaxtotal
        config["event"]["evtype"] = evtold
        FitRunner.config = config

        return FitRunner,Fit,ListOfAnalysisObjects

    # Standard (non-4FGL) analysis components
    oldxml = config['file']['xml']
    for k,evt in enumerate(evtnum):
        config['event']['evtype'] = evt
        if typeirfs[evt] != "" and typeirfs[evt]!="FRONTBACK":
            config["file"]["xml"] = oldxml.replace(".xml","_"+typeirfs[evt]+".xml")

        if EUnBinned>emintotal and EUnBinned<emaxtotal:
            mes.info("Breaking the analysis in Binned (low energy) and Unbinned (high energies)")
            analysestorun = ["lowE","highE"]

            for j,TYPE in enumerate(analysestorun):
                tag = TYPE
                if typeirfs[evt] != "" : tag += "_"+typeirfs[evt]# handle name of fits file
                config["file"]["xml"] = oldxml.replace(".xml","_"+tag+".xml")

                # Tune parameters
                if TYPE == "lowE":
                    config['energy']['emin'] = emintotal
                    config['energy']['emax'] = min(config['energy']['emax'],EUnBinned)
                    config['analysis']['likelihood'] = "binned"
                    config['analysis']['ComputeDiffrsp'] = "no"
                elif TYPE == "highE":
                    config['energy']['emin'] = max(config['energy']['emin'],EUnBinned)
                    config['energy']['emax'] = emaxtotal
                    config['analysis']['likelihood'] = "unbinned"
                    config['analysis']['ComputeDiffrsp'] = "yes"

                Analyse = Analysis(folder, config, \
                    configgeneric=config,\
                    tag=tag,\
                    verbose=verbose)
                ListOfAnalysisObjects.append(Analyse)

                mes.info('Creating Likelihood object for component.')
                Fit_component = Analyse.CreateLikeObject()
                mes.info('Adding component to the summed likelihood.')
                Fit.addComponent(Fit_component)
            FitRunner = Analyse
            FitRunner.obs.Emin = emintotal
            FitRunner.obs.Emax = emaxtotal
            config["energy"]["emin"] = emintotal
            config["energy"]["emax"] = emaxtotal

        else:
            Analyse = Analysis(folder, config, \
                configgeneric=config,\
                tag=typeirfs[evt], verbose = verbose)
            
            ListOfAnalysisObjects.append(Analyse)
            if not(xmlfile ==""): Analyse.obs.xmlfile = xmlfile
            mes.info('Creating Likelihood object for component.')
            Fit_component = Analyse.CreateLikeObject()
            mes.info('Adding component to the summed likelihood.')
            Fit.addComponent(Fit_component)
   
    FitRunner = Analyse
    config["event"]["evtype"] = evtold
    FitRunner.config = config

    return FitRunner,Fit,ListOfAnalysisObjects

def run(infile):
    from enrico import utils
    from enrico import energybin
    from enrico.config import get_config
    from enrico import Loggin
    mes = Loggin.Message()

    """Run an entire Fermi analysis (spectrum) by reading a config file"""
    config = get_config(infile)
    folder = config['out']
    utils.mkdir_p(folder)

    FitRunner,Fit,ListOfAnalysisObjects = GenAnalysisObjects(config)
    # create all the fit files and run gtlike
    FitRunner.PerformFit(Fit)
    sedresult = None

    #plot the SED and model map if possible and asked
    if float(config['UpperLimit']['TSlimit']) < Fit.Ts(config['target']['name']):
        if config['Spectrum']['ResultPlots'] == 'yes':
            from enrico.constants import SpectrumPath
            utils.mkdir_p("%s/%s/" %(config['out'],SpectrumPath))
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
                try:
                    spectrum.getParam(varscale).setBounds(20,3e6)
                    spectrum.getParam(varscale).setValue(sedresult.decE)
                    FitRunner.PerformFit(Fit)
                except RuntimeError:
                    mes.warning("Error occurred while setting decorrelation energy.") 
            
    #Get and dump the target specific results
    Result = FitRunner.GetAndPrintResults(Fit)
    utils.DumpResult(Result, config)
    
    #  Make energy bins by running a *new* analysis
    Nbin = config['Ebin']['NumEnergyBins']
    
    if (FitRunner.config['file']['parent_config']==""):
        FitRunner.config['file']['parent_config'] = infile
    
    if config['Spectrum']['ResultParentPlots'] == "yes":
        print(config['file']['parent_config'])
        plot_sed_fromconfig(config['file']['parent_config'],ignore_missing_bins=True) 
    
    if config['Spectrum']['ResultPlots'] == 'yes' :
        outXml = utils._dump_xml(config)
        # the possibility of making the model map is checked inside the function
        for AnalysisComponent in ListOfAnalysisObjects:
            AnalysisComponent.obs.ModelMap(outXml) 
        if Nbin>0:
            FitRunner.config['Spectrum']['ResultParentPlots'] = "yes"
        plot_sed_fromconfig(infile,ignore_missing_bins=True)
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
        print(('Usage: '+sys.argv[0]+' <config file name>'))
        mes.error('Config file not found.')

    run(infile)
