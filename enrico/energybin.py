import os
import numpy as np
import traceback
from enrico import environ
from enrico.constants import EbinPath
from enrico.submit import call
from enrico.config import get_config
from enrico import utils, Loggin

def ChangeModel(comp, E1, E2, name, Pref, Gamma):
    """Change the spectral model of a source called name
    to allow a fit between E1 and E2
    If the spectral model is PowerLaw, the prefactor is updated
    if not the model is change to PowerLaw.
    The index is frozen in all cases"""
    
    # if approximated Gamma is outside of bounds set it to limit.
    # Same for the prefix, do not allow crazy values (>1 or <1e-25, e.g. 0.)
    Gamma_min=-5
    Gamma_max=-0.501
    Gamma=max(min(Gamma_max,Gamma),Gamma_min)
    Pref =max(min(1,Pref),1e-25)

    Eav = utils.GetE0(E1, E2)

    spectrum = comp.logLike.getSource(name).getSrcFuncs()['Spectrum']
    spectrum.getParam('Prefactor').setScale(utils.fluxScale(Pref))
    spectrum.getParam('Prefactor').setValue(utils.fluxNorm(Pref))
    spectrum.getParam('Prefactor').setBounds(1e-5,1e5)
    spectrum.getParam('Index').setValue(Gamma)
    spectrum.getParam('Index').setBounds(Gamma_min,Gamma_max)
    spectrum.getParam('Index').setFree(False)
    spectrum.getParam('Scale').setValue(Eav)
    spectrum.getParam('Scale').setBounds(20,3e6)

    return comp

def PrepareEbin(Fit, FitRunner,sedresult=None):
    """ Prepare the computation of spectral point in energy bins by
    i) removing the weak sources (TS<1) # not true
    ii) updating the config file (option and energy)
    and save it in a new ascii file
    iii) changing the spectral model and saving it in a new xml file.
    A list of the ascii files is returned"""
        
    mes = Loggin.Message()

    NEbin = int(FitRunner.config['Ebin']['NumEnergyBins'])

    config = FitRunner.config

    config['verbose'] ='no' #Be quiet

    #Replace the evt file with the fits file produced before
    #in order to speed up the production of the fits files
    config['file']['event'] = FitRunner.obs.eventcoarse
    #update the config to allow the fit in energy bins
    config['UpperLimit']['envelope'] = 'no'
    config['Ebin']['NumEnergyBins'] = '0'#no new bin in energy!
    config['target']['redshift']    = '0'#Disable EBL correction
    config['out'] = FitRunner.config['out'] + '/'+EbinPath + str(NEbin)
    config['Spectrum']['ResultPlots'] = 'no' #no SED plot/modelmap
    config['Spectrum']['WriteCovMatrix'] = 'no' #no cov matrix
    #copy the chose of the user for the enery bin computing
    config['Spectrum']['FitsGeneration'] = config['Ebin']['FitsGeneration']
    config['UpperLimit']['TSlimit'] = config['Ebin']['TSEnergyBins']
    tag = FitRunner.config['file']['tag']
    Emax = float(FitRunner.config['energy']['emax'])
    Emin = float(FitRunner.config['energy']['emin'])
    lEmax = np.log10(Emax)
    lEmin = np.log10(Emin)
    utils._log("Preparing submission of fit into energy bins")
    print("Emin = {0} MeV".format(Emin),
          "Emax = {0} MeV".format(Emax),
          "Nbins = {0}".format(NEbin))

    ener = utils.string_to_list(config['Ebin']['DistEbins'])

    if ener is None:
        if (config['ComponentAnalysis']['FGL4'] == 'yes' or config['Ebin']['DistEbins']=='FGL4'):
            ener  = np.asarray([50,1e2,3e2,1e3,3e3,1e4,3e4,3e5])
            NEbin = len(ener)-1
        elif config['Ebin']['DistEbins'] in ['TS','mix'] and sedresult!=None:
            # Make the bins equispaced in sum(SED/SEDerr) - using the butterfly
            ipo = 0
            iTS = sedresult.SED/sedresult.Err
            TScumula = 0
            TSperbin = 1.*sum(iTS)/NEbin
            ener = [10**lEmin]
            while ipo<len(sedresult.E)-1:
                TScumula += iTS[ipo]
                if TScumula/TSperbin > 1:
                    ener.append(sedresult.E[ipo])
                    TScumula -= TSperbin
                ipo += 1
            ener.append(10**lEmax)
            ener = np.array(ener)
            # intermediate approach (between both TS-spaced and logE spaced)
            if config['Ebin']['DistEbins'] == 'mix':
                ener = 0.5*(ener + np.logspace(lEmin, lEmax, NEbin + 1))
        else:
            # Make the bins equispaced in logE (standard)
            ener = np.logspace(lEmin, lEmax, NEbin + 1)
    
    # 1. Remove bins that are out of the range covered by the data
    # 2. Limit the bin extend to the range covered by the data. 
    # Get elements strictly above threshold +1 element to the left for the left side
    # Get elements strictly below limit +1 element to the right side.
    # example. [1,2,3,4,5] -> if Emin=3.4, Emax=3.9 we want to keep [3.4,3.9].
    ener = np.asarray(ener)
    print(("Energy bins (before energy cuts): {0}".format(str(ener))))
    if len(ener)==0: 
        print("** Warning: energy bin array is empty")
        return(None)
    available_left = ener>Emin              # In the example FFFTT -> [4,5]
    for k,use in enumerate(available_left[:-1]):
        if not use and available_left[k+1]: 
            available_left[k] = True        # In the example FFTTT -> [3,5]
    available_right = ener<Emax             # In the example TTTFF -> [1,3]
    for k,use in enumerate(available_right[1:]):
        if not use and available_right[k]:
            available_right[k+1] = True     # In the example TTTTF -> [1,4]
    available = available_left*available_right
    ener  = ener[available]                 # In the example FFTTF -> [3,4]
    # Limit the range to the real energies that are covered by our data
    # If the energy bins are well placed this should not do anything.
    ener[0]  = np.max([Emin,ener[0]])
    ener[-1] = np.min([Emax,ener[-1]])
    NEbin = len(ener)-1
    print(("Energy bins (after energy cuts): {0}".format(str(ener))))
    if len(ener)==0: 
        print("** Warning: energy bin array is empty")
        return(None)

    utils.mkdir_p(config['out'])
    paramsfile = []

    srcname = FitRunner.config['target']['name']
    try:
        TSsrc = Fit.Ts(srcname)
    except RuntimeError:
        TSsrc = 0

    if config['UpperLimit']['TSlimit']>TSsrc:
        utils._log('Re-optimize', False)
        print("An upper limit has been computed. The fit need to be re-optimized")
        Fit.optimize(0)
    
    Pref = utils.ApproxPref(Fit, ener, srcname)
    Gamma = utils.ApproxGamma(Fit, ener, srcname)

    Model_type = Fit.model.srcs[srcname].spectrum().genericName()
    # if the model is not PowerLaw : change the model
    if not(Model_type == 'PowerLaw'):
        for comp in Fit.components:
            comp.logLike.getSource(srcname).setSpectrum("PowerLaw") #Change model
        config['target']['spectrum'] = "PowerLaw"

    xmltag_list = [""] #handle summed like analysis
    if config['ComponentAnalysis']['FrontBack'] == "yes":
        xmltag_list = ["_FRONT","_BACK"]
        mes.info("Splitting Front/Back events")
    elif config['ComponentAnalysis']['PSF'] == "yes":
        xmltag_list = ["_PSF0","_PSF1","_PSF2","_PSF3"]
        mes.info("Splitting PSF events")
    elif config['ComponentAnalysis']['EDISP'] == "yes":
        xmltag_list = ["_EDISP0","_EDISP1","_EDISP2","_EDISP3"]
        mes.info("Splitting EDISP events")
    elif config['ComponentAnalysis']['FGL4'] == "yes":
        from enrico.catalogComponents import evtnum, energybins, pixelsizes
        xmltag_list = []
        for ebin_i in energybins:
            for k,evt in enumerate(evtnum):
                #if pixelsizes[ebin_i][k] > 0:
                try:
                    xmltag_list.append("_{0}_En{1}".format(utils.typeirfs[k],ebin_i))
                except KeyError:
                    continue

    for ibin in range(NEbin):#Loop over the energy bins
        E = utils.GetE0(ener[ibin + 1],ener[ibin])
        mes.info("Submitting # "+str(ibin)+" at energy "+str(E))
        #Update the model for the bin
        for comp,xmltag in zip(Fit.components, xmltag_list):
            NewFitObject = ChangeModel(comp, ener[ibin], ener[ibin + 1], srcname, Pref[ibin] ,Gamma[ibin])
            Xmlname = (config['out'] + "/" + srcname + "_" + str(ibin) +xmltag+ ".xml")

            NewFitObject.writeXml(Xmlname)# dump the corresponding xml file
            config['file']['xml'] = Xmlname.replace(xmltag,"")
        #update the energy bounds
        config['energy']['emin'] = str(ener[ibin])
        config['energy']['emax'] = str(ener[ibin + 1])
        config['energy']['decorrelation_energy'] = "no"
        # Change the spectral index to follow the Estimated Gamma 
        # if approximated Gamma is outside of bounds set it to limit
        Gamma_min=-5
        Gamma_max=-0.501
        Gamma_bin=-max(min(Gamma_max,Gamma[ibin]),Gamma_min)
        config['Spectrum']['FrozenSpectralIndex'] = Gamma_bin
        config['UpperLimit']['SpectralIndex'] = Gamma_bin
        
        # Force ebin jobs to be sequential ?
        if config['Ebin']['Submit'] == 'yes':
            config['Submit'] = 'yes'
        elif config['Ebin']['Submit'] == 'no':
            config['Submit'] = 'no'
        
        config['file']['tag'] = tag + '_Ebin' + str(NEbin) + '_' + str(ibin)
        filename =  config['target']['name'] + "_" + str(ibin) + ".conf"
        paramsfile.append(filename)
        config.filename = config['out'] + '/' +filename
        config.write()
        #config.write(open(config['out'] + '/' +filename, 'w')) #save the config file in a ascii file

    return paramsfile


def RunEbin(folder,Nbin,Fit,FitRunner,sedresult=None):
    if int(Nbin) > 0:
        configfiles = PrepareEbin(Fit, FitRunner,sedresult)
        ind = 0
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')
        for conf in configfiles:
             pathconf = folder + "/"+ EbinPath + str(Nbin) +"/" + conf
             Newconfig = get_config(pathconf)
             cmd = enricodir+"/enrico/RunGTlike.py "+pathconf
             if Newconfig['Submit'] == 'no' : #run directly
                 try:
                     os.system(cmd)
                 except Exception as exc:
                     print('**** Error running energy bin analysis ****')
                     print(traceback.format_exc())
                     print(exc)
                     #os.system(cmd)
            
             else : #submit a job to a cluster
                 prefix = Newconfig['out'] + "/"+ EbinPath + str(ind)
                 scriptname = prefix + "_Script.sh"
                 JobLog = prefix + "_Job.log"
                 JobName = (Newconfig['target']['name'] + "_" +
                           Newconfig['analysis']['likelihood'] +
                           "_Ebin_" + str(ind) + "_" + Newconfig['file']['tag'])
                 call(cmd, enricodir, fermidir, scriptname, JobLog, JobName)# submition
             ind+=1
