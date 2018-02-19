import os
import numpy as np
from enrico import environ
from enrico.constants import EbinPath
from enrico.submit import call
from enrico.config import get_config
from enrico import utils,Loggin

def ChangeModel(comp, E1, E2, name, Pref, Gamma):
    """Change the spectral model of a source called name
    to allow a fit between E1 and E2
    If the spectral model is PowerLaw, the prefactor is updated
    if not the model is change to PowerLaw.
    The index is frozen in all cases"""
    
    # if approximated Gamma is outside of bounds set it to limit
    Gamma_min=-5
    Gamma_max=0.5
    Gamma = min(max(Gamma_min,Gamma),Gamma_max)
    
    Eav = utils.GetE0(E1, E2)

    spectrum = comp.logLike.getSource(name).getSrcFuncs()['Spectrum']
    spectrum.getParam('Prefactor').setBounds(1e-5,1e5)
    spectrum.getParam('Prefactor').setScale(utils.fluxScale(Pref))
    spectrum.getParam('Prefactor').setValue(utils.fluxNorm(Pref))
    spectrum.getParam('Index').setBounds(Gamma_min,Gamma_max)
    spectrum.getParam('Index').setValue(Gamma)
    spectrum.getParam('Index').setFree(0)
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
    #copy the chose of the user for the enery bin computing
    config['Spectrum']['FitsGeneration'] = config['Ebin']['FitsGeneration']
    config['UpperLimit']['TSlimit'] = config['Ebin']['TSEnergyBins']
    tag = FitRunner.config['file']['tag']
    lEmax = np.log10(float(FitRunner.config['energy']['emax']))
    lEmin = np.log10(float(FitRunner.config['energy']['emin']))
    utils._log("Preparing submission of fit into energy bins")
    print(" Emin = ", float(FitRunner.config['energy']['emin']),
          " Emax = ", float(FitRunner.config['energy']['emax']),
          " Nbins = ", NEbin)

    if config['Ebin']['DistEbins'] in ['TS','mix'] and sedresult!=None:
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

    os.system("mkdir -p " + config['out'])
    paramsfile = []

    srcname = FitRunner.config['target']['name']
    if config['UpperLimit']['TSlimit']>Fit.Ts(srcname) :
        utils._log('Re-optimize', False)
        print "An upper limit has been computed. The fit need to be re-optmized"
        Fit.optimize(0)


    Pref = utils.ApproxPref(Fit, ener, srcname)
    Gamma = utils.ApproxGamma(Fit, ener, srcname)

    Model_type = Fit.model.srcs[srcname].spectrum().genericName()
    # if the model is not PowerLaw : change the model
    if not(Model_type == 'PowerLaw'):
        for comp in Fit.components:
            comp.logLike.getSource(srcname).setSpectrum("PowerLaw") #Change model
        config['target']['spectrum'] = "PowerLaw"

    xmltag_list = [""]#handle summed like analysis
    if config['ComponentAnalysis']['FrontBack'] == "yes":
        xmltag_list = ["_FRONT","_BACK"]
        mes.info("Splitting Front/Back events")
    elif config['ComponentAnalysis']['PSF'] == "yes":
        xmltag_list = ["_PSF0","_PSF1","_PSF2"]
        mes.info("Splitting PSF events")
    elif config['ComponentAnalysis']['EDISP'] == "yes":
        xmltag_list = ["_EDISP0","_EDISP1","_EDISP2","_EDISP3"]
        mes.info("Splitting EDISP events")


    for ibin in xrange(NEbin):#Loop over the energy bins
        E = utils.GetE0(ener[ibin + 1],ener[ibin])
        mes.info("Submitting # "+str(ibin)+" at energy "+str(E))
        #Update the model for the bin
        for  comp,xmltag in zip(Fit.components, xmltag_list):
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
        Gamma_max=0.5
        config['UpperLimit']['SpectralIndex'] = -min(max(Gamma_min,Gamma[ibin]),Gamma_max)

        config['file']['tag'] = tag + '_Ebin' + str(NEbin) + '_' + str(ibin)
        filename =  config['target']['name'] + "_" + str(ibin) + ".conf"
        paramsfile.append(filename)
        config.write(open(config['out'] + '/' +paramsfile[ibin], 'w')) #save the config file in a ascii file

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
                 os.system(cmd)
             else : #submit a job to a cluster
                 prefix = Newconfig['out'] + "/"+ EbinPath + str(ind)
                 scriptname = prefix + "_Script.sh"
                 JobLog = prefix + "_Job.log"
                 JobName = (Newconfig['target']['name'] + "_" +
                           Newconfig['analysis']['likelihood'] +
                           "_Ebin_" + str(ind) + "_" + Newconfig['file']['tag'])
                 call(cmd, enricodir, fermidir, scriptname, JobLog, JobName)# submition
             ind+=1
