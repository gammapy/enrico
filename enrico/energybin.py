#!/usr/bin/env python
import utils
import numpy as np
import environ,os
from submit import call
from config import get_config

def ChangeModel(Fit, E1, E2, name, Pref, Gamma):
    """Change the spectral model of a source called name
    to allow a fit between E1 and E2
    If the spectral model is PowerLaw, the prefactor is updated
    if not the model is change to PowerLaw.
    The index is frozen in all cases"""

    Eav = utils.GetE0(E1, E2)

    # Set Parameters
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-5,1e5)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(utils.fluxScale(Pref))
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(utils.fluxNorm(Pref))

    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setBounds(-5,0)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setValue(Gamma)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setFree(0)

    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setValue(Eav)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(20,3e6)

    return Fit

def PrepareEbin(Fit, FitRunner):
    """ Prepare the computation of spectral point in energy bins by
    i) removing the weak sources (TS<1) # not true
    ii) updating the config file (option and energy) 
    and save it in a new ascii file
    iii) changing the spectral model and saving it in a new xml file.
    A list of the ascii files is returned"""

    NEbin = int(FitRunner.config['Ebin']['NumEnergyBins'])

    config = FitRunner.config

    config['verbose'] ='no' #Be quiet
    #Replace the evt file with the fits file produced before
    #in order to speed up the production of the fits files
    config['file']['event'] = FitRunner.obs.eventfile
    #update the config to allow the fit in energy bins
    config['UpperLimit']['envelope'] = 'no' 
    config['Ebin']['NumEnergyBins'] = '0'#no new bin in energy!
    config['out'] = FitRunner.config['out'] + '/Ebin' + str(NEbin)
    config['Spectrum']['ResultPlots'] = 'no' #no SED plot/modelmap
    #copy the chose of the user for the enery bin computing
    config['Spectrum']['FitsGeneration'] = config['Ebin']['FitsGeneration'] 
    config['UpperLimit']['TSlimit'] = config['Ebin']['TSEnergyBins']
    tag = FitRunner.config['file']['tag']
    lEmax = np.log10(float(FitRunner.config['energy']['emax']))
    lEmin = np.log10(float(FitRunner.config['energy']['emin']))
    utils._log("Preparing submission of fit into energy bins")
    print("Emin = ", float(FitRunner.config['energy']['emin']),
          " Emax = ", float(FitRunner.config['energy']['emax']),
          " Nbins = ", NEbin)

    ener = np.logspace(lEmin, lEmax, NEbin + 1)
    os.system("mkdir -p " + config['out'])
    paramsfile = []

    srcname = FitRunner.config['target']['name']
    if config['UpperLimit']['TSlimit']>Fit.Ts(srcname) :
        utils._log('Re-optimize', False)
        print "An upper limit has been computed. The fit need to be re-optmized"
        Fit.optimize(0)

#    utils.RemoveWeakSources(Fit,srcname)#remove source with TS<1 to be sure that MINUIT will converge

    Pref = utils.ApproxPref(Fit, ener, srcname)
    Gamma = utils.ApproxGamma(Fit, ener, srcname)

    Model_type = Fit.model.srcs[srcname].spectrum().genericName()
    # if the model is not PowerLaw : change the model
    if not(Model_type == 'PowerLaw') :
      Fit.logLike.getSource(srcname).setSpectrum("PowerLaw") #Change model

    for ibin in xrange(NEbin):#Loop over the energy bins
        E = utils.GetE0(ener[ibin + 1],ener[ibin])
        print "Submition # ", ibin, " at energy ", E
        #Update the model for the bin
        NewFitObject = ChangeModel(Fit, ener[ibin], ener[ibin + 1], srcname, Pref[ibin] ,Gamma[ibin])
        Xmlname = (config['out'] + "/" + srcname +
                    "_" + str(ibin) + ".xml")
        NewFitObject.writeXml(Xmlname)# dump the corresponding xml file
        config['file']['xml'] = Xmlname
        #update the energy bounds
        config['energy']['emin'] = str(ener[ibin])
        config['energy']['emax'] = str(ener[ibin + 1])
        config['file']['tag'] = tag + '_Ebin' + str(NEbin) + '_' + str(ibin)
        filename =  config['target']['name'] + "_" + str(ibin) + ".conf"
        paramsfile.append(filename)
        config.write(open(config['out'] + '/' +paramsfile[ibin], 'w')) #save the config file in a ascii file

    return paramsfile


def RunEbin(folder,Nbin,Fit,FitRunner):
    if int(Nbin) > 0:
        configfiles = PrepareEbin(Fit, FitRunner)
        ind = 0
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')
        for conf in configfiles:
             pathconf = folder + "/Ebin" + str(Nbin) +"/" + conf
             Newconfig = get_config(pathconf)
             cmd = enricodir+"/enrico/RunGTlike.py "+pathconf
             if Newconfig['Submit'] == 'no' : #run directly
                 os.system(cmd)
             else : #submit a job to a cluster
                 prefix = Newconfig['out'] + "/Ebin" + str(ind) 
                 scriptname = prefix + "_Script.sh"
                 JobLog = prefix + "_Job.log"
                 JobName = (Newconfig['target']['name'] + "_" +
                           Newconfig['analysis']['likelihood'] +
                           "_Ebin_" + str(ind) + "_" + Newconfig['file']['tag'])
                 call(cmd, enricodir, fermidir, scriptname, JobLog, JobName)# submition
             ind+=1
