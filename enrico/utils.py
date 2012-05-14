"""Random collection of usefull functions"""
import os
import numpy as np
from math import log10

def _log(text, line=True):
    if line:
        print('# ' + '*' * 60)
    print("# *** %s ***" % text)
    if line:
        print('# ' + '*' * 60)

def fluxScale(flux_value):
    """Get the scale of the flux value
    ex : 1.4e-14 ---> 1e-14"""
    return 10 ** np.floor(np.log10(flux_value) + 0.5)

def fluxNorm(flux_value):
    """Return the norm from the flux_value
     ex : 1.4e-14 ---> 1.4"""
    return pow(10, np.log10(flux_value) - int(np.log10(flux_value))) * 10

def Prefactor(flux,index,emin,emax,escale):
    """Compute the prefactor at the energy escale knowing
    the flux and index between emin and emax"""

    Denomin = pow(emax,-abs(index)+1.) -pow(emin,-abs(index)+1.)
    return flux*(-abs(index)+1)*pow(escale,-abs(index)) / Denomin

def dNde(energy,Fit,name):
    '''Compute the dN/dE value at energy E fir the source name'''
    import pyLikelihood
    ptsrc = pyLikelihood.PointSource_cast(Fit[name].src)
    arg = pyLikelihood.dArg(energy)
    return ptsrc.spectrum()(arg)

def meanEnergy(emin, emax, index_value):
    """Get the mean energy, weighted with a power law of a given index"""
    x = emax / emin
    if index_value == -2.0:
        eflux = emax * np.log(x) / (x - 1)
    elif index_value == -1.0:
        eflux = emin * (x - 1) / np.log(x)
    else:
        factor1 = emin * (index_value + 1) / (index_value + 2)
        factor2 = (x ** (index_value + 2) - 1) / (x ** (index_value + 1) - 1)
        eflux = factor1 * factor2
    return eflux

def GetE0(a,b):
    """"Get the center of a bin in log space"""
    return int(pow(10, (np.log10(a) + np.log10(b)) / 2))

def calcAngSepDeg(ra0, dec0, ra1, dec1):
    '''Return the angular separation between two objects. Use the
    special case of the Vincenty formula that is accurate for all
    distances'''
    C = np.pi / 180
    d0 = C * dec0
    d1 = C * dec1
    r12 = C * (ra0 - ra1)
    cd0 = np.cos(d0)
    sd0 = np.sin(d0)
    cd1 = np.cos(d1)
    sd1 = np.sin(d1)
    cr12 = np.cos(r12)
    sr12 = np.sin(r12)
    num = np.sqrt((cd0 * sr12) ** 2 + (cd1 * sd0 - sd1 * cd0 * cr12) ** 2)
    den = sd0 * sd1 + cd0 * cd1 * cr12
    return np.arctan2(num, den) / C

def etag(emin, emax, fmt='%07d'):
    return ('emin_%s_emax_%s' % (fmt, fmt)) % (emin, emax)

def cube_to_image(cube, slicepos=None, mean=False):
    """ Make an image out of a cube.
    Both in- and output shoud by pyfits.HDUs"""
    from pyfits import PrimaryHDU
    header = cube.header.copy()
    header['NAXIS'] = 2
    del header['NAXIS3']
    del header['CRVAL3']
    del header['CDELT3']
    if slicepos:
        data = cube.data[slicepos]
    else:
        if mean:
            data = cube.data.mean(0).astype(cube.data.dtype)
        else:
            data = cube.data.sum(0).astype(cube.data.dtype)

    return PrimaryHDU(data, header)


def SubstracFits(infile1, infile2, config):
    """Create (absolute and relative) difference images"""
    import pyfits
    data1 = pyfits.getdata(infile1)
    data2 = pyfits.getdata(infile2)
    head = pyfits.getheader(infile2)
    filebase = config['out'] + "/" + config['target']['name']
    abs_diff_file = filebase + "_Substract_Model_cmap.fits"
    rel_diff_file = filebase + "_Residual_Model_cmap.fits"
    os.system("rm " + abs_diff_file)
    os.system("rm " + rel_diff_file)
    pyfits.writeto(abs_diff_file, data1 - data2, head)
    pyfits.writeto(rel_diff_file, (data1 - data2) / data2, head)

def Analysis(folder, config, tag="", convtyp='-1'):
    """ run an analysis"""
    from gtfunction import Observation
    from fitfunction import MakeFit
    Obs = Observation(folder, config, convtyp, tag=tag)
    _log('SUMMARY' + tag)
    Obs.printSum()
    runfit = MakeFit(Obs, config)
    if config['Spectrum']['FitsGeneration'] == 'yes':
        runfit.PreparFit()
    return runfit

def RemoveWeakSources(Fit, SourceName=None):
    """Remove the weak source after a fit and reoptimized
     weak mens TS<1"""
    _log('Remove all the weak (TS<1) sources')
    NoWeakSrcLeft = False
    while not(NoWeakSrcLeft):
        NoWeakSrcLeft = True
        for src in Fit.model.srcNames:
            if Fit.Ts(src) < 1 and not(src == SourceName):
                print "delete source : ", src
                NoWeakSrcLeft = False
                Fit.deleteSource(src)
        if not(NoWeakSrcLeft):
            _log('Re-optimize', False)
            Fit.optimize(0)
        print
    return Fit


def GetFluxes(Fit,Emin=1e2,Emax=3e5):
    """Print the integral flux and error for all the sources"""
    print "Source Flux  [%2.2e MeV, %2.2e MeV] : " %(Emin,Emax)
    for src in Fit.model.srcNames:
        try:
            print(src + "   Integral Flux : %2.2e +/-  %2.2e ph/cm2/s" %
                  (Fit.flux(src,Emin,Emax), Fit.fluxError(src,Emin,Emax)))
        except:
            pass
    print

def GetCovar(srcname, Fit, verbose = True):
    """Extract covariance matrix"""
    import pyLikelihood
    par_index_map = {}
    indx = 0
    for src in Fit.sourceNames():
        parNames = pyLikelihood.StringVector()
        Fit[src].src.spectrum().getFreeParamNames(parNames)
        for par in parNames:
            par_index_map["::".join((src, par))] = indx
            indx += 1
    if Fit.covariance is None:
        raise RuntimeError("Covariance matrix has not been computed.")
    covar = np.array(Fit.covariance)
    if len(covar) != len(par_index_map):
        raise RuntimeError("Covariance matrix size does not match the " +
               "number of free parameters.")
    my_covar = []
    srcpars = pyLikelihood.StringVector()
    Fit[srcname].src.spectrum().getFreeParamNames(srcpars)
    pars = ["::".join((srcname, x)) for x in srcpars]
    for xpar in pars:
        ix = par_index_map[xpar]
        my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
    if verbose :
      print "The covariance matrix is :\n", np.array(my_covar)
      print

    return my_covar

def getParamIndx(fit, name, parameter):
    """Get index for a specific parameter for a specific source
    from model in UnbinnedAnalysis object fit"""
    ID = -1
    spec = fit[name].funcs['Spectrum']
    for indx, parName in zip(spec._parIds, spec.paramNames):
        if(parName == parameter):
            ID = indx
    if(ID == -1):
        print('Parameter %s not found for source %s in file %s.' %
              (parameter, name, fit.srcModel))
    return ID

def ChangeModel(Fit, E1, E2, name, model_type):
    """Change the spectral model of a source called name
    to allow a fit between E1 and E2
    If the spectral model is PowerLaw, the prefactor is updated
    if not the model is change to PowerLaw.
    The index is frozen in all case"""

    Eav = GetE0(E2,E1)
#    flux = Fit.flux(name,E1,E2) #Source flux between E2 and Em2

    # if the model is not PowerLaw : change the model
    if model_type == 'PowerLaw' :
      IdGamma = getParamIndx(Fit, name, 'Index')
      Gamma = Fit[IdGamma].value()
      Pref = dNde(Eav,Fit,name)
    else :
      #Compute an approximation of an index
      dnde1 = log10(dNde(E1,Fit,name))
      dnde2 = log10(dNde(E2,Fit,name))
      Pref = dNde(Eav,Fit,name)
      Gamma = (dnde2-dnde1)/(log10(1.*E2)-log10(1.*E1))
      Fit.logLike.getSource(name).setSpectrum("PowerLaw") #Change model

    # Set Parameters
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-5,1e5)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(fluxScale(Pref))
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Prefactor').setValue(fluxNorm(Pref))

    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setBounds(-5,0)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setValue(Gamma)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Index').setFree(0)

    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setValue(Eav)
    Fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(20,3e6)

    return Fit

def PrepareEbin(Fit, runfit):
    """ Prepare the computation of spectral point in energy bins by
    i) removeing the weak sources (TS<1)
    ii) updating the config file (option and energy) 
    and save it in a new ascii file
    iii) changing the spectral model and saving it in a new xml file.
    A list of the ascii files is returned"""

    NEbin = int(runfit.config['Ebin']['NumEnergyBins'])

    config = runfit.config

    #Replace the evt file with the fits file produced before
    #in order to speed up the production of the fits files
    config['file']['event'] = runfit.obs.eventfile
    #update the config to allow the fit in energy bins
    config['UpperLimit']['envelope'] = 'no' 
    config['Ebin']['NumEnergyBins'] = '0'#no new bin in energy!
    config['out'] = runfit.config['out'] + '/Ebin' + str(NEbin)
    config['Spectrum']['ResultPlots'] = 'no' #no SED plot/modelmap
    #copy the chose of the user for the enery bin computing
    config['Spectrum']['FitsGeneration'] = config['Ebin']['FitsGeneration'] 
    config['UpperLimit']['TSlimit'] = config['Ebin']['TSEnergyBins']
    tag = runfit.config['file']['tag']
    lEmax = np.log10(float(runfit.config['energy']['emax']))
    lEmin = np.log10(float(runfit.config['energy']['emin']))
    print("Preparing submission of fit into energy bins")
    print("Emin = ", float(runfit.config['energy']['emin']),
          " Emax = ", float(runfit.config['energy']['emax']),
          " Nbins = ", NEbin)

    ener = np.logspace(lEmin, lEmax, NEbin + 1)
    os.system("mkdir -p " + config['out'])
    paramsfile = []

    RemoveWeakSources(Fit)#remove source with TS<1 to be sure that MINUIT will converge
    srcname = runfit.config['target']['name']
    Model_type = Fit.model.srcs[srcname].spectrum().genericName()

    for ibin in xrange(NEbin):#Loop over the energy bins
        E = GetE0(ener[ibin + 1],ener[ibin])
        print "Submition # ", ibin, " at energy ", E
        #Update the model for the bin
        ChangeModel(Fit, ener[ibin], ener[ibin + 1], srcname, Model_type)
        Xmlname = (config['out'] + "/" + srcname +
                    "_" + str(ibin) + ".xml")
        Fit.writeXml(Xmlname)# dump the corresponding xml file
        config['file']['xml'] = Xmlname
        #update the energy bounds
        config['energy']['emin'] = str(ener[ibin])
        config['energy']['emax'] = str(ener[ibin + 1])
        config['file']['tag'] = tag + '_Ebin' + str(NEbin) + '_' + str(ibin)
        filename =  config['target']['name'] + "_" + str(ibin) + ".conf"
        paramsfile.append(filename)
        config.write(open(config['out'] + '/' +paramsfile[ibin], 'w')) #save the config file in a ascii file

    return paramsfile

def _dump_xml(config) :
    """Give the name of the XML file where the results will be save by gtlike"""
    return (config['out'] + "/" + config['target']['name'] + "_" +
                  config['file']['tag'] + "_out.xml")

def _dump_filename(config):
    """Give the name of the file where the results will be dumped"""
    return (config['out'] + '/' + config['target']['name'] + '_' +
            str(config['target']['spectrum']) + '_' +
            str(int(config['time']['tmin'])) + '_' +
            str(int(config['time']['tmax'])) + '_' +
            str(int(config['energy']['emin'])) + '_' +
            str(int(config['energy']['emax'])) + ".results")


def DumpResult(Result, config):
    """Dump the result into an ascii file """
    Dumpfile = open(_dump_filename(config), "w")
    for key in Result.iterkeys():
        Dumpfile.write(key + '\t' + str(Result[key]) + '\n')
        print key + '\t' + str(Result[key])
    Dumpfile.close()


def ReadResult(config):
    """Read the result from an ascii file """
    lines = open(_dump_filename(config)).readlines()
    results = dict()
    for line in lines:
        key, value = line.split()[0:2]
        try:
            value = float(value)
        except:
            pass
        results[key] = value
    return results
