"""Random collection of usefull functions"""
import os
from math import log10
import numpy as np
from enrico.constants import met_ref,mjd_ref, jd_ref, DAY_IN_SECOND

typeirfs={1:"FRONT",2:"BACK",3:"",4:"PSF0",8:"PSF1",16:"PSF2",32:"PSF3",64:"EDISP0",
        128:"EDISP1",256:"EDISP2",512:"EDISP3"}

def hasKey(dictionary,key):
    try:             Item = dictionary[key]
    except KeyError: return(False)
    else:            return(True)

def isKey(dictionary,key):
    if not hasKey(dictionary,key): return(None)
    return(dictionary[key])

def _log(text, line=True):
    if line:
        print("\033[34m"+'# ' + '*' * 60)
    print("\033[34m"+"# *** %s ***" % text)
    if line:
        print "\033[34m"+'# ' + '*' * 60+"\033[0m"

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


def SubtractFits(infile1, infile2, config):
    """Create (absolute and relative) difference images"""
    import astropy.io.fits as pyfits
    data1 = pyfits.getdata(infile1)
    data2 = pyfits.getdata(infile2)
    head = pyfits.getheader(infile2)
    filebase = config['out'] + "/" + config['target']['name']
    abs_diff_file = filebase + "_Subtract_Model_cmap.fits"
    rel_diff_file = filebase + "_Residual_Model_cmap.fits"
    os.system("rm " + abs_diff_file)
    os.system("rm " + rel_diff_file)
    pyfits.writeto(abs_diff_file, data1 - data2, head)
    pyfits.writeto(rel_diff_file, (data1 - data2) / data2, head)


def GetFluxes(Fit,Emin=1e2,Emax=3e5):
    """Print the integral flux and error for all the sources"""
    print "\nSource Flux  [%2.2e MeV, %2.2e MeV] : " %(Emin,Emax)
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
    # Try to get the first component (if summed analysis)
    try:    fit = fit.components[0]
    except: pass
    spec = fit[name].funcs['Spectrum']
    for indx, parName in zip(spec._parIds, spec.paramNames):
        if(parName == parameter):
            ID = indx
    if(ID == -1):
        print('Parameter %s not found for source %s in file %s.' %
             (parameter, name, fit.srcModel))

    return ID

def FreezeParams(fit, name, parameter, value):
    try:
        fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam(parameter).setBounds(value,value)
        fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam(parameter).setValue(value)
        fit.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam(parameter).setFree(0)
    except:
        for comp in self.components:
            try:
                comp.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam(parameter).setBounds(value,value)
                comp.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam(parameter).setValue(value)
                comp.logLike.getSource(name).getSrcFuncs()['Spectrum'].getParam(parameter).setFree(0)
            except:
                continue

def ApproxPref(Fit, ener,name):
    Pref = np.zeros(len(ener)-1)
    for ibin in xrange(len(ener)-1):
      Eav = GetE0(ener[ibin+1],ener[ibin])
      Pref[ibin] = dNde(Eav,Fit,name)

    return Pref

def ApproxGamma(Fit, ener,name):
    """ Get an approximation of the index for different bin in energy"""
    Gamma = np.zeros(len(ener)-1)
    for ibin in xrange(len(ener)-1):
      #Compute an approximation of an index
      dnde1 = np.log10(dNde(ener[ibin],Fit,name))
      dnde2 = np.log10(dNde(ener[ibin+1],Fit,name))
      Gamma[ibin] = (dnde2-dnde1)/(np.log10(1.*ener[ibin+1])-np.log10(1.*ener[ibin]))

    return Gamma


def _SpecFileName(config):
    """return a generic name for the file related to the spectrum (plot, results...)"""
    from enrico.constants import SpectrumPath
    return  config['out'] + '/'+SpectrumPath+'/SED_' + config['target']['name'] +'_'+ config['target']['spectrum']

def _dump_xml(config) :
    """Give the name of the XML file where the results will be save by gtlike"""
    return (config['out'] + "/" + config['target']['name']
                  + "_" + config['target']['spectrum'] + "_"+
                  config['file']['tag'] + "_out.xml")


def _dump_reg(config):
    return config['out'] + "/Roi_model.reg"

def _dump_findsrcout(config):
    res = config['out'] + "/FindSource.out"
    os.system("touch "+res)
    return res


def _dump_filename(config):
    """Give the name of the file where the results will be dumped"""
    return (config['out'] + '/' + config['target']['name'] + '_' +
            str(config['target']['spectrum']) + '_' +
            str(int(0.5+config['time']['tmin'])) + '_' +
            str(int(0.5+config['time']['tmax'])) + '_' +
            str(int(0.5+config['energy']['emin'])) + '_' +
            str(int(0.5+config['energy']['emax'])) + "_"+
                  config['file']['tag'] +  ".results")


def DumpResult(Result, config):
    """Dump the result into an ascii file """
    Dumpfile = open(_dump_filename(config), "w")
    for key in Result.iterkeys():
        Dumpfile.write(key + '\t' + str(Result[key]) + '\n')
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

def time_selection_string(config,numbin0):
    """Convert file with start stop pairs to gtmktime filter string"""

    if numbin0==None:
        numbin0=0

    # Read MET_TSTART, MET_TSTOP pairs from file
    bins = np.loadtxt(config['time']['file'])

    if config['time']['type']=='MJD':
        bins = MJD_to_met(bins)
    elif config['time']['type']=='JD':
        bins = JD_to_met(bins)

    selstr=''
    last=True

    if np.shape(bins)==(2,): bins = [bins]
    for numbin in range(numbin0,len(bins)):
        tbin=bins[numbin]
        selstr+='((START>{0:.0f})&&(STOP<{1:.0f}))||'.format(tbin[0],tbin[1])
        if len(selstr)>=800:
            last=False
            break

    # remove last ||, and enclose in parens
    selstr='('+selstr[:-2]+')'
    # add one to numbin so that on next call it starts on the following bin to the last one that was included in selstr
    return selstr, numbin0+1, last

def met_to_MJD(met):
  return mjd_ref + (met-met_ref)/DAY_IN_SECOND

def MJD_to_met(mjd):
  return met_ref + (mjd-mjd_ref)*DAY_IN_SECOND

def JD_to_met(jd):
  return MJD_to_met(mjd)+2400000.5

def create_dir(path):
    import os
    import os.path

    if (not os.path.exists(path)):
        os.makedirs(path)

def Checkevtclass(evclass):
    classirfs = {1:"P8R3_TRANSIENT100A",2:"P8R3_TRANSIENT100E",4:"P8R3_TRANSIENT100",8:"P8R3_TRANSIENT020E",
			16:"P8R3_TRANSIENT020",32:"P8R3_TRANSIENT010E",64:"P8R3_TRANSIENT010",128:"P8R3_SOURCE",
			256:"P8R3_CLEAN",521:"P8R3_ULTRACLEAN",1024:"P8R3_ULTRACLEANVETO",32768:"P8R3_TRANSIENT100S",
			65536:"P8R3_TRANSIENT015S",16777216:"P8R3_LLE"}
    try :
        tmp = classirfs[evclass]
    except:
        from enrico import Loggin
        mes = Loggin.Message()
        mes.error("evclass value in config file not valid")

def GetSDC(val):
    deno = 0
    while val>=2:
        val = val/2
        deno += 1
    return deno

def GetIRFS(evtclass,evttype,addversion=True):
    from enrico import Loggin
    mes = Loggin.Message()
    classirfs = {1:"P8R3_TRANSIENT100A",2:"P8R3_TRANSIENT100E",4:"P8R3_TRANSIENT100",8:"P8R3_TRANSIENT020E",
			16:"P8R3_TRANSIENT020",32:"P8R3_TRANSIENT010E",64:"P8R3_TRANSIENT010",128:"P8R3_SOURCE",
			256:"P8R3_CLEAN",521:"P8R3_ULTRACLEAN",1024:"P8R3_ULTRACLEANVETO",32768:"P8R3_TRANSIENT100S",
			65536:"P8R3_TRANSIENT015S",16777216:"P8R3_LLE"}


    result = []
    val = evttype
    while val>0 :
       deno = GetSDC(val)
       result.append(2**deno)
       val = val - result[-1]

    typ = []
    for t in result:
        typ.append(typeirfs[t])

    # P8R3_SOURCE_V2 is the irf, but iso_P8R3_SOURCE_V2_()_V2.txt does not exist, 
    # instead it is _V6_()_V2.txt. We need to get around this inconsistency.
    if (addversion):
        classirf = classirfs[evtclass]+"_V2"
    else:
        classirf = classirfs[evtclass]
    #mes.info("Using IRFS for: class %s and type %s" %(str(classirf),str(typ)))
    return classirf,typ


def GetIso(evtclass,evttype):
    irfs = GetIRFS(evtclass,evttype,addversion=False)
    import enrico.environ as e
    if len(irfs[1])> 1:
        res = os.path.join(e.DIFFUSE_DIR,'iso_'+str(irfs[0])+'_V2.txt')
    else:
        res = os.path.join(e.DIFFUSE_DIR,'iso_'+str(irfs[0])+'_'+str(irfs[1][0])+'_V2.txt')
    return res


def GetZenithCut(evtclass,evttype,emin):
    irfs = GetIRFS(evtclass,evttype)
    ener = np.array([50,100,200,300,500])
    emin_ind = sum(ener-0.1<emin)-1
    All = [80,90 ,95 ,100 ,100 ]
    #FRONT+BACK, EDISP0-EDISP3	80	90 	95 	100 	100
    FRONT = [85 ,95 ,100 ,100 ,100]
    BACK = [75,85 ,90 ,95 ,100]
    PSF0 = [70,80 ,85 ,90 ,95]
    PSF1 = [75,85 ,95 ,100 ,100]
    PSF2 = [85,95,100,100,100]
    PSF3 = [90,100,100,100,100]
    if len(irfs[1])>1 :
       return All[emin_ind]
    elif irfs[1][0] == "FRONT":
       return FRONT[emin_ind]
    elif irfs[1][0] == "BACK":
       return BACK[emin_ind]
    elif irfs[1][0] == "PSF0":
       return PSF0[emin_ind]
    elif irfs[1][0] == "PSF1":
       return PSF1[emin_ind]
    elif irfs[1][0] == "PSF2":
       return PSF2[emin_ind]
    elif irfs[1][0] == "PSF3":
       return PSF3[emin_ind]
