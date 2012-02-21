import gtfunction
import fitfunction
import pyfits
import string
import numpy as np
import os
import pyLikelihood as pyLike


def SubstracFits(infile1, infile2, config):
    """Create (absolute and relative) difference images"""
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


def _log(text, line=True):
    if line:
        print('# ' + '*' * 60)
    print("# *** %s ***\n" % text)


def PrintResult(Fit, Current_Obs):
    """@todo: document me"""
    _log('Model Result')
    Result = {}
    print Fit.model
    print
    print "Source Name\tNpred\tTS"
    for src in Fit.model.srcNames:
        if Fit.Ts(src) > 5:
            print src, "\t%2.3f\t%2.3f" % (Fit.NpredValue(src), Fit.Ts(src))
    print
    print '# ' + '*' * 60
    print
    Result['Npred'] = Fit.NpredValue(Current_Obs.srcname)
    Result['TS'] = Fit.Ts(Current_Obs.srcname)
    print "Values and (MINOS) errors for " + Current_Obs.srcname
    print "TS : ", Fit.Ts(Current_Obs.srcname)

    stype = Fit.model.srcs[Current_Obs.srcname].spectrum().genericName()
    Result['ModelType'] = stype
    if stype == 'PowerLaw2':
        integral = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Integral')
        Flux = integral.value()
        ErrFlux = integral.error()
        Scale = integral.getScale()
        Result['Flux'] = Flux * Scale
        Result['dFlux'] = ErrFlux * Scale
        if Fit.Ts(Current_Obs.srcname) > 5:
            try:
                Interror = Fit.minosError(Current_Obs.srcname, 'Integral')
                print(" Integral:  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ] %2.0e" % 
                      (Flux, ErrFlux, Interror[0], Interror[1], Scale))
                Result.update({'dFlux-': Interror[0] * Scale})
                Result.update({'dFlux+': Interror[1] * Scale})
            except:
                print(" Integral:  %2.2f +/-  %2.2f  %2.0e" %
                      (Flux, ErrFlux, Scale))
        else:
            print(" Integral:  %2.2f +/-  %2.2f  %2.0e" %
                  (Flux, ErrFlux, Scale))
    elif stype == 'PowerLaw':
        spectrum = Fit[Current_Obs.srcname].funcs['Spectrum']
        prefactor = spectrum.getParam('Prefactor')
        Flux = prefactor.value()
        ErrFlux = prefactor.error()
        Scale = prefactor.getScale()
        Escale = spectrum.getParam('Scale').value()
        Result['Prefactor'] = Flux * Scale
        Result['dPrefactor'] = ErrFlux * Scale
        Result['Escale'] = Escale
        if Fit.Ts(Current_Obs.srcname) > 5:
            try:
                Interror = Fit.minosError(Current_Obs.srcname, 'Prefactor')
                print(" Prefactor:  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ] %2.0e" % 
                      (Flux, ErrFlux, Interror[0], Interror[1], Scale))
                Result['dPrefactor-'] = Interror[0] * Scale
                Result['dPrefactor+'] = Interror[1] * Scale
            except:
                print(" Prefactor:  %2.2f +/-  %2.2f %2.0e" % 
                      (Flux, ErrFlux, Scale))
        else:
            print(" Prefactor:  %2.2f +/-  %2.2f %2.0e" % 
                  (Flux, ErrFlux, Scale))

    if stype == 'PowerLaw2' or stype == 'PowerLaw':
        index = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Index')
        Gamma = index.value()
        ErrGamma = index.error()
        Result['Index'] = Gamma
        Result['dIndex'] = ErrGamma
        if Fit.Ts(Current_Obs.srcname) > 5:
            try:
                Gamerror = Fit.minosError(Current_Obs.srcname, 'Index')
                print(" Index:  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ]" %
                      (Gamma, ErrGamma, Gamerror[0], Gamerror[1]))
                Result['dIndex-'] = Gamerror[0] * Scale
                Result['dIndex+'] = Gamerror[1] * Scale
            except:
                print " Index:  %2.2f +/-  %2.2f" % (Gamma, ErrGamma)
        else:
            print " Index:  %2.2f +/-  %2.2f" % (Gamma, ErrGamma)
    return Result


def RemoveWeakSources(Fit, SourceName=None):
    """@todo: document me"""
    _log('Remove all the weak sources')
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


def GetFlux(Fit):
    """@todo: document me"""
    print "Source Flux : "
    for src in Fit.model.srcNames:
        try:
            print(src + " : %2.2e +/-  %2.2e" % 
                  (Fit.flux(src), Fit.fluxError(src)))
        except:
            pass
    print
    return Fit.flux(src)


def GetCovar(srcname, Fit):
    """Extract covariance matrix"""
    # @todo: Unused variable 'ptsrc'
    ptsrc = pyLike.PointSource_cast(Fit[srcname].src)
    par_index_map = {}
    indx = 0
    for src in Fit.sourceNames():
        parNames = pyLike.StringVector()
        Fit[src].src.spectrum().getFreeParamNames(parNames)
        for par in parNames:
            par_index_map["::".join((src, par))] = indx
            indx += 1
            # @todo: ???
            # Build the source-specific covariance matrix.
            #
    if Fit.covariance is None:
        raise RuntimeError("Covariance matrix has not been computed.")
    covar = np.array(Fit.covariance)
    if len(covar) != len(par_index_map):
        raise RuntimeError("Covariance matrix size does not match the " + 
               "number of free parameters.")
    my_covar = []
    srcpars = pyLike.StringVector()
    Fit[srcname].src.spectrum().getFreeParamNames(srcpars)
    pars = ["::".join((srcname, x)) for x in srcpars]
    for xpar in pars:
        ix = par_index_map[xpar]
        my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
    print "The covariance matrix is :\n", np.array(my_covar)
    print


def getParamIndx(fit, name, NAME):
    """Get index for a specific parameter for a specific source
    from model in UnbinnedAnalysis object fit"""
    ID = -1
    spec = fit[name].funcs['Spectrum']
    for indx, parName in zip(spec._parIds, spec.paramNames):
        if(parName == NAME):
            ID = indx
    if(ID == -1):
        print('Parameter %s not found for source %s in file %s.' %
              (NAME, name, fit.srcModel))
    return ID


def ChangeModel(Fit, Em0, Em1):
    """@todo: document me"""
    # @todo: call utility function
    E0 = int(pow(10, (np.log10(Em1) + np.log10(Em0)) / 2))

    for name in Fit.model.srcNames:
        generic_name = Fit.model.srcs[name].spectrum().genericName()
        if generic_name == 'PowerLaw':
            IdPref = getParamIndx(Fit, name, 'Prefactor')
            IdEScale = getParamIndx(Fit, name, 'Scale')
            Flux = Fit[IdPref].value()
            Scale = Fit[IdPref].getScale()
            Escale = Fit[IdEScale].value()
            IdGamma = getParamIndx(Fit, name, 'Index')
            Gamma = Fit[IdGamma].value()
            Fit[IdGamma].setFree(0)
            NewFlux = Flux * pow(E0 / Escale, Gamma) * Scale
            NormFlux = pow(10, np.log10(NewFlux) - int(np.log10(NewFlux))) * 10
            NewScale = pow(10, int(np.log10(NewFlux)) - 1)
            print "NormFlux, Scale, Gamma:"
            print NormFlux, " ", Scale, " ", Gamma
            Fit[IdEScale].setBounds(0.0, 4e5)
            Fit[IdPref].setBounds(0, Flux * 1000)
            Fit[IdPref].setScale(NewScale)
            Fit[IdPref] = NormFlux
            Fit[IdPref].setBounds(NormFlux * 0.05, NormFlux * 50)
            Fit[IdEScale] = E0
            Fit[IdEScale].setBounds(E0 * 0.05, E0 * 50)
        elif generic_name == 'PowerLaw2':
            IdInt = getParamIndx(Fit, name, 'Integral')
            IdEmin = getParamIndx(Fit, name, 'LowerLimit')
            IdEmax = getParamIndx(Fit, name, 'UpperLimit')
            IdGamma = getParamIndx(Fit, name, 'Index')
            Gamma = Fit[IdGamma].value()
            Fit[IdGamma].setFree(0)
            Flux = Fit[IdInt].value()
            Scale = Fit[IdInt].getScale()
            Emin = Fit[IdEmin].value()
            Emax = Fit[IdEmax].value()
            D = pow(Em1, Gamma + 1) - pow(Em0, Gamma + 1)
            N = pow(Emax, Gamma + 1) - pow(Emin, Gamma + 1)
            NewFlux = Flux * D / N * Scale
            NormFlux = pow(10, np.log10(NewFlux) - int(np.log10(NewFlux))) * 10
            NewScale = pow(10, int(np.log10(NewFlux)) - 1)
            Fit[IdInt].setBounds(0, Flux * 1000)
            Fit[IdInt].setScale(NewScale)
            Fit[IdInt] = NormFlux
            Fit[IdInt].setBounds(NormFlux * 0.05, NormFlux * 50)
            Fit[IdEmin] = Em0
            Fit[IdEmax] = Em1
            print "NormFlux, NewScale, Gamma"
            print NormFlux, " ", NewScale, " ", Gamma
    return Fit


def Analysis(folder, config, tag="", convtyp='-1'):
    Obs = gtfunction.Observation(folder, config, convtyp, tag=tag)
    _log('SUMMARY' + tag)
    Obs.printSum()
    runfit = fitfunction.MakeFit(Obs, config)
    if config['Spectrum']['FitsGeneration'] == 'yes':
        runfit.PreparFit()
    return runfit, Obs


def PrepareEbin(Fit, runfit):
    NEbin = int(runfit.config['Ebin']['NumEnergyBins'])
    NewConfig = runfit.config
    NewConfig['UpperLimit']['envelope'] = 'no'
    NewConfig['Ebin']['NumEnergyBins'] = '0'
    NewConfig['out'] = runfit.config['out'] + '/Ebin'
    NewConfig['Spectrum']['ResultPlots'] = 'no'
    NewConfig['Spectrum']['FitsGeneration'] = 'yes'
    NewConfig['UpperLimit']['TSlimit'] = NewConfig['Ebin']['TSEnergyBins']
    tag = runfit.config['file']['tag']
    lEmax = np.log10(float(runfit.config['energy']['emax']))
    lEmin = np.log10(float(runfit.config['energy']['emin']))
    print("Preparing submission of fit into energy bins")
    print("Emin = ", float(runfit.config['energy']['emin']),
          " Emax = ", float(runfit.config['energy']['emax']),
          " Nbins = ", NEbin)
    ener = np.logspace(lEmin, lEmax, NEbin + 1)
    os.system("mkdir -p " + runfit.config['out'] + '/Ebin')
    paramsfile = []

    RemoveWeakSources(Fit)

    for ibin in xrange(NEbin):
        E = int(pow(10, (np.log10(ener[ibin + 1]) + np.log10(ener[ibin])) / 2))
        print "Submition # ", ibin, " at energy ", E
        ChangeModel(Fit, ener[ibin], ener[ibin + 1])
        NewModel = NewConfig['out'] + "/" + runfit.config['target']['name'] + "_" + str(E) + ".xml"
        Fit.writeXml(NewModel)
        NewConfig['file']['xml'] = NewModel
        NewConfig['Spectrum']['FitsGeneration'] = NewConfig['Ebin']['FitsGeneration']
        NewConfig['energy']['emin'] = str(ener[ibin])
        NewConfig['energy']['emax'] = str(ener[ibin + 1])
        NewConfig['file']['tag'] = tag + '_Ebin_' + str(ibin)
        filename = (runfit.config['out'] + '/' + 
                    runfit.config['target']['name'] + "_" + str(E) + ".conf")
        paramsfile.append(filename)
        NewConfig.write(open(paramsfile[ibin], 'w'))

    return paramsfile


def _dump_filename(config):
    """@todo: document me"""
    return (config['out'] + '/' + config['target']['name'] + '_' +
            str(int(config['time']['tmin'])) + '_' +
            str(int(config['time']['tmax'])) + '_' +
            str(int(config['energy']['emin'])) + '_' +
            str(int(config['energy']['emax'])) + ".results")


def DumpResult(Result, config):
    """@todo: Use configobj, not this hand-written utility function!"""
    Dumpfile = open(_dump_filename(config), "w")
    for key in Result.iterkeys():
        Dumpfile.write(key + '\t' + str(Result[key]) + '\n')
    Dumpfile.close()


def ReadResult(config):
    """@todo: Use configobj, not this hand-written utility function!"""
    lines = open(_dump_filename(config)).readlines()
    results = dict()
    for line in lines:
        key, value = string.split(line)[0:2]
        try:
            value = float(value)
        except:
            pass
        results[key] = value
    return results
