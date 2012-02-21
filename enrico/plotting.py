import os
import string
import array  # @todo: get rid of array, use numpy instead
import ROOT
import numpy as np
import pyfits
import root_style
import pyLikelihood
from config import get_config
import utils


class Params:
    """Collection of plot parameters"""
    def __init__(self, srcname, Emin=100, Emax=3e5, extend=False,
                 PlotName="LAT_SED", LineColor=1, AreaColor=2, plotRes=False):
        self.Emin = Emin
        self.Emax = Emax
        self.N = 500
        self.srcname = srcname
        self.extend = extend
        self.PlotName = PlotName
        self.LineColor = LineColor
        self.AreaColor = AreaColor
        self.plotRes = plotRes


class Result:
    """Helper class to get at results"""
    Prefac_UL = 0.

    def __init__(self, like, pars):
        self.like = like
        self.ra = like[pars.srcname].funcs['Position'].getParam('RA').value()
        self.dec = like[pars.srcname].funcs['Position'].getParam('DEC').value()
        self.Model = like[pars.srcname].funcs['Spectrum'].genericName()
        try:
            self.TS = like.Ts(pars.srcname)
        except RuntimeError:
            self.TS = -1
        self.ptsrc = pyLikelihood.PointSource_cast(like[pars.srcname].src)
        par_index_map = {}
        indx = 0
        for src in like.sourceNames():
            parNames = pyLikelihood.StringVector()
            like[src].src.spectrum().getFreeParamNames(parNames)
            for par in parNames:
                par_index_map["::".join((src, par))] = indx
                indx += 1
        # Build the source-specific covariance matrix.
        if like.covariance is None:
            raise RuntimeError("Covariance matrix has not been computed.")
        covar = np.array(like.covariance)
        if len(covar) != len(par_index_map):
            raise RuntimeError("Covariance matrix size does not match the " +
                   "number of free parameters.")
        my_covar = []
        srcpars = pyLikelihood.StringVector()
        like[pars.srcname].src.spectrum().getFreeParamNames(srcpars)
        pars = ["::".join((pars.srcname, x)) for x in srcpars]
        for xpar in pars:
            ix = par_index_map[xpar]
            my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
        self.covar = np.array(my_covar)
        self.srcpars = srcpars


def MakeError(result, pars):
    """@todo: document me"""
    estep = np.log(pars.Emax / pars.Emin) / (pars.N - 1)
    energies = pars.Emin * np.exp(estep * np.arange(np.float(pars.N)))
    err = np.array(pars.N * [0.])
    j = 0
    for ene in energies:
        arg = pyLikelihood.dArg(ene)
        partials = np.array(len(result.srcpars) * [0.])
        for i in xrange(len(result.srcpars)):
            x = result.srcpars[i]
            partials[i] = result.ptsrc.spectrum().derivByParam(arg, x)
        err[j] = np.sqrt(np.dot(partials, np.dot(result.covar, partials)))
        j += 1
    # @todo: document constant
    return 1.602e-6 * energies ** 2 * err


def dNde(result, energy):
    arg = pyLikelihood.dArg(energy)
    return result.ptsrc.spectrum()(arg)


def MakeFlux(result, Params):
    """Compute differential Flux distribution and
    corresponding energy and return a numpy array"""
    E = np.logspace(np.log10(Params.Emin), np.log10(Params.Emax), Params.N)
    Flux = np.array(Params.N * [0.])
    for i in xrange(Params.N):
        Flux[i] = dNde(result, E[i])
    return E, Flux


def MakeFluxInBin(result, pars, Emin_bin, Emax_bin):
    """Compute differential Flux distribution and
    corresponding energy and return a numpy array"""
    E = np.logspace(np.log10(Emin_bin), np.log10(Emax_bin), pars.N)
    Flux = np.array(pars.N * [0.])
    for i in xrange(pars.N):
        Flux[i] = dNde(result, E[i])
    return E, Flux


def MakeSED(result, Params):
    """Compute Spectral energy distribution and corresponding energy
    and return a numpy array"""
    E = np.logspace(np.log10(Params.Emin), np.log10(Params.Emax), Params.N)
    nuFnu = np.array(Params.N * [0.])
    for i in xrange(Params.N):
        nuFnu[i] = 1.602e-6 * E[i] ** 2 * dNde(result, E[i])
    return E, nuFnu


def MakeSEDInBin(result, Params, Emin_bin, Emax_bin):
    """Compute Spectral energy distribution and corresponding energy
    and return a numpy array"""
    E = np.logspace(np.log10(Emin_bin), np.log10(Emax_bin), Params.N)
    nuFnu = np.array(Params.N * [0.])
    for i in xrange(Params.N):
        nuFnu[i] = 1.602e-6 * E[i] ** 2 * dNde(result, E[i])
    return E, nuFnu


def CountsPlot(Result, Parameter):
    """@todo: document me"""
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    imName = "tmp.fits"
    Result.like.writeCountsSpectra(imName)
    image = pyfits.open(imName)
    #loop on the source names to find the good one
    j = 0
    for ID in image[1].data.names:
        if ID == Parameter.srcname:
            indice = j
        j += 1
    emax = image[3].data.field(0)
    emin = image[3].data.field(1)

    E = array.array('f', ((emax + emin) / 2.))
    err_E = array.array('f', ((-emax + emin) / 2.))

    src = array.array('f', (image[1].data.field(indice)))
    Nbin = len(src)

    obs = array.array('f', (image[1].data.field(0)))
    obs_err = array.array('f', np.sqrt(image[1].data.field(0)))
    summ = 0
    for i in xrange(len(image[1].data.names) - 1):
        summ = summ + image[1].data.field(i + 1)

    other = array.array('f', summ - image[1].data.field(indice))
    summ = array.array('f', summ)

    residual = array.array('f', Nbin * [0.])
    Dres = array.array('f', Nbin * [0.])

    for i in xrange(Nbin):
        try:
            residual[i] = (obs[i] - summ[i]) / summ[i]
        except:
            residual[i] = 0.
        try:
            Dres[i] = (obs_err[i] / summ[i])
        except:
            Dres[i] = 0.

    cmap = ROOT.TCanvas("cmap")
    cmap.SetLogy()
    cmap.SetLogx()
    ghcount = ROOT.TH2F("ghcount", "", 80, 91, 499e3, 100, 0.1, max(obs) * 2)
    ghcount.SetStats(000)
    ghcount.SetXTitle("E (MeV) ")
    ghcount.SetYTitle("Counts / bin")
    ghcount.Draw()
    tgrobs = ROOT.TGraphErrors(Nbin, E, obs, err_E, obs_err)
    tgrobs.SetLineColor(1)
    tgrobs.SetMarkerColor(1)
    tgrobs.SetMarkerStyle(1)
    tgrobs.Draw("P")

    tgrother = ROOT.TGraph(Nbin, E, other)
    tgrother.SetLineWidth(2)
    tgrother.SetLineStyle(2)
    tgrother.Draw("L")

    tgr = ROOT.TGraph(Nbin, E, src)
    tgr.SetLineColor(1)
    tgr.SetLineWidth(2)
    tgr.Draw("L")

    tgrsum = ROOT.TGraph(Nbin, E, summ)
    tgrsum.SetLineStyle(3)
    tgrsum.SetLineWidth(2)
    tgrsum.Draw("L")

    cres = ROOT.TCanvas("cres")
    cres.SetLogx()
    ghres = ROOT.TH2F("ghres", "", 80, 91, 499e3, 100,
                      min(residual) - max(Dres), max(residual) + max(Dres))
    ghres.SetStats(000)
    ghres.SetXTitle("E (MeV) ")
    ghres.SetYTitle("(counts -model)/model")
    ghres.Draw()
    tgres = ROOT.TGraphErrors(Nbin, E, residual, err_E, Dres)
    tgres.Draw("P")

    zero = array.array('f', 2 * [0.])
    Ezero = array.array('f', 2 * [0.])
    Ezero[1] = 1e10
    tg0 = ROOT.TGraph(2, Ezero, zero)
    tg0.SetLineStyle(2)
    tg0.Draw("L")

    # Save the plots in different formats
    cmap.Print(Parameter.PlotName + "_CountsPlot.eps")
    cmap.Print(Parameter.PlotName + "_CountsPlot.C")
    cmap.Print(Parameter.PlotName + "_CountsPlot.png")
    cres.Print(Parameter.PlotName + "_ResPlot.eps")
    cres.Print(Parameter.PlotName + "_ResPlot.C")
    cres.Print(Parameter.PlotName + "_ResPlot.png")
    os.system("rm " + imName)
    image.close()


def Tgraph(like, par):
    """@todo: document me"""
    Res = Result(like, par)
    E, SED = MakeSED(Res, par)
    Err = MakeError(Res, par)
    try:
        CountsPlot(Res, par)
    except:
        pass
    # Save all in ascii file
    # log(E)  log (E**2*dN/dE)   log(E**2*dN/dE_err)  is_dot (0,1) is_upper (0,1)
    save_file = open(par.PlotName + '.dat', 'w')
    save_file.write("log(E)  log (E**2*dN/dE)   log(E**2*dN/dE_err)   \n")
    for i in xrange(par.N):
        save_file.write("%12.4e  %12.4e  %12.4e \n" % (E[i], SED[i], Err[i]))


def PlotTS(Time, TimeErr, TS):
    """@todo: document me"""
    Zero = np.array(len(TS) * [0.])
    gh = ROOT.TH2F("ghts", "", 80, min(Time) - max(TimeErr) * 3,
                   max(Time) + max(TimeErr) * 3, 100, 0, max(TS) * 1.2)
    gh.SetStats(000)
    gh.SetXTitle("Time ")
    gh.SetYTitle("Test Statistic ")
    tgraph = ROOT.TGraphErrors(len(TS), array.array('f', Time),
                               array.array('f', TS), array.array('f', TimeErr),
                               array.array('f', Zero))
    tgraph.SetMarkerColor(1)
    tgraph.SetMarkerStyle(5)
    return gh, tgraph


def PlotNpred(Npred, Flux, FluxErr):
    """@todo: document me"""
    FdF = Flux / (FluxErr + 0.0001)
    NpredErr = np.sqrt(Npred)
    gh = ROOT.TH2F("ghnpred", "", 80, min(Npred / NpredErr) * 0.8, 
                   max(Npred / NpredErr) * 1.2, 100, 
                   min(FdF) * 0.8, max(FdF) * 1.2)
    gh.SetStats(000)
    gh.SetXTitle("Npred/sqrt(Npred) ")
    gh.SetYTitle("Flux/#Delta Flux ")
    tgraph = ROOT.TGraph(len(Npred), array.array('f', Npred / NpredErr),
                         array.array('f', FdF))
    tgraph.SetMarkerColor(1)
    tgraph.SetMarkerStyle(5)
    return gh, tgraph


def PlotLC(Time, TimeErr, Flux, FluxErr):
    """@todo: document me"""
    ArrowSize = (max(Flux) + max(FluxErr) * 1.3 -
                 (min(Flux) - max(FluxErr) * 1.3)) * 0.1
    Arrow = []
    for i in xrange(len(Time)):
        if FluxErr[i] == 0:
            Arrow.append(ROOT.TArrow(Time[i], Flux[i], Time[i],
                                     Flux[i] - ArrowSize, 0.015, "|>"))

    gh = ROOT.TH2F("ghflux", "", 80, min(Time) - max(TimeErr) * 10,
                   max(Time) + max(TimeErr) * 10, 100,
                   min(Flux) - max(FluxErr) * 1.3,
                   max(Flux) + max(FluxErr) * 1.3)
    gh.SetStats(000)
    gh.SetXTitle("Time ")
    gh.SetYTitle("Flux (photon cm^{-2} s^{-1})")

    tgraph = ROOT.TGraphErrors(len(Time), array.array('f', Time),
                               array.array('f', Flux),
                               array.array('f', TimeErr),
                               array.array('f', FluxErr))
    tgraph.SetMarkerColor(1)
    tgraph.SetMarkerStyle(20)
    return gh, tgraph, Arrow


def PlotDataPoints(config):
    """@todo: document me"""
    Arrow = []
    NEbin = int(config['Ebin']['NumEnergyBins'])
    lEmax = np.log10(float(config['energy']['emax']))
    lEmin = np.log10(float(config['energy']['emin']))
    Epoint = array.array('f', NEbin * [0])
    EpointErrp = array.array('f', NEbin * [0])
    EpointErrm = array.array('f', NEbin * [0])
    Fluxpoint = array.array('f', NEbin * [0])
    FluxpointErrp = array.array('f', NEbin * [0])
    FluxpointErrm = array.array('f', NEbin * [0])
    ener = np.logspace(lEmin, lEmax, NEbin + 1)
    for i in xrange(NEbin):
        E = int(pow(10, (np.log10(ener[i + 1]) + np.log10(ener[i])) / 2))
        filename = (config['out'] + '/Ebin/' + config['target']['name'] +
                    "_" + str(E) + ".conf")
        CurConf = get_config(filename)
        try:
            results = utils.ReadResult(CurConf)
        except:
            print "cannot read the Results of energy ", E
            continue
        Epoint[i] = E
        EpointErrm[i] = E - results.get("Emin")
        EpointErrp[i] = results.get("Emax") - E
        try:
            Fluxpoint[i] = 1.6022e-6 * results.get("Ulvalue") * Epoint[i] ** 2
            Arrow.append(ROOT.TArrow(Epoint[i], Fluxpoint[i], Epoint[i],
                                     Fluxpoint[i] * 0.7, 0.02, "|>"))
        except:
            Fluxpoint[i] = 1.6022e-6 * results.get("Prefactor") * Epoint[i] ** 2
            try:
                down = abs(results.get("dPrefactor-"))
                up = results.get("dPrefactor+")
                FluxpointErrp[i] = 1.6022e-6 * up * Epoint[i] ** 2
                FluxpointErrm[i] = 1.6022e-6 * down * Epoint[i] ** 2
            except:
                try:
                    prefactor = results.get("dPrefactor")
                    err = 1.6022e-6 * prefactor * Epoint[i] ** 2
                    FluxpointErrp[i] = err
                    FluxpointErrm[i] = err
                except:
                    pass
        print Epoint[i]
        print Fluxpoint[i]
        print FluxpointErrp[i]
        print FluxpointErrm[i]
    tgpoint = ROOT.TGraphAsymmErrors(NEbin, Epoint, Fluxpoint, EpointErrm,
                                     EpointErrp, FluxpointErrm, FluxpointErrp)
    tgpoint.SetMarkerStyle(20)
    return tgpoint, Arrow


def PlotSED(infile):
    """@todo: document me"""
    config = get_config(infile)
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    root_style.RootStyle()
    FilesName = config['out'] + '/SED_' + config['target']['name']
    AsciiFile = FilesName + '.dat'

    fascii = open(AsciiFile, 'r')
    lines = fascii.readlines()
    ilen = len(lines) - 1
    fascii.close()

    SED = np.array(ilen * [0.])
    E = np.array(ilen * [0.])
    Err = np.array(ilen * [0.])

    for i in xrange(ilen):
        words = string.split(lines[i + 1])
        E[i] = float(words[0])
        SED[i] = float(words[1])
        Err[i] = float(words[2])

    Fluxp = SED + Err
    Fluxm = SED - Err
    ErrorFlux = array.array('f', (2 * ilen + 1) * [0])
    ErrorE = array.array('f', (2 * ilen + 1) * [0])

    for i in xrange(ilen):
        ErrorFlux[i] = Fluxp[i]
        ErrorE[i] = E[i]
    for i in xrange(ilen):
        ErrorFlux[ilen + i] = Fluxm[ilen - i - 1]
        ErrorE[ilen + i] = E[ilen - i - 1]
    ErrorFlux[-1] = Fluxp[0]
    ErrorE[-1] = E[0]

    c_plot = ROOT.TCanvas("c_plot")
    c_plot.SetLogx()
    c_plot.SetLogy()

    ghSED = ROOT.TH2F("ghSED", "", 10000, E[0] * 0.8, E[-1] * 1.5, 100,
                      min(SED[0] - Err[0], SED[-1] - Err[-1]) * 0.2,
                      max(SED[0] + Err[0], SED[-1] + Err[-1]) * 3)
    ghSED.SetStats(000)
    ghSED.SetTitle("Fermi SED")
    ghSED.SetXTitle("E [MeV]")
    ghSED.SetYTitle("E^{2}dN/dE [ erg cm^{-2} s^{-1} ] ")
    ghSED.Draw()

    tgr = ROOT.TGraph(ilen, array.array('f', E), array.array('f', SED))
    tgr.SetLineWidth(2)
    tgr.Draw("L")

    tgerr = ROOT.TGraph(2 * ilen + 1, ErrorE, ErrorFlux)
    tgerr.SetLineColor(2)
    tgerr.Draw("L")

    #Plot points
    NEbin = int(config['Ebin']['NumEnergyBins'])
    if NEbin > 0:
        tgpoint, Arrow = PlotDataPoints(config)
        tgpoint.Draw("pz")

        for i in xrange(len(Arrow)):
            Arrow[i].Draw()

    c_plot.Print(FilesName + '.C')
    c_plot.Print(FilesName + '.eps')
    c_plot.Print(FilesName + '.png')
