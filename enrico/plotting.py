import os
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
    err = np.zeros(pars.N)
    j = 0
    for ene in energies:
        arg = pyLikelihood.dArg(ene)
        partials = np.zeros(len(result.srcpars))
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


def MakeFlux(result, params):
    """Compute differential Flux distribution and
    corresponding energy and return a numpy array"""
    E = np.logspace(np.log10(params.Emin), np.log10(params.Emax), params.N)
    Flux = np.zeros(params.N)
    for i in xrange(params.N):
        Flux[i] = dNde(result, E[i])
    return E, Flux


def MakeFluxInBin(result, pars, Emin_bin, Emax_bin):
    """Compute differential Flux distribution and
    corresponding energy and return a numpy array"""
    E = np.logspace(np.log10(Emin_bin), np.log10(Emax_bin), pars.N)
    Flux = np.zeros(pars.N)
    for i in xrange(pars.N):
        Flux[i] = dNde(result, E[i])
    return E, Flux


def MakeSED(result, pars):
    """Compute Spectral energy distribution and corresponding energy
    and return a numpy array"""
    E = np.logspace(np.log10(pars.Emin), np.log10(pars.Emax), pars.N)
    nuFnu = np.zeros(pars.N)
    for i in xrange(pars.N):
        nuFnu[i] = 1.602e-6 * E[i] ** 2 * dNde(result, E[i])
    return E, nuFnu


def MakeSEDInBin(result, pars, Emin_bin, Emax_bin):
    """Compute Spectral energy distribution and corresponding energy
    and return a numpy array"""
    E = np.logspace(np.log10(Emin_bin), np.log10(Emax_bin), pars.N)
    nuFnu = np.zeros(pars.N)
    for i in xrange(pars.N):
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

    E = np.array((emax + emin) / 2.)
    err_E = np.array((-emax + emin) / 2.)

    src = np.array(image[1].data.field(indice))
    Nbin = len(src)

    obs = np.array(image[1].data.field(0))
    obs_err = np.array(np.sqrt(image[1].data.field(0)))
    total = 0
    for i in xrange(len(image[1].data.names) - 1):
        total = total + image[1].data.field(i + 1)

    other = np.array(total - image[1].data.field(indice))
    total = np.array(total)
    residual = np.zeros(Nbin)
    Dres = np.zeros(Nbin)

    for i in xrange(Nbin):
        try:
            residual[i] = (obs[i] - total[i]) / total[i]
        except:
            residual[i] = 0.
        try:
            Dres[i] = (obs_err[i] / total[i])
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

    tgrsum = ROOT.TGraph(Nbin, E, total)
    tgrsum.SetLineStyle(3)
    tgrsum.SetLineWidth(2)
    tgrsum.Draw("L")

    cres = ROOT.TCanvas("cres")
    cres.SetLogx()
    ymin = min(residual) - max(Dres)
    ymax = max(residual) + max(Dres)
    ghres = ROOT.TH2F("ghres", "", 80, 91, 499e3, 100, ymin, ymax)
    ghres.SetStats(000)
    ghres.SetXTitle("E (MeV) ")
    ghres.SetYTitle("(counts -model)/model")
    ghres.Draw()
    tgres = ROOT.TGraphErrors(Nbin, E, residual, err_E, Dres)
    tgres.Draw("P")

    zero = np.zeros(2)
    Ezero = np.array([0, 1e10])
    tg0 = ROOT.TGraph(2, Ezero, zero)
    tg0.SetLineStyle(2)
    tg0.Draw("L")

    # Save the plots in different formats
    filebase = Parameter.PlotName
    cmap.Print(filebase + "_CountsPlot.eps")
    cmap.Print(filebase + "_CountsPlot.C")
    cmap.Print(filebase + "_CountsPlot.png")
    cres.Print(filebase + "_ResPlot.eps")
    cres.Print(filebase + "_ResPlot.C")
    cres.Print(filebase + "_ResPlot.png")
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
    """Scatter plot TS(Time)"""
    Time = np.asarray(Time)
    TimeErr = np.asarray(TimeErr)
    TS = np.asarray(TS)
    zero = np.zeros_like(TS)
    gh = ROOT.TH2F("ghts", "", 80, min(Time) - max(TimeErr) * 3,
                   max(Time) + max(TimeErr) * 3, 100, 0, max(TS) * 1.2)
    gh.SetStats(000)
    gh.SetXTitle("Time")
    gh.SetYTitle("Test Statistic")
    tgraph = ROOT.TGraphErrors(len(TS), Time, TS, TimeErr, zero)
    tgraph.SetMarkerColor(1)
    tgraph.SetMarkerStyle(5)
    return gh, tgraph


def PlotNpred(Npred, Flux, FluxErr):
    """Scatter plot Flux(Npred)"""
    Npred = np.asarray(Npred)
    Flux = np.asarray(Flux)
    FluxErr = np.asarray(FluxErr)
    FdF = Flux / (FluxErr + 0.0001)
    NpredErr = np.sqrt(Npred)
    xmin = min(Npred / NpredErr) * 0.8
    xmax = max(Npred / NpredErr) * 1.2
    ymin, ymax = min(FdF) * 0.8, max(FdF) * 1.2
    gh = ROOT.TH2F("ghnpred", "", 80, xmin, xmax, 100, ymin, ymax)
    gh.SetStats(000)
    gh.SetXTitle("Npred/sqrt(Npred)")
    gh.SetYTitle("Flux/#Delta Flux")
    tgraph = ROOT.TGraph(len(Npred), Npred / NpredErr, FdF)
    tgraph.SetMarkerColor(1)
    tgraph.SetMarkerStyle(5)
    return gh, tgraph


def PlotLC(Time, TimeErr, Flux, FluxErr):
    """Scatter plot Flux(Time)"""
    ArrowSize = (max(Flux) + max(FluxErr) * 1.3 -
                 (min(Flux) - max(FluxErr) * 1.3)) * 0.1
    arrows = []
    for i in xrange(len(Time)):
        if FluxErr[i] == 0:
            arrows.append(ROOT.TArrow(Time[i], Flux[i], Time[i],
                                     Flux[i] - ArrowSize, 0.015, "|>"))
    xmin = min(Time) - max(TimeErr) * 10
    xmax = max(Time) + max(TimeErr) * 10
    ymin = min(Flux) - max(FluxErr) * 1.3
    ymax = max(Flux) + max(FluxErr) * 1.3
    gh = ROOT.TH2F("ghflux", "", 80, xmin, xmax, 100, ymin, ymax)
    gh.SetStats(000)
    gh.SetXTitle("Time ")
    gh.SetYTitle("Flux (photon cm^{-2} s^{-1})")
    tgraph = ROOT.TGraphErrors(len(Time), Time, Flux, TimeErr, FluxErr)
    tgraph.SetMarkerColor(1)
    tgraph.SetMarkerStyle(20)
    return gh, tgraph, arrows


def PlotDataPoints(config):
    """@todo: document me"""
    arrows = []
    NEbin = int(config['Ebin']['NumEnergyBins'])
    lEmax = np.log10(float(config['energy']['emax']))
    lEmin = np.log10(float(config['energy']['emin']))
    Epoint = np.zeros(NEbin)
    EpointErrp = np.zeros(NEbin)
    EpointErrm = np.zeros(NEbin)
    Fluxpoint = np.zeros(NEbin)
    FluxpointErrp = np.zeros(NEbin)
    FluxpointErrm = np.zeros(NEbin)
    ener = np.logspace(lEmin, lEmax, NEbin + 1)
    for i in xrange(NEbin):
        E = int(pow(10, (np.log10(ener[i + 1]) + np.log10(ener[i])) / 2))
        filename = (config['out'] + '/Ebin/' + config['target']['name'] +
                    "_" + str(E) + ".conf")
        CurConf = get_config(filename)
        try:
            print "Reading ",filename
            results = utils.ReadResult(CurConf)
        except:
            print "cannot read the Results of energy ", E
            continue
        Epoint[i] = E
        EpointErrm[i] = E - results.get("Emin")
        EpointErrp[i] = results.get("Emax") - E
        try:
            Fluxpoint[i] = 1.6022e-6 * results.get("Ulvalue") * Epoint[i] ** 2
            arrows.append(ROOT.TArrow(Epoint[i], Fluxpoint[i], Epoint[i],
                                     Fluxpoint[i] * 0.7, 0.02, "|>"))
        except:
            Fluxpoint[i] = 1.6022e-6 * results.get("Prefactor") * Epoint[i] ** 2
            try:
                down = abs(results.get("dPrefactor-"))
                up = results.get("dPrefactor+")
                if down==0 or  up ==0 :
                  raise RuntimeError("error")
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
        print "Energy = ",Epoint[i]
        print "E**2. dN/dE = ",Fluxpoint[i]," + ",FluxpointErrp[i]," - ",FluxpointErrm[i]

    tgpoint = ROOT.TGraphAsymmErrors(NEbin, Epoint, Fluxpoint, EpointErrm,
                                     EpointErrp, FluxpointErrm, FluxpointErrp)
    tgpoint.SetMarkerStyle(20)
    return tgpoint, arrows


def PlotSED(infile):
    """@todo: document me"""
    config = get_config(infile)
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    root_style.RootStyle()
    filebase = config['out'] + '/SED_' + config['target']['name']
    lines = open(filebase + '.dat', 'r').readlines()
    ilen = len(lines) - 1

    SED = np.zeros(ilen)
    E = np.zeros(ilen)
    Err = np.zeros(ilen)

    for i in xrange(ilen):
        words = lines[i + 1].split()
        E[i] = float(words[0])
        SED[i] = float(words[1])
        Err[i] = float(words[2])

    Fluxp = SED*np.exp(Err/SED)#SED + Err
    Fluxm =  SED*np.exp(-Err/SED)#SED - Err
    ErrorFlux = np.zeros(2 * ilen + 1)
    ErrorE = np.zeros(2 * ilen + 1)

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

    xmin, xmax = E[0] * 0.8, E[-1] * 1.5
    ymin = min(SED[0] - Err[0], SED[-1] - Err[-1]) * 0.2
    ymax = max(SED[0] + Err[0], SED[-1] + Err[-1]) * 3
    ghSED = ROOT.TH2F("ghSED", "", 10000, xmin, xmax, 100, ymin, ymax)
    ghSED.SetStats(000)
    ghSED.SetTitle("Fermi SED")
    ghSED.SetXTitle("E [MeV]")
    ghSED.SetYTitle("E^{2}dN/dE [ erg cm^{-2} s^{-1} ] ")
    ghSED.Draw()

    tgr = ROOT.TGraph(ilen, E, SED)
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

#TODO add a writeTOASCII

    c_plot.Print(filebase + '.C')
    c_plot.Print(filebase + '.eps')
    c_plot.Print(filebase + '.png')
