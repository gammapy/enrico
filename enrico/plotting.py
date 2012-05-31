import os
import ROOT
import numpy as np
import pyfits
import root_style
import pyLikelihood
from config import get_config
import utils


class Params:
    """Collection of Plotting parameters like Energy bounds,
    colors, file name, etc...."""
    def __init__(self, srcname, Emin=100, Emax=3e5,
                 PlotName="LAT_SED", LineColor=2,
                 PointColor = 1, N = 500):
        self.Emin = Emin #Energy bounds
        self.Emax = Emax
        self.N = N #Number of points for the TGraph
        self.srcname = srcname # Source of interest
        self.PlotName = PlotName #file name
        #color options
        self.LineColor = LineColor
        self.PointColor = PointColor

class Result:
    """Helper class to get the results from a (Un)BinnedAnalysis object
    and compute the SED and errors"""
    def __init__(self, Fit, pars):
        self.Fit = Fit
        self.Model = Fit[pars.srcname].funcs['Spectrum'].genericName()
        self.ptsrc = pyLikelihood.PointSource_cast(Fit[pars.srcname].src)
        self.covar = np.array(utils.GetCovar(pars.srcname, self.Fit, False))
        self.srcpars = pyLikelihood.StringVector()
        Fit[pars.srcname].src.spectrum().getFreeParamNames(self.srcpars)

    def _DumpSED(self,par):
        """Save the energy, E2.dN/dE, and corresponding  error in an ascii file
        The count and residuals plot vs E is also made"""
        E, SED = self.MakeSED(par)
        Err = self.MakeSEDError(par)

        print
        for i in xrange(par.N):
          if Err[i] == min(Err):
            print "Decorrelation energy : %4.2e MeV"%E[i]
            print "SED at the Decorrelation energy : %2.2e +/-  %2.2e erg/cm2/s" %(SED[i],Err[i])
        print

        try:
            self.CountsPlot(par)
        except:
            pass
        # Save all in ascii file
        # log(E)  log (E**2*dN/dE)   log(E**2*dN/dE_err)  is_dot (0,1) is_upper (0,1)
        save_file = open(par.PlotName + '.dat', 'w')
        save_file.write("log(E)  log (E**2*dN/dE)   Error on log(E**2*dN/dE)   \n")
        for i in xrange(par.N):
            save_file.write("%12.4e  %12.4e  %12.4e \n" % (E[i], SED[i], Err[i]))
        save_file.close()

    def MakeFlux(self, params):
        """Compute differential Flux distribution and
        corresponding energy and return a numpy array"""
        E = np.logspace(np.log10(params.Emin), np.log10(params.Emax), params.N)
        Flux = np.zeros(params.N)
        for i in xrange(params.N):
            Flux[i] = self.dNde(E[i])
        return E, Flux

    def MakeSED(self, pars):
        """Compute Spectral energy distribution and corresponding energy
        and return a numpy array"""
        E = np.logspace(np.log10(pars.Emin), np.log10(pars.Emax), pars.N)
        nuFnu = np.zeros(pars.N)
        for i in xrange(pars.N):
            nuFnu[i] = 1.602e-6 * E[i] ** 2 * self.dNde(E[i]) #Mev to Ergs
        return E, nuFnu

    def MakeSEDError(self, pars):
        """@todo: document me"""
        estep = np.log(pars.Emax / pars.Emin) / (pars.N - 1)
        energies = pars.Emin * np.exp(estep * np.arange(np.float(pars.N)))
        err = np.zeros(pars.N)
        j = 0
        for ene in energies:
            arg = pyLikelihood.dArg(ene)
            partials = np.zeros(len(self.srcpars))
            for i in xrange(len(self.srcpars)):
                x = self.srcpars[i]
                partials[i] = self.ptsrc.spectrum().derivByParam(arg, x)
            err[j] = np.sqrt(np.dot(partials, np.dot(self.covar, partials)))
            j += 1

        return 1.602e-6 * energies ** 2 * err #Mev to Ergs

    def dNde(self, energy):
        arg = pyLikelihood.dArg(energy)
        return self.ptsrc.spectrum()(arg)

    def CountsPlot(self, Parameter):
        """@todo: document me"""
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        imName = "tmp.fits"
        self.Fit.writeCountsSpectra(imName)
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
                Dres[i] = (obs_err[i] / total[i])
            except:
                residual[i] = 0.
                Dres[i] = 0.

        cplot = ROOT.TCanvas("Count Plot")
        cplot.SetLogy()
        cplot.SetLogx()
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

        cres = ROOT.TCanvas("Residuals plot")
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
        cplot.Print(filebase + "_CountsPlot.eps")
        cplot.Print(filebase + "_CountsPlot.C")
        cploy.Print(filebase + "_CountsPlot.png")
        cres.Print(filebase + "_ResPlot.eps")
        cres.Print(filebase + "_ResPlot.C")
        cres.Print(filebase + "_ResPlot.png")
        os.system("rm " + imName)
        image.close()

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
    """Collect the data points/UL and generate a TGraph for the points
    and a list of TArrow for the UL. All is SED format"""

    ErgsToMeV = 1.6022e-6

    #Preparation + declaration of arrays
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

    for i in xrange(NEbin):#Loop over the energy bins
        E = int(pow(10, (np.log10(ener[i + 1]) + np.log10(ener[i])) / 2))
        filename = (config['out'] + '/Ebin'+str(NEbin)+'/' + config['target']['name'] +
                    "_" + str(i) + ".conf")
        try:#read the config file of each data points
            CurConf = get_config(filename)
            print "Reading ",filename
            results = utils.ReadResult(CurConf)
        except:
            print "cannot read the Results of energy ", E
            continue
        #fill the energy arrays
        Epoint[i] = E
        EpointErrm[i] = E - results.get("Emin")
        EpointErrp[i] = results.get("Emax") - E

        #Compute the flux or the UL (in SED format)
        if results.has_key('Ulvalue'):
            PrefUl = utils.Prefactor(results.get("Ulvalue"),results.get("Index"),
                                    results.get("Emin"),results.get("Emax"),E)
            Fluxpoint[i] = ErgsToMeV * PrefUl * Epoint[i] ** 2
            arrows.append(ROOT.TArrow(Epoint[i], Fluxpoint[i], Epoint[i],
                                     Fluxpoint[i] * 0.5, 0.02, "|>"))
        else : #Not an UL : compute points + errors
            Fluxpoint[i] = ErgsToMeV * results.get("Prefactor") * Epoint[i] ** 2
            try:
                down = abs(results.get("dPrefactor-"))
                up = results.get("dPrefactor+")
                if down==0 or  up ==0 :
                  raise RuntimeError("cannot get Error value")
                FluxpointErrp[i] = ErgsToMeV * up * Epoint[i] ** 2
                FluxpointErrm[i] = ErgsToMeV * down * Epoint[i] ** 2
            except:
                try:
                    prefactor = results.get("dPrefactor")
                    err = ErgsToMeV * prefactor * Epoint[i] ** 2
                    FluxpointErrp[i] = err
                    FluxpointErrm[i] = err
                except:
                    pass
        print "Energy = ",Epoint[i]
        print "E**2. dN/dE = ",Fluxpoint[i]," + ",FluxpointErrp[i]," - ",FluxpointErrm[i]

    #create a TGraph for the points
    tgpoint = ROOT.TGraphAsymmErrors(NEbin, Epoint, Fluxpoint, EpointErrm,
                                     EpointErrp, FluxpointErrm, FluxpointErrp)
    tgpoint.SetMarkerStyle(20)
    return tgpoint, arrows


def PlotSED(infile,pars):
    """plot a nice SED with a butterfly and points"""
    config = get_config(infile)
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    root_style.RootStyle()

    # Read the ascii file where the butterfly is stored
    filebase = utils._SpecFileName(config)

    lines = open(filebase + '.dat', 'r').readlines()
    SED = []
    E = []
    Err = []

    for i in xrange(len(lines) - 1):
        words = lines[i + 1].split()
        if float(words[0])<pars.Emax :
            E.append(float(words[0]))
            SED.append(float(words[1]))
            Err.append(float(words[2]))
    ilen = len(SED)

    #From dN/dE to SED
    Fluxp = np.array(SED)*np.exp(np.array(Err)/np.array(SED))
    Fluxm =  np.array(SED)*np.exp(-np.array(Err)/np.array(SED))
    ErrorFlux = np.zeros(2 * ilen + 1)
    ErrorE = np.zeros(2 * ilen + 1)

    #Compute the butterfly and close it
    for i in xrange(ilen):
        ErrorFlux[i] = Fluxp[i]
        ErrorE[i] = E[i]
    for i in xrange(ilen):
        ErrorFlux[ilen + i] = Fluxm[ilen - i - 1]
        ErrorE[ilen + i] = E[ilen - i - 1]
    ErrorFlux[-1] = Fluxp[0]
    ErrorE[-1] = E[0]

    #Actually make the plot
    c_plot = ROOT.TCanvas(pars.PlotName)
    c_plot.SetLogx()
    c_plot.SetLogy()

    xmin, xmax = E[0] * 0.8, E[-1] * 1.5
    ymin = min(SED[0] - Err[0], SED[-1] - Err[-1]) * 0.2
    ymax = max(SED[0] + Err[0], SED[-1] + Err[-1]) * 3
    ghSED = ROOT.TH2F("ghSED", "", 10000, xmin, xmax, 100, ymin, ymax)
    ghSED.SetStats(000)
    ghSED.SetTitle(pars.PlotName)
    ghSED.SetXTitle("E [MeV]")
    ghSED.SetYTitle("E^{2}dN/dE [ erg cm^{-2} s^{-1} ] ")
    ghSED.Draw()

    tgr = ROOT.TGraph(ilen, np.array(E), np.array(SED))
    tgr.SetLineWidth(2)
    tgr.SetLineColor(pars.LineColor)
    tgr.Draw("L")

    tgerr = ROOT.TGraph(2 * ilen + 1, ErrorE, ErrorFlux)
    tgerr.SetLineColor(pars.LineColor)
    tgerr.Draw("L")

    #Plot points
    NEbin = int(config['Ebin']['NumEnergyBins'])
    if NEbin > 0:
        tgpoint, Arrow = PlotDataPoints(config) #collect data points
        tgpoint.SetLineColor(pars.PointColor)
        tgpoint.SetMarkerColor(pars.PointColor)
        tgpoint.Draw("pz")
        for i in xrange(len(Arrow)):
            Arrow[i].SetLineColor(pars.PointColor)
            Arrow[i].SetFillColor(pars.PointColor)
            Arrow[i].Draw()

    #save the canvas
    c_plot.Print(filebase + '.C')
    c_plot.Print(filebase + '.eps')
    c_plot.Print(filebase + '.png')

def PlotUL(pars,config,ULFlux,Index):

    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    root_style.RootStyle()

    #Compute the SED
    E = np.logspace(np.log10(pars.Emin), np.log10(pars.Emax), pars.N)
    SED = 1.602e-6 * E ** 2 * (-Index+1)*ULFlux* np.power(E,-Index)/(np.power(pars.Emax,-Index+1)-np.power(pars.Emin,-Index+1))

    #Actually make the plot
    c_plot = ROOT.TCanvas(pars.PlotName)
    c_plot.SetLogx()
    c_plot.SetLogy()

    xmin, xmax = E[0] * 0.7, E[-1] * 1.6
    ymin = min(SED[0], SED[-1]) * 0.15
    ymax = max(SED[0], SED[-1]) * 3
    ghSED = ROOT.TH2F("ghSED", "", 10000, xmin, xmax, 100, ymin, ymax)
    ghSED.SetStats(000)
    ghSED.SetTitle(pars.PlotName)
    ghSED.SetXTitle("E [MeV]")
    ghSED.SetYTitle("E^{2}dN/dE [ erg cm^{-2} s^{-1} ] ")
    ghSED.Draw()

    tgr = ROOT.TGraph(pars.N, np.array(E), np.array(SED))
    tgr.Draw("L")

    Ar_1=ROOT.TArrow(E[0],SED[0]*0.2,E[0],SED[0],0.02)
    Ar_2=ROOT.TArrow(E[-1],SED[-1]*0.2,E[-1],SED[-1],0.02)
    Ar_2.Draw("<|")
    Ar_1.Draw("<|")

    #save the canvas
    filebase = utils._SpecFileName(config)
    c_plot.Print(filebase + '.C')
    c_plot.Print(filebase + '.eps')
    c_plot.Print(filebase + '.png')

