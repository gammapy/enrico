from UnbinnedAnalysis import UnbinnedAnalysis, UnbinnedObs
from BinnedAnalysis import BinnedAnalysis, BinnedObs
import UpperLimits
import IntegralUpperLimit
from math import log10
import Utility
import numpy

class MakeFit:
    def __init__(self, Observ, config):
        self.Observation = Observ
        self.Configuration = config
        self.OperationNum = 1
        self.logLike = 0

    def _log_message(self, description):
        print
        print '# ' + '*' * 60
        print '# *** ' + str(self.OperationNum) + ' --- ' + description
        print '# ' + '*' * 60
        self.OperationNum += 1

    def PreparFit(self):
        self._log_message('Select data from library - GTSELECT')
        self.Observation.FirstCut()
        self._log_message('Update the GTI and cut data based on ROI - GTMKTIME')
        self.Observation.MkTime()
        self._log_message('Compute Diffuse response - GTDIFFRP')
        self.Observation.DiffResps()
        self._log_message('Create a count map - GTBIN')
        self.Observation.Gtbin()
        self._log_message('Make live time cube - GTLTCUBE')
        self.Observation.ExpCube()

        if self.Configuration['analysis']['likelihood'] == 'binned':
            self._log_message('Make count map CCUBE - GTBIN')
            self.Observation.GtCcube()
            self._log_message('Make Binned exposure map - GTEXPCUBE2')
            self.Observation.GtBinnedMap()
            self._log_message('Make Source Map - GTSRCMAP')
            self.Observation.SrcMap()

        if self.Configuration['analysis']['likelihood'] == 'unbinned':
            self._log_message('Make exposure map - GTEXPMAP')
            self.Observation.ExpMap()

    def CreateFit(self):
        if self.Configuration['analysis']['likelihood'] == 'binned':
            Obs = BinnedObs(srcMaps=self.Observation.scrMap,
                            expCube=self.Observation.Cubename,
                            binnedExpMap=self.Observation.BinnedMapfile,
                            irfs=self.Observation.irfs)
            Fit = BinnedAnalysis(Obs, self.Observation.xmlfile,
                                 optimizer='DRMNGB')

        if self.Configuration['analysis']['likelihood'] == 'unbinned':
            Obs = UnbinnedObs(self.Observation.eventfile,
                              self.Observation.ft2,
                              expMap=self.Observation.Mapname,
                              expCube=self.Observation.Cubename,
                              irfs=self.Observation.irfs)
            Fit = UnbinnedAnalysis(Obs, self.Observation.xmlfile,
                                   optimizer='DRMNGB')
        # @todo: Casting to float should not be necessary
        if float(self.Configuration['Spectrum']['FreeSpectralIndex']) > 0:
            PhIndex = Fit.par_index(self.Observation.srcname, 'Index')
            Fit[PhIndex] = -float(self.Configuration['Spectrum']['FreeSpectralIndex'])
            Fit.freeze(PhIndex)

        return Fit

    def PerformFit(self, Fit):
        self._log_message('Likelihood analysis - GTLIKE')
        try:
            Fit.fit()
        except:
            pass
        Fit.ftol = float(self.Configuration['fitting']['ftol'])
        self.logLike = Fit.fit(covar=True, optimizer=self.Configuration['fitting']['optimizer'])
        Fit.writeXml(self.Configuration['out'] + "/" +
                     self.Observation.srcname + "_" +
                     self.Configuration['file']['tag'] + "_out.xml")

        self._log_message('Result of the fit')
        Result = Utility.PrintResult(Fit, self.Observation)
        Result['logLike'] = self.logLike
        Result['Emin'] = self.Observation.Emin
        Result['Emax'] = self.Observation.Emax
        try:
            Utility.GetCovar(self.Observation.srcname, Fit)
        except:
            pass

        Utility.GetFlux(Fit)

        if float(self.Configuration['UpperLimit']['TSlimit']) > Fit.Ts(self.Observation.srcname):
            if self.Configuration['UpperLimit']['envelope'] == 'yes':
                self.EnvelopeUL(Fit)
            else:
                Ulval = self.ComputeUL(Fit)
                Result['Ulvalue'] = Ulval
        Result['tmin'] = self.Configuration['time']['tmin']
        Result['tmax'] = self.Configuration['time']['tmax']
        return Result

    def ComputeUL(self, Fit):
        self._log_message('Make Upper Limit')
        print "Assumed index is ", self.Configuration['UpperLimit']['SpectralIndex']
        PhIndex = Fit.par_index(self.Observation.srcname, 'Index')
        Fit[PhIndex] = -self.Configuration['UpperLimit']['SpectralIndex']
        Fit.freeze(PhIndex)
        if self.Configuration['UpperLimit']['Method'] == "Profile":
            ul = UpperLimits.UpperLimits(Fit)
            source_ul = ul[self.Observation.srcname]
            ul, _ = source_ul.compute(emin=self.Observation.Emin,
                                      emax=self.Observation.Emax,
                                      delta=2.71 / 2)
            print "Upper limit using Profile method: ", ul
        if self.Configuration['UpperLimit']['Method'] == "Integral":
            ul, _ = IntegralUpperLimit.calc_int(Fit, self.Observation.srcname,
                                                verbosity=0)
            print "Upper limit using Integral method: ", ul
        return ul

    def EnvelopeUL(self, Fit):
        self._log_message('Make Upper Limit Envelope')
        PhIndex = Fit.par_index(self.Observation.srcname, 'Index')
        Nbp = 20
        Npgraph = 100
        ener = numpy.logspace(log10(self.Observation.Emin), log10(self.Observation.Emax), Npgraph)
        Ulenv = numpy.array(Npgraph * [0.])

        for i in xrange(Nbp):
            indx = -1.5 - i / (Nbp - 1.)
            Fit[PhIndex] = indx
            Fit.freeze(PhIndex)
            if self.Configuration['UpperLimit']['Method'] == "Profile":
                ul = UpperLimits.UpperLimits(Fit)
                source_ul = ul[self.Observation.srcname]
                ul_val, _ = source_ul.compute(emin=self.Observation.Emin,
                                              emax=self.Observation.Emax,
                                              delta=2.71 / 2)
            if self.Configuration['UpperLimit']['Method'] == "Integral":
                ul_val, _ = IntegralUpperLimit.calc_int(Fit, self.Observation.srcname,
                                                        verbosity=0)
            print "Index = ", indx, " UL = ", ul_val
            for j in xrange(Npgraph):
                model_name = Fit.model.srcs[self.Observation.srcname].spectrum().genericName()
                if model_name == 'PowerLaw2':
                    newUl = ul_val * (indx + 1) * pow(ener[j], indx + 2) / (pow(self.Observation.Emax, indx + 1) - pow(self.Observation.Emin, indx + 1))
                elif model_name == 'PowerLaw':
                    IdEScale = Utility.getParamIndx(Fit, self.Observation.srcname, 'Scale')
                    Escale = Fit[IdEScale].value()
                    newUl = ul_val * pow(ener[j] / Escale, indx + 2)
                Ulenv[j] = max(Ulenv[j], newUl)
        print
        print "Result of the UL envelope"
        for j in xrange(Npgraph):
            print ener[j], " ", Ulenv[j]

    def PlotSED(self, Fit):
        self._log_message('Generating SED plot')
        import pyPlot
        filename = self.Configuration['out'] + '/SED_' + self.Observation.srcname
        Par = pyPlot.Params(self.Observation.srcname,
                            Emin=self.Observation.Emin,
                            Emax=3e5, extend=False, PlotName=filename)
        pyPlot.Tgraph(Fit, Par)

    def ModelMap(self, xml):
        if self.Configuration['analysis']['likelihood'] == 'binned':
            self._log_message('Make Model Map - GTMODEL')
            self.Observation.ModelMaps(xml)
