import numpy as np
from UnbinnedAnalysis import UnbinnedAnalysis, UnbinnedObs
from BinnedAnalysis import BinnedAnalysis, BinnedObs
import UpperLimits
import IntegralUpperLimit
import utils


class MakeFit(object):
    """@todo: document me"""

    def __init__(self, obs, config):
        self.obs = obs
        self.config = config
        self.task_number = 1
        self.log_like = 0

    def _log(self, task='', description=''):
        print
        print('# ' + '*' * 60)
        if task:
            task = '%10s --- ' % task
        print('# *** %3d %s%s' %
              (self.task_number, task, description))
        print '# ' + '*' * 60
        self.task_number += 1

    def PreparFit(self):
        """@todo: document me"""
        self._log('gtselect', 'Select data from library')
        self.obs.FirstCut()
        self._log('gtmktime', 'Update the GTI and cut data based on ROI')
        self.obs.MkTime()
        self._log('gtdiffrsp', 'Compute Diffuse response')
        self.obs.DiffResps()
        self._log('gtbin', 'Create a count map')
        self.obs.Gtbin()
        self._log('gtltcube', 'Make live time cube')
        self.obs.ExpCube()

        if self.config['analysis']['likelihood'] == 'binned':
            self._log('gtbin', 'Make count map CCUBE')
            self.obs.GtCcube()
            self._log('gtexpcube2', 'Make binned exposure cube')
            self.obs.GtBinnedMap()
            self._log('gtsrcmap', 'Make a source map')
            self.obs.SrcMap()

        if self.config['analysis']['likelihood'] == 'unbinned':
            self._log('gtexpmap', 'Make an exposure map')
            self.obs.ExpMap()

    def CreateFit(self):
        """@todo: document me"""
        if self.config['analysis']['likelihood'] == 'binned':
            Obs = BinnedObs(srcMaps=self.obs.scrMap,
                            expCube=self.obs.Cubename,
                            binnedExpMap=self.obs.BinnedMapfile,
                            irfs=self.obs.irfs)
            Fit = BinnedAnalysis(Obs, self.obs.xmlfile,
                                 optimizer='DRMNGB')

        if self.config['analysis']['likelihood'] == 'unbinned':
            Obs = UnbinnedObs(self.obs.eventfile,
                              self.obs.ft2,
                              expMap=self.obs.Mapname,
                              expCube=self.obs.Cubename,
                              irfs=self.obs.irfs)
            Fit = UnbinnedAnalysis(Obs, self.obs.xmlfile,
                                   optimizer='DRMNGB')
        # @todo: Casting to float should not be necessary
        if float(self.config['Spectrum']['FreeSpectralIndex']) > 0:
            PhIndex = Fit.par_index(self.obs.srcname, 'Index')
            Fit[PhIndex] = -float(self.config['Spectrum']['FreeSpectralIndex'])
            Fit.freeze(PhIndex)

        return Fit

    def PerformFit(self, Fit):
        """@todo: document me"""
        self._log('gtlike', 'Run likelihood analysis')
        try:
            Fit.fit()
        except:
            pass
        Fit.ftol = float(self.config['fitting']['ftol'])
        self.log_like = Fit.fit(covar=True, optimizer=self.config['fitting']['optimizer'])
        Fit.writeXml(self.config['out'] + "/" +
                     self.obs.srcname + "_" +
                     self.config['file']['tag'] + "_out.xml")

        self._log('Results', 'Print results of the fit')
        Result = utils.PrintResult(Fit, self.obs)
        Result['log_like'] = self.log_like
        Result['Emin'] = self.obs.Emin
        Result['Emax'] = self.obs.Emax
        try:
            utils.GetCovar(self.obs.srcname, Fit)
        except:
            pass

        utils.GetFlux(Fit)

        if float(self.config['UpperLimit']['TSlimit']) > Fit.Ts(self.obs.srcname):
            if self.config['UpperLimit']['envelope'] == 'yes':
                self.EnvelopeUL(Fit)
            else:
                Ulval = self.ComputeUL(Fit)
                Result['Ulvalue'] = Ulval
        Result['tmin'] = self.config['time']['tmin']
        Result['tmax'] = self.config['time']['tmax']
        return Result

    def ComputeUL(self, Fit):
        """@todo: document me"""
        self._log('UpperLimit', 'Compute upper Limit')
        print "Assumed index is ", self.config['UpperLimit']['SpectralIndex']
        PhIndex = Fit.par_index(self.obs.srcname, 'Index')
        Fit[PhIndex] = -self.config['UpperLimit']['SpectralIndex']
        Fit.freeze(PhIndex)
        if self.config['UpperLimit']['Method'] == "Profile":
            ul = UpperLimits.UpperLimits(Fit)
            source_ul = ul[self.obs.srcname]
            ul, _ = source_ul.compute(emin=self.obs.Emin,
                                      emax=self.obs.Emax,
                                      delta=2.71 / 2)
            print "Upper limit using Profile method: ", ul
        if self.config['UpperLimit']['Method'] == "Integral":
            ul, _ = IntegralUpperLimit.calc_int(Fit, self.obs.srcname,
                                                verbosity=0)
            print "Upper limit using Integral method: ", ul
        return ul

    def EnvelopeUL(self, Fit):
        """@todo: document me"""
        self._log('EnvelopeUL', 'Compute upper limit envelope')
        PhIndex = Fit.par_index(self.obs.srcname, 'Index')
        Nbp = 20
        Npgraph = 100
        ener = np.logspace(np.log10(self.obs.Emin),
                           np.log10(self.obs.Emax), Npgraph)
        Ulenv = np.array(Npgraph * [0.])

        for i in xrange(Nbp):
            indx = -1.5 - i / (Nbp - 1.)
            Fit[PhIndex] = indx
            Fit.freeze(PhIndex)
            if self.config['UpperLimit']['Method'] == "Profile":
                ul = UpperLimits.UpperLimits(Fit)
                source_ul = ul[self.obs.srcname]
                ul_val, _ = source_ul.compute(emin=self.obs.Emin,
                                              emax=self.obs.Emax,
                                              delta=2.71 / 2)
            if self.config['UpperLimit']['Method'] == "Integral":
                ul_val, _ = IntegralUpperLimit.calc_int(Fit, self.obs.srcname,
                                                        verbosity=0)
            print "Index = ", indx, " UL = ", ul_val
            for j in xrange(Npgraph):
                model_name = Fit.model.srcs[self.obs.srcname].spectrum().genericName()
                if model_name == 'PowerLaw2':
                    newUl = ul_val * (indx + 1) * pow(ener[j], indx + 2) / (pow(self.obs.Emax, indx + 1) - pow(self.obs.Emin, indx + 1))
                elif model_name == 'PowerLaw':
                    IdEScale = utils.getParamIndx(Fit, self.obs.srcname, 'Scale')
                    Escale = Fit[IdEScale].value()
                    newUl = ul_val * pow(ener[j] / Escale, indx + 2)
                Ulenv[j] = max(Ulenv[j], newUl)
        print
        print "Result of the UL envelope"
        for j in xrange(Npgraph):
            print ener[j], " ", Ulenv[j]

    def PlotSED(self, Fit):
        """@todo: document me"""
        self._log('PlotSED', 'Generate SED plot')
        import plotting
        filename = self.config['out'] + '/SED_' + self.obs.srcname
        Par = plotting.Params(self.obs.srcname,
                              Emin=self.obs.Emin,
                              Emax=3e5, extend=False, PlotName=filename)
        plotting.Tgraph(Fit, Par)

    def ModelMap(self, xml):
        """@todo: document me"""
        if self.config['analysis']['likelihood'] == 'binned':
            self._log('gtmodel', 'Make model map')
            self.obs.ModelMaps(xml)
