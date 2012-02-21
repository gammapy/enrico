import os
from math import sqrt, log10
import numpy as np
import gt_apps
import utils
import pyLikelihood


class Observation(object):
    """@todo: document me"""

    def __init__(self, folder, config, convtyp=-1, tag=""):
        self.config = config
        inttag = "_" + config['file']['tag']
        if tag != "":
            inttag += "_" + tag
        self.srcname = config['target']['name']
        self.ft1 = config['file']['event']
        self.ft2 = config['file']['spacecraft']
        self.xmlfile = config['file']['xml']
        filebase = folder + '/' + self.srcname + inttag
        self.eventfile = filebase + "_Evt.fits"
        self.Cubename = filebase + "_ltCube.fits"
        self.Mapname = filebase + "_ExpMap.fits"
        self.BinnedMapfile = filebase + "_BinnedMap.fits"
        self.cmapfile = filebase + "_CountMap.fits"
        self.ccube = filebase + "_CCUBE.fits"
        self.scrMap = filebase + "_srcMap.fits"
        self.ModelMap = filebase + "_ModelMap.fits"
        self.t1 = float(config['time']['tmin'])
        self.t2 = float(config['time']['tmax'])
        self.Emin = float(config['energy']['emin'])
        self.Emax = float(config['energy']['emax'])
        self.ra = float(config['target']['ra'])
        self.dec = float(config['target']['dec'])
        self.roi = float(config['space']['rad'])
        self.irfs = config['analysis']['irfs']
        if convtyp == 0:
            self.irfs += "::FRONT"
        if convtyp == 1:
            self.irfs += "::BACK"

        self.binsz = self.config['space']['binsz']
        self.npix = int(2 * self.roi / sqrt(2.) / self.binsz)
        self.convtyp = convtyp

    def printSum(self):
        print "Source = ", self.srcname
        print "Center RA = ", self.ra, " degrees"
        print "Center Dec = ", self.dec, " degrees"
        print "Start Time = ", self.t1, "  MET (s)"
        print "Stop Time = ", self.t2, "  MET (s)"
        print "ROI size = ", self.roi, " degrees"
        print "E min = ", self.Emin, " MeV"
        print "E max = ", self.Emax, " MeV"
        print "IRFs = ", self.irfs

    def Gtbin(self):
        tool = gt_apps.evtbin
        tool['evfile'] = self.eventfile
        tool['scfile'] = self.ft2
        tool['outfile'] = self.cmapfile
        tool['algorithm'] = "CMAP"
        tool['nxpix'] = self.npix
        tool['nypix'] = self.npix
        tool['binsz'] = self.binsz
        tool['coordsys'] = 'CEL'
        tool["emin"] = self.Emin
        tool["emax"] = self.Emax
        tool['xref'] = self.ra
        tool['yref'] = self.dec
        tool['axisrot'] = 0
        tool['proj'] = self.config['space']['proj']
        tool.run()

    def _enumbins(self):
        bins_per_decade = self.config['energy']['enumbins_per_decade']
        number_of_decades = log10(self.Emax) - log10(self.Emin)
        enumbins = int(bins_per_decade * number_of_decades)
        return enumbins

    def GtCcube(self):
        tool = gt_apps.evtbin
        tool['evfile'] = self.eventfile
        tool['scfile'] = self.ft2
        tool['outfile'] = self.ccube
        tool['algorithm'] = "CCUBE"
        tool['nxpix'] = self.npix
        tool['nypix'] = self.npix
        tool['binsz'] = self.binsz
        tool['coordsys'] = 'CEL'
        tool['xref'] = self.ra
        tool['yref'] = self.dec
        tool["emin"] = self.Emin
        tool["emax"] = self.Emax
        tool['ebinalg'] = "LOG"
        tool['axisrot'] = 0
        tool['proj'] = "AIT"
        tool["enumbins"] = self._enumbins()
        tool.run()

    def GtBinnedMap(self):
        tool = gt_apps.GtApp('gttool', 'Likelihood')
        tool['infile'] = self.Cubename
        tool['outfile'] = self.BinnedMapfile
        tool['cmap'] = self.ccube
        tool['irfs'] = self.irfs
        tool['emin'] = self.Emin
        tool['emax'] = self.Emax
        tool.run()

    def FirstCut(self):
        tool = gt_apps.filter
        tool['infile'] = self.ft1
        tool['outfile'] = self.eventfile
        tool['ra'] = self.ra
        tool['dec'] = self.dec
        tool['rad'] = self.roi
        tool['emin'] = self.Emin
        tool['emax'] = self.Emax
        tool['tmin'] = self.t1
        tool['tmax'] = self.t2
        tool['zmax'] = self.config['analysis']['zmax']
        tool['evclsmin'] = self.config['analysis']['evclass']
        tool['evclsmax'] = 4
        tool['convtype'] = self.convtyp
        tool.run()

    def MkTime(self):
        tool = gt_apps.filter
        tool['scfile'] = self.ft2
        tool['filter'] = self.config['analysis']['filter']
        tool['roicut'] = 'yes'
        tool['tstart'] = self.t1
        tool['tstop'] = self.t2
        tool['evfile'] = self.eventfile
        tempfile = self.eventfile + ".tmp"
        tool['outfile'] = tempfile
        tool.run()
        os.system("mv %s %s" % (tempfile, self.eventfile))

    def DiffResps(self):
        tool = gt_apps.diffResps
        tool['evfile'] = self.eventfile
        tool['scfile'] = self.ft2
        tool['srcmdl'] = self.xmlfile
        tool['irfs'] = self.irfs
        tool['convert'] = "no"
        tool['evclsmin'] = self.config['analysis']['evclass']
        tool.run()

    def ExpCube(self):
        tool = gt_apps.expCube
        tool['evfile'] = self.eventfile
        tool['scfile'] = self.ft2
        tool['outfile'] = self.Cubename
        tool['dcostheta'] = 0.025
        tool['binsz'] = 1
        tool.run()

    def ExpMap(self):
        tool = gt_apps.expMap
        tool['evfile'] = self.eventfile
        tool['scfile'] = self.ft2
        tool['expcube'] = self.Cubename
        tool['outfile'] = self.Mapname
        tool['irfs'] = self.irfs
        tool['srcrad'] = self.roi + 10
        tool['nenergies'] = self._enumbins() + 1
        tool.run()

    def SrcMap(self):
        tool = gt_apps.srcMaps
        tool['scfile'] = self.ft2
        tool['expcube'] = self.Cubename
        tool['cmap'] = self.ccube
        tool['bexpmap'] = self.BinnedMapfile
        tool['srcmdl'] = self.xmlfile
        tool['irfs'] = self.irfs
        tool['outfile'] = self.scrMap
        tool['emapbnds'] = 'no'
        tool.run()

    def ModelMaps(self, xml):
        tool = gt_apps.model_map
        tool['expcube'] = self.Cubename
        tool['srcmaps'] = self.scrMap
        tool['bexpmap'] = self.BinnedMapfile
        tool['srcmdl'] = xml
        tool["irfs"] = self.irfs
        tool['outfile'] = self.ModelMap
        tool.run()
        utils.SubstracFits(self.cmapfile, self.ModelMap, self.config)

    def GetCovar(self, Fit):
        # @todo: unused variable?
        # ptsrc = pyLikelihood.PointSource_cast(Fit[self.srcname].src)
        par_index_map = {}
        indx = 0
        for src in Fit.sourceNames():
            parNames = pyLikelihood.StringVector()
            Fit[src].src.spectrum().getFreeParamNames(parNames)
            for par in parNames:
                par_index_map["::".join((src, par))] = indx
                indx += 1
            #
            # Build the source-specific covariance matrix.
            #
        if Fit.covariance is None:
            raise RuntimeError("Covariance matrix has not been computed.")
        covar = np.array(Fit.covariance)
        if len(covar) != len(par_index_map):
            raise RuntimeError("Covariance matrix size does not match the " +
                   "number of free parameters.")
        my_covar = []
        srcpars = pyLikelihood.StringVector()
        Fit[self.srcname].src.spectrum().getFreeParamNames(srcpars)
        pars = ["::".join((self.srcname, x)) for x in srcpars]
        for xpar in pars:
            ix = par_index_map[xpar]
            my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
        print "The covariance matrix is :\n", np.array(my_covar)
        print
