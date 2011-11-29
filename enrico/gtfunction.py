# gtfunction.py written by David Sanchez : dsanchez@poly.in2p3.fr
#October 2010

from gt_apps import *
from UnbinnedAnalysis import *
from BinnedAnalysis import *
from GtApp import GtApp
import ROOT
import numpy
import Utility

import os,string
from math import *

class Observation:
    def __init__(self,folder,Configuration,convtyp=-1,tag=""):
	
#	folder = dictionnary.get('folder')+'/'
	self.Configuration = Configuration
	inttag = "_"+Configuration['file']['tag']
	if not(tag==""):
		inttag+="_"+tag

	self.srcname   = Configuration['target']['name']
	self.ft1       = Configuration['file']['event']
	self.ft2       = Configuration['file']['spacecraft']
	self.xmlfile   = Configuration['file']['xml']

	self.eventfile = folder+'/'+self.srcname+inttag+"_Evt.fits"
	self.Cubename  = folder+'/'+self.srcname+inttag+"_ltCube.fits"
	self.Mapname   = folder+'/'+self.srcname+inttag+"_ExpMap.fits"
	self.BinnedMapfile = folder+'/'+self.srcname+inttag+"_BinnedMap.fits"
	self.cmapfile  = folder+'/'+self.srcname+inttag+"_CountMap.fits"
	self.ccube     = folder+'/'+self.srcname+inttag+"_CCUBE.fits"
	self.scrMap    = folder+'/'+self.srcname+inttag+"_srcMap.fits"
	self.ModelMap  = folder+'/'+self.srcname+inttag+"_ModelMap.fits"
	self.t1        = float(Configuration['time']['tmin'])
	self.t2        = float(Configuration['time']['tmax'])
	self.Emin      = float(Configuration['energy']['emin'])
	self.Emax      = float(Configuration['energy']['emax'])
	self.ra        = float(Configuration['target']['ra'])
	self.dec       = float(Configuration['target']['dec'])
	self.roi       = float(Configuration['space']['rad'])
	self.irfs      = Configuration['analysis']['irfs']
	if convtyp==0 :
			self.irfs += "::FRONT"
	if convtyp==1 :
			self.irfs += "::BACK"

	self.binsz     = self.Configuration['space']['binsz']
	self.npix      = int(2*self.roi/sqrt(2.)/self.binsz)
	self.convtyp   = convtyp

    def printSum(self):
	print "Source = ",self.srcname
	print "Center RA = ",self.ra," degrees"
	print "Center Dec = ",self.dec," degrees"
	print "Start Time = ",self.t1,"  MET (s)"
	print "Stop Time = ",self.t2,"  MET (s)"
	print "ROI size = ",self.roi," degrees"
	print "E min = ",self.Emin," MeV"
	print "E max = ",self.Emax," MeV"
	print "IRFs = ",self.irfs


    def Gtbin(self):
#	SizePixel = 2.*self.roi/sqrt(2.)/self.npix
#	np = int(self.npix)
#	npix = int(self.roi*2/SizePixel)

	evtbin['evfile'] = self.eventfile
	evtbin['scfile'] = self.ft2
	evtbin['outfile'] = self.cmapfile
	evtbin['algorithm'] = "CMAP"
	evtbin['nxpix'] = self.npix # 160 #
	evtbin['nypix'] = self.npix #160
	evtbin['binsz'] = self.binsz
	evtbin['coordsys'] = 'CEL'
	evtbin["emin"] = self.Emin
	evtbin["emax"] = self.Emax
	evtbin['xref'] = self.ra 
	evtbin['yref'] = self.dec
	evtbin['axisrot'] = 0
	evtbin['proj'] = self.Configuration['space']['proj'] #"AIT"
	evtbin.run()

    def GtCcube(self):
#	npix    = 200
#	binsize = 2.*self.roi/sqrt(2.)/self.npix

	Nbdecade = log10(self.Emax)-log10(self.Emin)

	evtbin['evfile'] = self.eventfile
	evtbin['scfile'] = self.ft2
	evtbin['outfile'] = self.ccube
	evtbin['algorithm'] = "CCUBE"
	evtbin['nxpix'] = self.npix
	evtbin['nypix'] = self.npix
	evtbin['binsz'] = self.binsz
	evtbin['coordsys'] = 'CEL'
	evtbin['xref'] = self.ra 
	evtbin['yref'] = self.dec
	evtbin["emin"] = self.Emin
	evtbin["emax"] = self.Emax
	evtbin['ebinalg'] = "LOG"
	evtbin['axisrot'] = 0
	evtbin['proj'] = "AIT"
	evtbin["enumbins"] = int(Nbdecade*self.Configuration['energy']['enumbins_per_decade'])
	evtbin.run()


    def GtBinnedMap(self):
#
# Likelihood applications
#
	expcube2 = GtApp('gtexpcube2', 'Likelihood')
	expcube2['infile'] = self.Cubename
	expcube2['outfile'] = self.BinnedMapfile
	expcube2['cmap'] = self.ccube
	expcube2['irfs'] = self.irfs
	expcube2['emin'] = self.Emin
	expcube2['emax'] = self.Emax
	expcube2.run()

    def FirstCut(self):

	filter['infile'] = self.ft1
	filter['outfile'] = self.eventfile
	filter['ra'] = self.ra 
	filter['dec'] = self.dec 
	filter['rad'] = self.roi
	filter['emin'] = self.Emin
	filter['emax'] = self.Emax
	filter['tmin'] = self.t1
	filter['tmax'] = self.t2
	filter['zmax'] = self.Configuration['analysis']['zmax'] 
	filter['evclsmin'] = self.Configuration['analysis']['evclass']
	filter['evclsmax'] = 4
	filter['convtype'] = self.convtyp
	filter.run()

    def MkTime(self):
	maketime['scfile']=self.ft2
	maketime['filter']=self.Configuration['analysis']['filter']
	maketime['roicut']='yes'
	maketime['tstart'] = self.t1
	maketime['tstop'] = self.t2
	maketime['evfile']= self.eventfile
	maketime['outfile']=self.eventfile+".tmp"
	maketime.run()
	os.system("mv "+self.eventfile+".tmp "+self.eventfile)
 
    def DiffResps(self):
	diffResps['evfile']=self.eventfile
	diffResps['scfile']=self.ft2
	diffResps['srcmdl']=self.xmlfile
	diffResps['irfs']=self.irfs
	diffResps['convert']="no"
	diffResps['evclsmin']=self.Configuration['analysis']['evclass']
	diffResps.run()
	print "\ndone"

    def ExpCube(self): 
	"Run gtltcube"
  	expCube['evfile']=self.eventfile
  	expCube['scfile']=self.ft2
  	expCube['outfile'] =self.Cubename
 	expCube['dcostheta']=0.025
  	expCube['binsz']=1
	expCube.run()

    def ExpMap(self):
	"Run gtexpmap"
	expMap['evfile'] = self.eventfile
	expMap['scfile'] = self.ft2
	expMap['expcube'] = self.Cubename
	expMap['outfile'] = self.Mapname
	expMap['irfs']= self.irfs
	expMap['srcrad']=self.roi+10
	expMap['nenergies']=int(self.Configuration['energy']['enumbins_per_decade']*(log10(self.Emax)-log10(self.Emin)))+1
	expMap.run() 

    def SrcMap(self):
	srcMaps['scfile'] = self.ft2
	srcMaps['expcube'] = self.Cubename
	srcMaps['cmap'] = self.ccube
	srcMaps['bexpmap'] = self.BinnedMapfile
	srcMaps['srcmdl']=self.xmlfile
	srcMaps['irfs']= self.irfs
	srcMaps['outfile'] = self.scrMap
	srcMaps['emapbnds']='no'
	srcMaps.run()

    def ModelMaps(self,xml):
	model_map['expcube'] = self.Cubename
	model_map['srcmaps'] = self.scrMap
	model_map['bexpmap'] = self.BinnedMapfile
	model_map['srcmdl'] = xml
	model_map["irfs"]=self.irfs
	model_map['outfile'] = self.ModelMap
	model_map.run()
	Utility.SubstracFits(self.cmapfile,self.ModelMap,self.Configuration)


    def GetCovar(self,Fit):
	ptsrc = pyLike.PointSource_cast(Fit[self.srcname].src)
        par_index_map = {}
        indx = 0
        for src in Fit.sourceNames():
            parNames = pyLike.StringVector()
            Fit[src].src.spectrum().getFreeParamNames(parNames)
            for par in parNames:
               par_index_map["::".join((src, par))] = indx
               indx += 1
       	 	#
       	 	# Build the source-specific covariance matrix.
        	#
        if Fit.covariance is None:
          		raise RuntimeError("Covariance matrix has not been computed.")
        covar = num.array(Fit.covariance)
        if len(covar) != len(par_index_map):
         		raise RuntimeError("Covariance matrix size does not match the " +
                               "number of free parameters.")
        my_covar = []
        srcpars = pyLike.StringVector()
        Fit[self.srcname].src.spectrum().getFreeParamNames(srcpars)
        pars = ["::".join((self.srcname, x)) for x in srcpars]
        for xpar in pars:
           	ix = par_index_map[xpar]
           	my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
        print "The covariance matrix is :\n",numpy.array(my_covar)
	print 


