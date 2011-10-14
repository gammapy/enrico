# fitfunction.py written by David Sanchez : dsanchez@poly.in2p3.fr
#September 2011

from UnbinnedAnalysis import *
from BinnedAnalysis import *
import UpperLimits
import IntegralUpperLimit
from math import *
import string
import Utility
import gtfunction
import numpy

class MakeFit:
    def __init__(self,Observ,Configuration):
	self.Observation = Observ
	self.Configuration = Configuration
	self.OperationNum = 1

    def PreparFit(self):
	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Select data from library - GTSELECT'
	print '# *********************************************************************'
	self.Observation.FirstCut()
	self.OperationNum+=1

	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Update the GTI and cut data based on ROI - GTMKTIME'
	print '# *********************************************************************'
	self.Observation.MkTime()
	self.OperationNum+=1

	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Compute Diffuse response - GTDIFFRP'
	print '# *********************************************************************'
	self.Observation.DiffResps()
	self.OperationNum+=1

	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Create a count map - GTBIN'
	print '# *********************************************************************'
	self.Observation.Gtbin()
	self.OperationNum+=1

	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Make live time cube - GTLTCUBE'
	print '# *********************************************************************'
	self.Observation.ExpCube()
	self.OperationNum+=1


	
	if self.Configuration['analysis']['likelihood'] == 'binned':

		print
		print '# *********************************************************************'
		print '# * '+str(self.OperationNum)+' - Make count map CCUBE - GTBIN'
		print '# *********************************************************************'
		self.Observation.GtCcube()
		self.OperationNum+=1

		print
		print '# *********************************************************************'
		print '# * '+str(self.OperationNum)+' - Make Binned exposure map - GTEXPCUBE2'
		print '# *********************************************************************'
		self.Observation.GtBinnedMap()
		self.OperationNum+=1

		print
		print '# *********************************************************************'
		print '# * '+str(self.OperationNum)+' - Make Source Map - GTSRCMAP'
		print '# *********************************************************************'
		self.Observation.SrcMap()
		self.OperationNum+=1

	if self.Configuration['analysis']['likelihood'] == 'unbinned':

		print
		print '# *********************************************************************'
		print '# * '+str(self.OperationNum)+' - Make exposure map - GTEXPMAP'
		print '# *********************************************************************'
		self.Observation.ExpMap()
		self.OperationNum+=1

    def CreateFit(self) :
	if self.Configuration['analysis']['likelihood'] == 'binned':
		Obs = BinnedObs(srcMaps=self.Observation.scrMap,expCube=self.Observation.Cubename,binnedExpMap=self.Observation.BinnedMapfile,irfs=self.Observation.irfs)

		Fit = BinnedAnalysis(Obs,self.Observation.xmlfile,optimizer='DRMNGB')

	if self.Configuration['analysis']['likelihood'] == 'unbinned':

		Obs=UnbinnedObs(self.Observation.eventfile,self.Observation.ft2,expMap=self.Observation.Mapname,expCube=self.Observation.Cubename,irfs=self.Observation.irfs)

		Fit = UnbinnedAnalysis(Obs,self.Observation.xmlfile,optimizer='DRMNGB')

	if float(self.Configuration['enricobehavior']['FreeSpectralIndex']) >0 :
		PhIndex = Fit.par_index(self.Observation.srcname, 'Index')
		Fit[PhIndex] = -float(self.Configuration['enricobehavior']['FreeSpectralIndex'])
		Fit.freeze(PhIndex)

	return Fit

    def PerformFit(self,Fit):

	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Likelihood analysis - GTLIKE'
	print '# *********************************************************************'
	self.OperationNum+=1
	try :
		Fit.fit()
	except:
		pass
	Fit.ftol=float(self.Configuration['fitting']['ftol'])
	Fit.fit(covar=True,optimizer=self.Configuration['fitting']['optimizer'])
	Fit.writeXml(self.Configuration['out']+"/"+self.Observation.srcname+"_out.xml")


	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Result of the fit'
	print '# *********************************************************************'
	self.OperationNum+=1
	try :
		Utility.GetCovar(self.Observation.srcname,Fit)
	except :
		pass
	Utility.PrintResult(Fit,self.Observation)
	Utility.GetFlux(Fit)
	if float(self.Configuration['UpperLimit']['TSlimit'])>Fit.Ts(self.Observation.srcname) :
		if self.Configuration['UpperLimit']['envelope'] =='yes' :
			self.EnvelopeUL(Fit)
		else :
			self.ComputeUL(Fit)


    def ComputeUL(self,Fit) :
	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Make Upper Limit'
	print '# *********************************************************************'
	self.OperationNum+=1


	print "Assumed index is ",self.Configuration['UpperLimit']['SpectralIndex']
	PhIndex = Fit.par_index(self.Observation.srcname, 'Index')
	Fit[PhIndex] = -self.Configuration['UpperLimit']['SpectralIndex']
	Fit.freeze(PhIndex)

	if self.Configuration['UpperLimit']['Method'] == "Profile" :
		ul = UpperLimits.UpperLimits(Fit)
		ul_prof, par_prof = ul[self.Observation.srcname].compute(emin=self.Observation.Emin, emax=self.Observation.Emax, delta=2.71/2)

		print "Upper limit using Profile method: ", ul_prof

	if self.Configuration['UpperLimit']['Method'] == "Integral" :
		ul, results = IntegralUpperLimit.calc_int(Fit, self.Observation.srcname, verbosity=0)

		print "Upper limit using Integral method: ", ul

    def EnvelopeUL(self,Fit) :
	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Make Upper Limit Envelope'
	print '# *********************************************************************'
	self.OperationNum+=1

	PhIndex = Fit.par_index(self.Observation.srcname, 'Index')
	Nbp = 20
	Npgraph = 100
	ener = numpy.logspace(log10(self.Observation.Emin),log10(self.Observation.Emax),Npgraph)
	Ulenv = numpy.array(Npgraph*[0.]) 

	for i in xrange(Nbp):
		indx = -1.5-i/(Nbp-1.)
		Fit[PhIndex] = indx
		Fit.freeze(PhIndex)
		if self.Configuration['UpperLimit']['Method'] == "Profile" :
			ul = UpperLimits.UpperLimits(Fit)
			ul_val, par_prof = ul[self.Observation.srcname].compute(emin=self.Observation.Emin, emax=self.Observation.Emax, delta=2.71/2)

		if self.Configuration['UpperLimit']['Method'] == "Integral" :
			ul_val, results = IntegralUpperLimit.calc_int(Fit, self.Observation.srcname, verbosity=0)

		print "Index = ",indx," UL = ",ul_val

		for j in xrange(Npgraph):
			if Fit.model.srcs[self.Observation.srcname].spectrum().genericName() == 'PowerLaw2':
				newUl = ul_val*(indx+1)*pow(ener[j],indx+2)/(pow(self.Observation.Emax,indx+1)-pow(self.Observation.Emin,indx+1))

			if Fit.model.srcs[self.Observation.srcname].spectrum().genericName() == 'PowerLaw':
				IdEScale  = getParamIndx(Fit,self.Observation.srcname,'Scale')
				Escale  = Fit[IdEScale].value()
				newUl = ul_val*pow(ener[j]/Escale,indx+2)

			Ulenv[j] = max(Ulenv[j],newUl)

	print 
	print "Resultat of the UL envelope"
	for j in xrange(Npgraph):
		print ener[j]," ",Ulenv[j]

	

    def PlotSED(self,Fit) :
	print
	print '# *********************************************************************'
	print '# * '+str(self.OperationNum)+' - Generating SED plot'
	print '# *********************************************************************'
	self.OperationNum+=1
	import pyPlot as P
#	try :
	Par=P.Params(self.Observation.srcname,Emin=self.Observation.Emin,Emax=1e7,extend=True,PlotName=self.Configuration['out']+'/SED_'+self.Observation.srcname)
	P.Tgraph(Fit,Par)
#		print "... Done"
#	except :
#		print "Error in the generation of the SED"

    def ModelMap(self,xml) :
	if self.Configuration['analysis']['likelihood'] == 'binned':
		print
		print '# *********************************************************************'
		print '# * '+str(self.OperationNum)+' - Make Model Map - GTMODEL'
		print '# *********************************************************************'
		self.OperationNum+=1
		self.Observation.ModelMaps(xml)

