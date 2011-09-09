#Define some classes

#this class takes the very basic inputs and gathers the necessary information that one isn't likley to change
#also define class functions to easily access the other info gathering functions
class likeInput:
	#the arguments are: 
	#LIKE (an UnbinnedAnalysis object),
	#model (xml model for fitting the energy bands), and
	#SrcName (name of your source of interest)
	#nbins (integer number of bins)
	#model and SrcNames are assumed to be strings
	
	def __init__(self,LIKE,SrcName,model='',nbins=20,phCorr=1.0):
		if(model!=''): #if given a value for model, make sure file exists
			if(not fileCheck(model)):
				return
		if(model==''): #if no value given for model, set model to that in LIKE
			model=LIKE.srcModel
		self.ubAn=LIKE
		self.ft1=LIKE.observation.eventFiles[0]
		self.ft2=LIKE.observation.scFiles[0]
		self.expMap=LIKE.observation.expMap
		self.expCube=LIKE.observation.expCube
		self.IRFs=LIKE.observation.irfs
		self.source=SrcName
		self.bandModel=model
		self.NBins=nbins
		self.ft1ROI=self.ubAn.observation.roiCuts().roiCone()
		self.ft1EBounds=self.ubAn.observation.roiCuts().getEnergyCuts()
		self.phCorr=phCorr #allow a phase correction factor for pulsar analyses
	
	#print out the parameters
	def Print(self):
		print self.ubAn
		print 'Model for energy bands:', self.bandModel
		print 'Source of Interest:', self.source
	
	#direct acces to maxEnergy function
	def srcMaxE(self):
		file=pyfits.open(self.ft1)
		MAX=maxEnergy(file,self)
		file.close()
		return MAX
	
	#makes bins for plotting taking max energy into account, makes NBins bands for full range and then accesses plotEbounds function
	def plotBins(self,MaxE=0,evclsmin=3,evclsmax=4,evclass=None):
		if(MaxE==0):
			MaxE=self.srcMaxE()
		ebounds=log_array(self.NBins+1,self.ft1EBounds[0],self.ft1EBounds[1])
		self.bins=plotEbounds(ebounds,MaxE)
		self.nbins=len(self.bins[0])
		self.makeFiles(evclsmin,evclsmax,evclass)
		self.custBins=0
	
	#function for those who want to define their own plot bins in a manner different from the plotBins function
	def customBins(self,MIN,MAX,evclsmin=3,evclsmax=4,evclass=None):#here MIN and MAX are the bin low and high edges (respectively) as you've defined them
		self.bins=[MIN,MAX]
		self.nbins=len(MIN)
		self.makeFiles(evclsmin,evclsmax,evclass)
		self.custBins=1
	
	#runs both makeFT1Bands and makeExpMaps
	def makeFiles(self,evclsmin,evclsmax,evclass):
		makeFT1Bands(self,evclsmin,evclsmax,evclass)
		makeExpMaps(self)
	
	#allows one to access the full energy range without redoing the fit if it has already been done
	def getfullFit(self,Emin=0,Emax=0,expCorrect=False,writeXML=False):
		if(Emin==0):
			try:
				Emin=self.bins[0][0]
			except:
				self.plotBins()
				Emin=self.bins[0][0]
		if(Emax==0):
			try:
				Emax=self.bins[1][-1]
			except:
				self.plotBins()
				Emax=self.bins[1][-1]
		self.FitBounds=[Emin,Emax]
		self.Fit=getFit(self,Emin,Emax,51,1,expCorrect=expCorrect,wx=writeXML)
	
	#run the full fit on LIKE, objective to get model and bowtie plot info
	def fullFit(self,ftol=1e-3,CoVar=False,toltype='ABS',Emin=0,Emax=0,writeXML=False): #must say CoVar=True if one wants bowtie info
		if(toltype!='ABS' and toltype!='REL'):
			print 'Invalid toltype parameter, valid choices are "REL" for relative tolerance and "ABS" for absolute tolerance.'
			return
		if(toltype=='ABS'):
			self.ubAn.setFitTolType(1)
		else:
			self.ubAn.setFitTolType(0)
		if(Emin==0):
			try:
				Emin=self.bins[0][0]
			except:
				self.plotBins()
				Emin=self.bins[0][0]
		if(Emax==0):
			try:
				Emax=self.bins[1][-1]
			except:
				self.plotBins()
				Emax=self.bins[1][-1]
		self.FitBounds=[Emin,Emax]
		self.Fit=fullLike(self,CoVar,ftol,Emin,Emax,wx=writeXML) #returns bowtie info and best fit model (in that order)

class likeSED:
	def __init__(self,Input):
		#the arguements are:
		#Input (a likeInput object)
		self.likeIn=Input
		try:
			self.Fit=Input.Fit[1]
			self.BT=Input.Fit[0]
			self.FitBounds=Input.FitBounds
		except:
			Input.getfullFit()
			self.Fit=Input.Fit[1]
			self.BT=Input.Fit[0]
			self.FitBounds=Input.FitBounds
	
	#get center energies based on best fit model, basically just gives access to specWeightedCenters function.
	def getECent(self):
		self.centers=specWeightedCenters(self.likeIn)
	
	#define custom center energies (i.e. if you want different or no weighting)
	def customECent(self,Cent):
		self.centers=Cent
	
	#do likelihood fits in each band
	def fitBands(self,ftol=1e-3,tslim=25,toltype='ABS',opt='NewMinuit',lastbinUL=False,rescaleAll=False,writeXML=False):
		if(toltype!='ABS' and toltype!='REL'):
			print 'Invalid toltype parameter, valid choices are "REL" for relative tolerance and "ABS" for absolute tolerance.'
			return
		if(toltype=='ABS'):
			ttype=1
		else:
			ttype=0
		self.data=runLike(self.likeIn,self.centers,ftol,tslim,ttype,opt,rescaleAll,lastbinUL,writeXML)
	
	#plot the results and create and save the histograms, .eps, and fits results files.
	def Plot(self,plot=True):
		try:
			len(self.data)
		except:
			self.fitBands
		plotSED(self,plot)
	
	#make a residuals plot
	def residuals(self,plot=True):
		try:
			len(self.data)
		except:
			self.fitBands()
		try:
			self.likeIn.ubAn()
			residualsPlot(self,plot)
		except:
			print 'Full energy range UnbinnedAnalysis object has been deleted (probably with likeInput.delFullLike() function), can not make residuals plot.'
			return

from GtApp import GtApp
import numpy as num
import pyfits
import os
import pyIrfLoader
from UnbinnedAnalysis import *
from ROOT import TCanvas, TH1F, gStyle, TLegend, TFile, TArrow, TLine, TPostScript, TF1, SetOwnership, TGraph, TGraphAsymmErrors
from UpperLimits import UpperLimits
import pyLikelihood as pyLike
from math import sin,cos,log10,pi,log,ceil,floor
print 'This is likeSED version 12.1, modified to handle Pass 7 selections.'
#this function creates an array with npts-1 logarithmically spaced bins
#borrowed from macro Jim sent to do profile likelihood
def log_array(npts, xmin, xmax):
	xstep = num.log(xmax/xmin)/(npts - 1)
	return xmin*num.exp(num.arange(npts, dtype=num.float)*xstep)

#this function creates the different energy band ft1 files from the original ft1
def makeFT1Bands(likeIn,evmin,evmax,evclass):
	gtselect=GtApp('gtselect') #note, checks for existence of ft1 bands, if found won't remake 
	Src=likeIn.source
	NBins=likeIn.NBins
	ft1=likeIn.ft1
	ebins=likeIn.bins
	for i in range(0,likeIn.nbins):
		out='%s_%ibins_band%i.fits' %(Src.replace(' ','_'),NBins,i) #ft1 has name of source and number of bins to be specific about what was done
		if(not os.access(out,os.F_OK)):
			if evclass==None:
				gtselect.run(infile=ft1,outfile=out,ra=0,dec=0,rad=180,tmin=0,tmax=0,emin=ebins[0][i],emax=ebins[1][i],evclsmin=evmin,evclsmax=evmax,zmax=0,chatter=0) #zmax=0 makes no zenith selection
			else:
				gtselect.run(infile=ft1,outfile=out,ra=0,dec=0,rad=180,tmin=0,tmax=0,emin=ebins[0][i],emax=ebins[1][i],evclass=evclass,zmax=0,chatter=0)
		else:
			print '  -%s' %out,'already exists, skipping gtselect for this energy band.-'
	return

#this function will calcluate the exposure maps for each energy band, one exposure cube (same as that used for the full fit) is sufficient
def makeExpMaps(likeIn):
	gtexpmap=GtApp('gtexpmap') #note, checks for existence of band expMaps, won't overwrite
	Src=likeIn.source
	NBins=likeIn.NBins
	IRFs=likeIn.IRFs
	ft2=likeIn.ft2
	expCube=likeIn.expCube
	radius=likeIn.ft1ROI[2]+10.
	for i in range(0,likeIn.nbins):
		ev='%s_%ibins_band%i.fits' %(Src.replace(' ','_'),NBins,i)
		out='%s_%ibins_band%i_%s_em.fits' %(Src.replace(' ','_'),NBins,i,IRFs) #add IRFs name to file to be specific
		if(not os.access(out,os.F_OK)):
			gtexpmap.run(evfile=ev,scfile=ft2,expcube=expCube,outfile=out, irfs=IRFs,srcrad=radius,nlong=120,nlat=120,nenergies=20,chatter=0)
		else:
			print '  -%s' %out,'already exists, skipping gtexpmap for this energy band.-'
	return 

#this checks for the maximum energy event within the 95% PSF of the source
def maxEnergy(file,likeIn):
	data=file[1].data
	nrows=int(file[1].header['NAXIS2'])
	ENERGY=data.field('ENERGY')
	THETA=data.field('THETA')
	PHI=data.field('PHI')
	RA=data.field('RA')
	DEC=data.field('DEC')
	CNVTYPE=data.field('CONVERSION_TYPE')
	
	pyIrfLoader.Loader_go() #get available IRFs
	irfs={} #make front and back instances
	irfs[0]=pyIrfLoader.IrfsFactory.instance().create(likeIn.IRFs+'::FRONT')
	irfs[1]=pyIrfLoader.IrfsFactory.instance().create(likeIn.IRFs+'::BACK')
	 #make psf instances
	psf={}
	psf[0]=irfs[0].psf()
	psf[1]=irfs[1].psf()
	ra=likeIn.ubAn.model.srcs[likeIn.source].src.getSrcFuncs()['Position'].getParamValue('RA')
	dec=likeIn.ubAn.model.srcs[likeIn.source].src.getSrcFuncs()['Position'].getParamValue('DEC')
	EMax=[]
	Sep=angsep(ra,dec,RA,DEC)
	fangInt=psf[0].angularIntegral
	bangInt=psf[1].angularIntegral
	for s,c,r,d,e,t,p in zip(Sep,CNVTYPE,RA,DEC,ENERGY,THETA,PHI): #check containment of each event and find max energy
		if(c==0):
			cont=fangInt(float(e),float(t),float(p),s)
			if(cont<=0.95):
				EMax+=[e]
		if(c==1):
			cont=bangInt(float(e),float(t),float(p),s)
			if(cont<=0.95):
				EMax+=[e]
	return max(EMax)

#calculates the angular separation between two points on the sky
def angsep(ra1,dec1,ra2,dec2):
	ra1*=pi/180.
	dec1*=pi/180.
	ra2r=[r*pi/180. for r in ra2]
	dec2r=[d*pi/180. for d in dec2]
	cdec1=cos(dec1)
	sdec1=sin(dec1)
	cdec2=num.cos(dec2r)
	sdec2=num.sin(dec2r)
	diff=[ra1-r for r in ra2r]
	cdiff=num.cos(diff)
	SEP=num.arccos([cdec1*cd2*cdf+sdec1*sd2 for cd2,cdf,sd2 in zip(cdec2,cdiff,sdec2)])
	#SEP=num.arccos(num.cos(dec1)*num.cos(dec2)*num.cos(ra1-ra2)+num.sin(dec1)*num.sin(dec2)) #returns values between 0 and pi radians
	return [s*180./pi for s in SEP]

#creates lower and upper bounds for energy bins, takes into account maximum energy found wihtin 95% containment of the source
def plotEbounds(energies,maxE):
	i=len(energies)-1
	j=0
	while i>=0: #start at the top, check to see which energy bin contains maxE, corresponds to highest energy bin for plotting
		if(energies[i]>=maxE and energies[i-1]<=maxE):
			j=i
			i=-1
		i-=1
	#mins=energies[0:j]
	#maxs=energies[1:j+1]
	return energies[0:j],energies[1:j+1]

def scaleEst(likeIn,myEnergy):
	modEs=log_array(51,likeIn.FitBounds[0],likeIn.FitBounds[1])
	centEs=[0.5*(e1+e2) for e1,e2 in zip(modEs[0:-1],modEs[1:])]
	for i in range(0,len(centEs)-1):
		if centEs[i]<=myEnergy and centEs[i+1]>myEnergy:
			myindex=i
	return int(log10(likeIn.Fit[1][myindex]))


#runs likelihood for each band and records prefactors, values, and TS values for each, computes upperlimits if necessary
def runLike(likeIn,ecent,ftol,tslim,ttype,opt,rescaleAll,lastbinUL,wx):
	pts=[]
	errs=[]
	tsPts=[]
	gamma=[]
	Src=likeIn.source
	NBins=likeIn.NBins
	IRFs=likeIn.IRFs
	expCube=likeIn.expCube
	ft2=likeIn.ft2
	bandModel=likeIn.bandModel
	nbins=likeIn.nbins
	phCorr=likeIn.phCorr
	flux=likeIn.ubAn.flux
	for i in range(0,nbins): #first, set up for a likelihood run
		print '  -Runnng Likelihood for band%i-' %i
		ev='%s_%ibins_band%i.fits' %(Src.replace(' ','_'),NBins,i)
		em='%s_%ibins_band%i_%s_em.fits' %(Src.replace(' ','_'),NBins,i,IRFs)
		band_obs=UnbinnedObs((ev),ft2,irfs=IRFs,expMap=em,expCube=expCube)
		band_like=UnbinnedAnalysis(band_obs,bandModel,opt)
		band_like.setFitTolType(ttype)
		stype=band_like.model.srcs[Src].spectrum().genericName()
		emin,emax=band_like.observation.roiCuts().getEnergyCuts()
		if phCorr!=1:
			for src in band_like.sourceNames():
				par=band_like.normPar(src)
				par.setValue(par.getValue()*phCorr)
		#then set the scale factor to the center of the energy band, make sure it's frozen, and get index for prefactor while you're at it
		Ts=band_like.Ts
		freeze=band_like.freeze
		fit=band_like.fit
		DO=1
		if stype=='PowerLaw':
			scale=getParamIndx(band_like,Src,'Scale') #this is where you have to use PowerLaw, for PowerLaw2 these parameters don't exist and this will cause problems
			pref=getParamIndx(band_like,Src,'Prefactor')
			freeze(scale)
			band_like[scale].setBounds(20,5e5)
			band_like[scale].setScale(1)
			band_like[scale]=1000.*ecent[i] #put center energies in units of GeV but xml files use MeV
			#multiplier=band_like[pref].getScale() #need to get the scale of the prefactor so the values will not be too large
			try:
				logFlux=log10(flux(Src,emin=emin,emax=emax)/(emax-emin))
			except:
				logFlux=-14
			newScale=max(int(floor(logFlux)),-14)
			print newScale
			band_like[pref].setScale(10.**newScale)
			multiplier=10.**newScale
			#cycle through the other point sources and adjust parameters of free point sources with PowerLaw2 models
			for src in band_like.sourceNames():
				spec=band_like[src].funcs['Spectrum']
				par=spec.normPar()
				if par.isFree()==True and band_like.model.srcs[src].spectrum().genericName()=='PowerLaw2':
					HIGH=getParamIndx(band_like,src,'UpperLimit')
					LOW=getParamIndx(band_like,src,'LowerLimit')
					band_like[HIGH].setBounds(20,5e5) #just in case, make sure no out of range error gets thrown
					band_like[LOW].setBounds(20,5e5)
					band_like[HIGH].setScale(1.) #just in case, make sure scale is MeV
					band_like[LOW].setScale(1.)
					band_like[HIGH]=emax
					band_like[LOW]=emin
					freeze(HIGH) #just in case, make sure these aren't fit values
					freeze(LOW)
					if rescaleAll==True:
						try:
							logFlux=log10(flux(src,emin=emin,emax=emax))
						except:
							logFlux=-14
						newScale=max(int(floor(logFlux)),-14)
						par.setScale(10**newScale)
				if rescaleAll==True and src!=Src and par.isFree()==True and band_like.model.srcs[src].spectrum().genericName()=='PowerLaw':
					try:
						logFlux=log10(flux(src,emin=emin,emax=emax)/(emax-emin))
					except:
						logFlux=-14
					newScale=max(int(floor(logFlux)),-14)
					par.setScale(10**newScale)
			try:
				fit(tol=1,verbosity=0,optimizer=opt)
				fail=0			
			except:
				try:
					fit(tol=1*10,verbosity=0,optimizer=opt)
					fail=0
				except:
					try:
						fit(tol=1./10,verbosity=0,optimizer=opt)
						fail=0
					except:
						print "Fit with optimizer %s with tolerances ~1 to look for negative or zero TS sources failed, if error bars are unrealistically small you may need to redo the fit for energy band %i manually" %(opt,i)
						fail=1
			if fail==0:
				for src in band_like.sourceNames():
					if src!=Src:
						par=band_like.normPar(src)
						if band_like[src].type=='PointSource' and Ts(src)<=0 and par.isFree()==True:
							band_like.deleteSource(src)
							print "   -Removing %s from the model" %src
			scale=getParamIndx(band_like,Src,'Scale') 
			pref=getParamIndx(band_like,Src,'Prefactor')
			#do the actual fit
			try:
				fit(tol=ftol,verbosity=0)
			except:
				try:
					print 'Trying lower tolerance of %s for band%i.' %(ftol/10,i)
					fit(tol=ftol/10,verbosity=0)
				except:
					try:
						print 'Trying higher tolerance of %s for band%i.' %(ftol*10,i)
						fit(tol=ftol*10,verbosity=0)
					except:
						print 'No convergence for band%i, skipping.' %i
						pts+=[0]
						errs+=[0]
						tsPts+=[0]
						pass
			#get the prefactor, error, and ts values
			val=band_like[pref].value()
			err=band_like[pref].error()
			TS=Ts(Src)
			if(TS<tslim or (i==(nbins-1) and lastbinUL==True)): #calculate 95% upperlimit if source TS<tslim, 25 by default (corresponds to ~5sigma)
				try:
					freeze(pref)
					ul=UpperLimits(band_like)
					UL=ul[Src].compute(emin=emin,emax=emax)
					val=UL[1]
					err=0
					print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
				except:
					try:
						print '   Tyring higher tolerance of %s for band %i to get good starting point for upper limit calculations.' %(ftol*10,i)
						band_like[pref].setFree(1)
						fit(tol=ftol*10,verbosity=0)
						TS=Ts(Src)
						if(TS<tslim or (i==(nbins-1) and lastbinUL==True)):
							freeze(pref)
							ul=UpperLimits(band_like)
							UL=ul[Src].compute(emin=emin,emax=emax)
							val=UL[1]
							err=0
							print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
						else:
							val=band_like[pref].value()
							err=band_like[pref].error()
					except:
						try:
							print '   Tyring lower tolerance of %s for band %i to get good starting point for upper limit calculations.' %(ftol/10,i)
							band_like[pref].setFree(1)
							fit(tol=ftol/10,verbosity=0)
							TS=Ts(Src)
							if(TS<tslim or (i==(nbins-1) and lastbinUL==True)):
								freeze(pref)
								ul=UpperLimits(band_like)
								UL=ul[Src].compute(emin=emin,emax=emax)
								val=UL[1]
								err=0
								print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
							else:
								val=band_like[pref].value()
								err=band_like[pref].error()
						except:
							err=0
							print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, TS<%s' %tslim,'but UpperLimits computation failed.'
							print '          Quoting best fit value with zero error.'
			tsPts+=[TS]
			pts+=[val*multiplier/phCorr]
			errs+=[err*multiplier/phCorr]
		elif stype=='PowerLaw2':
			Upper=getParamIndx(band_like,Src,'UpperLimit')
			Lower=getParamIndx(band_like,Src,'LowerLimit')
			Integral=getParamIndx(band_like,Src,'Integral')
			Index=getParamIndx(band_like,Src,'Index')
			freeze(Upper)
			freeze(Lower)
			band_like[Upper].setBounds(20,5e5)
			band_like[Lower].setBounds(20,5e5)
			band_like[Lower].setScale(1)
			band_like[Lower].setScale(1)
			band_like[Upper]=emax
			band_like[Lower]=emin
			#multiplier=band_like[Integral].getScale()
			try:
				logFlux=log10(flux(Src,emin=emin,emax=emax))
			except:
				logFlux=-14
			newScale=max(int(floor(logFlux)),-14)
			band_like[Integral].setScale(10**newScale)
			multiplier=10**newScale
			indxMult=band_like[Index].getScale()
			for src in band_like.sourceNames():
				if src!=Src:
					spec=band_like[src].funcs['Spectrum']
					par=spec.normPar()
					if par.isFree()==True and band_like.model.srcs[src].spectrum().genericName()=='PowerLaw2':
						HIGH=getParamIndx(band_like,src,'UpperLimit')
						LOW=getParamIndx(band_like,src,'LowerLimit')
						band_like[HIGH].setBounds(20,5e5) #just in case, make sure no out of range error gets thrown
						band_like[LOW].setBounds(20,5e5)
						band_like[HIGH].setScale(1.) #just in case, make sure scale is MeV
						band_like[LOW].setScale(1.)
						band_like[HIGH]=emax
						band_like[LOW]=emin
						freeze(HIGH) #just in case, make sure these aren't fit values
						freeze(LOW)
						if rescaleAll==True:
							try:
								logFlux=log10(flux(src,emin=emin,emax=emax))
							except:
								logFlux=-14
							newScale=max(int(floor(logFlux)),-14)
							par.setScale(10**newScale)
				if rescaleAll==True and src!=Src and par.isFree()==True and band_like.model.srcs[src].spectrum().genericName()=='PowerLaw':
					try:
						logFlux=log10(flux(src,emin=emin,emax=emax)/(emax-emin))
					except:
						logFlux=-14
					newScale=max(int(floor(logFlux)),-14)
					par.setScale(10**newScale)
			try:
				fit(tol=1,verbosity=0,optimizer=opt)
				fail=0			
			except:
				try:
					fit(tol=1*10,verbosity=0,optimizer=opt)
					fail=0
				except:
					try:
						fit(tol=1./10,verbosity=0,optimizer=opt)
						fail=0
					except:
						print "Fit with optimizer %s with tolerances ~1 to look for negative or zero TS sources failed, if error bars are unrealistically small you may need to redo the fit for energy band %i manually" %(opt,i)
						fail=1
			if not fail:
				for src in band_like.sourceNames():
					if src!=Src:
						par=band_like.normPar(src)
						if band_like[src].type=='PointSource' and Ts(src)<=0 and par.isFree()==True:
							band_like.deleteSource(src)
							print "   -Removing %s from the model" %src
			Upper=getParamIndx(band_like,Src,'UpperLimit')
			Lower=getParamIndx(band_like,Src,'LowerLimit')
			Integral=getParamIndx(band_like,Src,'Integral')
			Index=getParamIndx(band_like,Src,'Index')
			try:
				fit(tol=ftol,verbosity=0)
			except:
				try:
					print 'Trying lower tolerance of %s for band%i.' %(ftol/10,i)
					fit(tol=ftol/10,verbosity=0)
				except:
					try:
						print 'Trying higher tolerance of %s for band%i.' %(ftol*10,i)
						fit(tol=ftol*10,verbosity=0)
					except:
						print 'No convergence for band%i, skipping.' %i
						pts+=[0]
						errs+=[0]
						tsPts+=[0]
						gamma+=[0]
						DO=0
			if DO:
				val=band_like[Integral].value()
				err=band_like[Integral].error()
				TS=Ts(Src)
				gam=band_like[Index].value()*indxMult*-1.
				if(TS<tslim  or (i==(nbins-1) and lastbinUL==True)): #calculate 95% upperlimit if source TS<9, i.e. less than 3 sigma detection in each energy band
					try:
						freeze(Integral)
						ul=UpperLimits(band_like)
						UL=ul[Src].compute(emin=emin,emax=emax)
						val=UL[1]
						err=0
						#need to redo the fit with band_like set to Upper Limit value to get correct spectral index for that value
						band_like[Integral]=val
						freeze(Integral)
						fit(tol=ftol,verbosity=0)
						gam=band_like[Index].value()*indxMult*-1.
						print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
					except:
						try:
							print '   Tyring higher tolerance of %s for band %i to get good starting point for upper limit calculations.' %(ftol*10,i)
							band_like[Integral].setFree(1)
							fit(tol=ftol*10,verbosity=0)
							TS=Ts(Src)
							if(TS<tslim  or (i==(nbins-1) and lastbinUL==True)):
								freeze(Integral)
								ul=UpperLimits(band_like)
								UL=ul[Src].compute(emin=emin,emax=emax)
								val=UL[1]
								err=0
								band_like[Integral]=val
								freeze(Integral)
								fit(tol=ftol*10,verbosity=0)
								gam=band_like[Index].value()*indxMult*-1.
								print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
							else:
								val=band_like[Integral].value()
								err=band_like[Integral].error()
								gam=band_like[Index].value()*indxMult*-1.					
						except:
							try:
								print '   Tyring lower tolerance of %s for band %i to get good starting point for upper limit calculations.' %(ftol/10,i)
								band_like[Integral].setFree(1)
								fit(tol=ftol/10,verbosity=0)
								TS=Ts(Src)
								if(TS<tslim  or (i==(nbins-1) and lastbinUL==True)):
									freeze(Integral)
									ul=UpperLimits(band_like)
									UL=ul[Src].compute(emin=emin,emax=emax)
									val=UL[1]
									err=0
									band_like[Integral]=val
									freeze(Integral)
									fit(tol=ftol/10,verbosity=0)
									gam=band_like[Index].value()*indxMult*-1.
									print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
								else:
									val=band_like[Integral].value()
									err=band_like[Integral].error()
									gam=band_like[Index].value()*indxMult*-1.
							except:
								err=0
								print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, TS<%s' %tslim,'but UpperLimits computation failed.'
								print '          Quoting best fit value with zero error.'
				tsPts+=[TS]
				pts+=[val*multiplier/phCorr]
				errs+=[err*multiplier/phCorr]
				gamma+=[gam]
		else:
			print '%s needs to have PowerLaw or PowerLaw2 spectral model, not %s' %(Src,stype)
			print 'exiting without running likelihood in the energy bands'
			return None,None,None,None
		if wx:
			band_like.writeXml('%s_%ibins_band%i_fitmodel.xml'%(Src.replace(' ','_'),NBins,i))
		del band_like
		del band_obs
	return pts,errs,tsPts,gamma

#function gets index for a specific parameter for a specific source from model in UnbinnedAnalysis object fit
def getParamIndx(fit,name,NAME):
	ID=-1
	spec=fit[name].funcs['Spectrum']
	for indx, parName in zip(spec._parIds, spec.paramNames):
		if(parName==NAME):
			ID = indx
	if(ID==-1):
		print 'Parameter %s not found for source %s in file %s.' %(NAME,name,fit.srcModel)
	return ID

#will make and plot the actual plots without bowtie plot
def plotSED(sed,plot):
	gStyle.SetTextFont(132)
	gStyle.SetTitleFont(132,'xyz')
	gStyle.SetLabelFont(132,'xyz')
	minEs=sed.likeIn.bins[0]
	maxEs=sed.likeIn.bins[1]
	ecent=sed.centers
	fluxPts=sed.data[0]
	fluxErrs=sed.data[1]
	tsbands=sed.data[2]
	gammas=sed.data[3]
	gStyle.SetOptStat(0)
	gStyle.SetOptTitle(0)
	gStyle.SetPaperSize(14,10)
	XTitle='Energy (GeV)'
	YTitle=['dN/dE (cm^{-2} s^{-1} GeV^{-1} )','E^{2}dN/dE (erg cm^{-2} s^{-1} )','#nuF_{#nu} (erg cm^{-2} s^{-1} )']	
	#let's make the graphs for the model, both the counts spectrum and the E^2dN/dE spectrum
	modBins=log_array(51,sed.FitBounds[0]/1000.,sed.FitBounds[1]/1000.)
	modelEcent=[0.5*(e1+e2) for e1,e2 in zip(modBins[:-1],modBins[1:])]
	#for i in range(0,len(modBins)-1):
	#	modelEcent+=[0.5*(modBins[i]+modBins[i+1])]
	modEs=num.array(modelEcent)
	modFs=num.array(sed.Fit)
	modesq=[y*x**2*0.001602 for x,y in zip(modelEcent,sed.Fit)]
#	for x,y in zip(modelEcent,sed.Fit):
#		modesq+=[y*x**2*0.001602]
	modEsq=num.array(modesq)
	modcnt=TGraph(len(modEs),modEs,modFs)
	modflx=TGraph(len(modEs),modEs,modEsq)
	modcnt.GetXaxis().SetTitle(XTitle)
	modflx.GetXaxis().SetTitle(XTitle)
	modcnt.GetYaxis().SetTitle(YTitle[0])
	modflx.GetYaxis().SetTitle(YTitle[1])
	modcnt.GetXaxis().CenterTitle()
	modflx.GetXaxis().CenterTitle()
	modcnt.GetYaxis().CenterTitle()
	modflx.GetYaxis().CenterTitle()
	if(len(sed.BT)>1):
		up=[]
		down=[]
		for x,y,z in zip(modesq,sed.BT,modelEcent):
			up+=[x+y*z**2*0.001602]
			down+=[x-y*z**2*0.001602]
		Up=num.array(up)
		Down=num.array(down)
		top=TGraph(len(modEs),modEs,Up)
		bot=TGraph(len(modEs),modEs,Down)
		top.SetLineColor(2)
		bot.SetLineColor(2)
		front=TLine(modEs[0],down[0],modEs[0],up[0])
		back=TLine(modEs[-1],down[-1],modEs[-1],up[-1])
		front.SetLineColor(2)
		back.SetLineColor(2)
	plotMax=max(maxEs[-1],sed.FitBounds[1])
	plotMin=min(minEs[0],sed.FitBounds[0])
	cntdummy=TH1F('cntdummy',"",1000,plotMin*0.80/1000.,plotMax*1.2/1000.)
	cntdummy.SetXTitle(XTitle)
	cntdummy.SetYTitle(YTitle[0])
	cntdummy.GetXaxis().CenterTitle()
	cntdummy.GetYaxis().CenterTitle()
	modLow=min(sed.Fit)
	modHigh=max(sed.Fit)
	#do some prep for the data plots
	#if using PowerLaw 2 model, need to make a few adjustments to the data points and errors
	if len(sed.data[3])==0:
		GeVfluxPts=[float(x*1000.) for x in fluxPts if x!=0]#need to check for points with zero flux
		GeVfluxErrs=[float(x*1000.) for x,X in zip(fluxErrs,fluxPts) if X!=0]#means energy band fit failed
	else:
		newfluxPts=[]
		newfluxErrs=[]
		for m,M,f,e,g,c in zip(minEs,maxEs,fluxPts,fluxErrs,sed.data[3],ecent):
			if g!=1 and f!=0:
				newfluxPts+=[f*(-g+1)*(c*1000.)**(-g)/(M**(-g+1)-m**(-g+1))] #works out to be prefactor assuming Scale parameter (E0) of c (center energy of bin)
				newfluxErrs+=[e*(-g+1)*(c*1000.)**(-g)/(M**(-g+1)-m**(-g+1))]#similar for the errors
			else:
				if f!=0:
					newfluxPts+=[f/(c*log(M/m))]
					newfluxErrs+=[e/(c*log(M/m))]
		GeVfluxPts=[float(x*1000.) for x in newfluxPts]
		GeVfluxErrs=[float(x*1000.) for x in newfluxErrs]
	Fluxes=num.array(GeVfluxPts)
	FErrors=num.array(GeVfluxErrs)
	dataLow=min(GeVfluxPts)
	dataHigh=max(GeVfluxPts)
	cntdummy.GetYaxis().SetRangeUser(min(dataLow,modLow)*0.25,max(dataHigh,modHigh)*5)
	#now, lets make the data graphs
	cntArrows=[] #create lists to be filled with TArrows for upper limits points for dN/dE spectrum, TArrows not saved, only used if plot==True 
	flxArrows=[] #same as above but for E^2dN/dE spectrum
	newecent=[]
	for e,f in zip(ecent,fluxPts):
		if f!=0:
			newecent+=[e]
	energies=num.array(newecent)
	esqf=[]
	esqfErr=[]
	GeVBands=[float(x/1000.) for x in minEs]
	GeVBands+=[float(maxEs[-1]/1000.)]
	Bands=num.array(GeVBands) #make an array that TH1F will recognize for binning purposes
	if(len(sed.data[3])==0):
		for x,y,z in zip(energies,GeVfluxPts,GeVfluxErrs):
			if y!=0:
				esqf+=[y*x**2*0.001602]
				esqfErr+=[z*x**2*0.001602]
	else:
		for x,y,z,g,e,E in zip(ecent,fluxPts,fluxErrs,gammas,minEs,maxEs):
			if g!=1 and y!=0:
				esqf+=[y*x**(2-g)*(g-1)/((e/1000.)**(1-g)-(E/1000.)**(1-g))*0.001602]
				esqfErr+=[z*x**(2-g)*(g-1)/((e/1000.)**(1-g)-(E/1000.)**(1-g))*0.001602]
			else:
				if y!=0:
					esqf+=[y/(x*log(E/e))*(E*log(E)-e*log(e)-E+e)]
					esqfErr+=[z/(x*log(E/e))*(E*log(E)-e*log(e)-E+e)]
	EsqFluxes=num.array(esqf)
	EsqFErrors=num.array(esqfErr)
	eminus=[]
	eplus=[]
	newminEs=[]
	newmaxEs=[]
	for e,E,f in zip(minEs,maxEs,fluxPts):
		if f!=0:
			newminEs+=[e]
			newmaxEs+=[E]
	for x,y,z in zip(energies,newminEs,newmaxEs):
		eminus+=[x-(y/1000.)]
		eplus+=[(z/1000.)-x]
	Eminus=num.array(eminus)
	Eplus=num.array(eplus)
	cnt=TGraphAsymmErrors(len(energies),energies,Fluxes,Eminus,Eplus,FErrors,FErrors)
	flx=TGraphAsymmErrors(len(energies),energies,EsqFluxes,Eminus,Eplus,EsqFErrors,EsqFErrors)
	for x,y,z in zip(GeVfluxPts,GeVfluxErrs,ecent):
		if x>0 and y==0:
			cntArrows+=[TArrow(z,x,z,x*0.50,0.025,'|->')]
	for x,y,z in zip(esqf,esqfErr,ecent):
		if x>0 and y==0:
			flxArrows+=[TArrow(z,x,z,x*0.40,0.025,'|->')]
	cnt.GetXaxis().SetTitle(XTitle)
	flx.GetXaxis().SetTitle(XTitle)
	cnt.GetYaxis().SetTitle(YTitle[0])
	if len(sed.data[3])==0:
		flx.SetTitle(YTitle[1])
	else:
		flx.SetTitle(YTitle[2])
	cnt.GetXaxis().CenterTitle()
	flx.GetXaxis().CenterTitle()
	cnt.GetYaxis().CenterTitle()
	flx.GetYaxis().CenterTitle()
	cnt.SetMarkerStyle(21)
	flx.SetMarkerStyle(21)
	cnt.SetMarkerSize(0.7)
	flx.SetMarkerSize(0.7)
	flx.GetYaxis().SetTitleOffset(1.25)
	cnt.GetYaxis().SetTitleOffset(1.25)
	#also, make TS histogram, can be interesting
	tshist=TH1F("tshist","TS values",len(Bands)-1,Bands)
	for i in range(1,len(Bands)):
		if fluxPts[i-1]!=0:
			tshist.SetBinContent(i,float(tsbands[i-1]))
		else:
			tshist.SetBinContent(i,-1)
	tshist.SetXTitle(XTitle)
	tshist.SetYTitle(sed.likeIn.source+" TS")
	tshist.GetXaxis().CenterTitle()
	tshist.GetYaxis().CenterTitle()
	#make the canvases for these guys
	if(len(sed.BT)>1):
		high_mod=max(up)
		low_mod=min(down)
	else:
		hold=[x*z**2*0.001602 for x,z in zip(sed.Fit,modelEcent)]
		#for x,z in zip(sed.Fit,modelEcent):
		#	hold+=[x*z**2*0.001602]
		high_mod=max(hold)
		low_mod=min(hold)
	high_data=max(esqf)
	low_data=min(esqf)
	high=max(high_data,high_mod)
	low=min(low_data,low_mod)
	if low<=0:
		low=high*1e-3
	flxdummy=TH1F("flxdummy","",100,plotMin*0.80/1000.,plotMax*1.2/1000.) #make a dummy histogram with slightly larger to see more of bowtie shape
	flxdummy.GetYaxis().SetRangeUser(low*0.25,high*5) #set the  y range to avoid having zero as min, bad for log axis
	flxdummy.SetXTitle(XTitle)
	if len(sed.data[3])==0:
		flxdummy.SetYTitle(YTitle[1])
	else:
		flxdummy.SetYTitle(YTitle[2])
	flxdummy.GetXaxis().CenterTitle()
	flxdummy.GetYaxis().CenterTitle()
	flxdummy.GetXaxis().SetTitleOffset(1.25)
	flxdummy.GetYaxis().SetTitleOffset(1.25)
	cntCanvas=TCanvas("cntCanvas","Count Spectrum",700,500)
	flxCanvas=TCanvas("flxCanvas","E^{2}dN/dE Spectrum",700,500)
	tsCanvas=TCanvas("tsCanvas","Ts vs Energy",700,500)
	cntCanvas.SetFillColor(0)
	cntCanvas.SetFrameBorderMode(0)
	cntCanvas.SetBorderMode(0)
	cntCanvas.SetTicks(1,1)
	flxCanvas.SetFillColor(0)
	flxCanvas.SetFrameBorderMode(0)
	flxCanvas.SetBorderMode(0)
	flxCanvas.SetTicks(1,1)
	tsCanvas.SetFillColor(0)
	tsCanvas.SetFrameBorderMode(0)
	tsCanvas.SetBorderMode(0)
	tsCanvas.SetTicks(1,1)
	cntCanvas.SetLogx()
	cntCanvas.SetLogy()
	flxCanvas.SetLogx()
	flxCanvas.SetLogy()
	tsCanvas.SetLogx()
		#make TLegends for the two spectral plots, don't really need two but doesn't hurt to make one for each in the event that one wants to change one or the other
	cntLegend=TLegend(0.65,0.7,0.85,0.85) #the good thing about these is that they are movable on the graph in the even that the placement covers anything up
	flxLegend=TLegend(0.65,0.7,0.85,0.85)
	cntLegend.SetBorderSize(1)
	flxLegend.SetBorderSize(1)
	cntLegend.SetFillColor(0)
	flxLegend.SetFillColor(0)
	cntLegend.AddEntry(cnt,"Energy Band Fits",'p')
	cntLegend.AddEntry(modcnt,"Maximum Likelihood Model",'l')
	flxLegend.AddEntry(flx,"Energy Band Fits",'p')
	flxLegend.AddEntry(modflx,"Maximum Likelihood Model",'l')
	#now draw the histograms on the canvases, do data histograms first so in the likely event that model extends beyond data, won't be shown
	cntCanvas.cd()
	ctps=TPostScript('%s_%ibins_cntSpec.eps' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),113)
	cntdummy.Draw()
	cnt.Draw("psame") #draw as points with capped error bars
	modcnt.Draw("csame") #draw on same canvas as a curve
	cntLegend.Draw()
	for x in cntArrows:
		x.Draw()
	ctps.Close()
	flxCanvas.cd()
	flxps=TPostScript('%s_%ibins_EsqdNdE.eps' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),113)
	flxdummy.Draw()
	flx.Draw("psame")
	modflx.Draw("csame")
	if(len(sed.BT)>1):
		top.Draw("csame")
		bot.Draw("csame")
		front.Draw()
		back.Draw()
	flxLegend.Draw()
	for x in flxArrows:
		x.Draw()
	flxps.Close()
	tsCanvas.cd()
	tsps=TPostScript('%s_%ibins_TSvEnergy.eps' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),113)
	tshist.Draw()
	tsps.Close()
	if(plot==True):
		SetOwnership(cntCanvas,False)
		SetOwnership(modcnt,False)
		SetOwnership(cnt,False)
		for x in cntArrows:
			SetOwnership(x,False)
		SetOwnership(cntLegend,False)
		SetOwnership(flxCanvas,False)
		SetOwnership(flxdummy,False)
		SetOwnership(cntdummy,False)
		SetOwnership(modflx,False)
		SetOwnership(flx,False)
		for x in flxArrows:
			SetOwnership(x,False)
		SetOwnership(flxLegend,False)
		if(len(sed.BT)>1):
			SetOwnership(top,False)
			SetOwnership(bot,False)
			SetOwnership(front,False)
			SetOwnership(back,False)
		SetOwnership(tsCanvas,False)
		SetOwnership(tshist,False)
	histFile=TFile(sed.likeIn.source.replace(' ','_')+"_%ibins_histos.root" %sed.likeIn.NBins,"RECREATE") #This writes the histograms to a file for use later, note that it is set to overwrite existing files
	tshist.Write()
	cntCanvas.Write()
	flxCanvas.Write()
	tsCanvas.Write()
	writeSpecFits(sed,modBins,modelEcent,[x/1000. for x in newminEs],[x/1000. for x in newmaxEs],GeVfluxPts,GeVfluxErrs,newecent,EsqFluxes,EsqFErrors) #write the output to a sed.Fits file
	return 

#create fits file with spectral output
def writeSpecFits(sed,modelEbins,modelEcent,minEs,maxEs,GeVfluxPts,GeVfluxErrs,ecent,esqf,esqferr):
	#make arrays, columns, and table for the fit
	fitArray=num.array(sed.Fit)
	modEcentArray=num.array(modelEcent) #this is necessary to make E^2dN/dE spectrum later
	if(len(sed.BT)>1):
		btArray=num.array(sed.BT)
		btCol=pyfits.Column(name='BowTie',format='E',unit='ph/cm^2/s/GeV',array=btArray)
	fitCol=pyfits.Column(name='Model',format='E',unit='ph/cm^2/s/GeV',array=fitArray)
	modEcentCol=pyfits.Column(name='Center Energy',format='E',unit='GeV',array=modEcentArray)
	if(len(sed.BT)>1):
		fitTable=pyfits.new_table([modEcentCol,fitCol,btCol])
	else:
		fitTable=pyfits.new_table([modEcentCol,fitCol])
	fitTable.name='MODEL FLUX'
	#make arrays, colums, and table for the data
	EcentArray=num.array(ecent)
	fluxArray=num.array(GeVfluxPts)
	errArray=num.array(GeVfluxErrs)
	esqArray=num.array(esqf)
	esqerrArray=num.array(esqferr)
	EcentCol=pyfits.Column(name='Center Energy',format='E',unit='GeV',array=EcentArray)
	fluxCol=pyfits.Column(name='dN/dE',format='E',unit='ph/cm^2/s/GeV',array=fluxArray)
	errCol=pyfits.Column(name='dN/dE_Error',format='E',unit='ph/cm^2/s/GeV',array=errArray)
	esqCol=pyfits.Column(name='E^2dN/dE',format='E',unit='erg/cm^2/s',array=esqArray)
	esqerrCol=pyfits.Column(name='E^2dN/dE_Error',format='E',unit='erg/cm^2/s',array=esqerrArray)
	dataTable=pyfits.new_table([EcentCol,fluxCol,errCol,esqCol,esqerrCol])
	dataTable.name='DATA FLUX'
	#make arrays, columns, and table for the model energy bins
	modLowArray=num.array(modelEbins[:-1])
	modHighArray=num.array(modelEbins[1:])
	modLowCol=pyfits.Column(name='Bin Low Edge',format='E',unit='GeV',array=modLowArray)
	modHighCol=pyfits.Column(name='Bin High Edge',format='E',unit='GeV',array=modHighArray)
	modBinsTable=pyfits.new_table([modLowCol,modHighCol])
	modBinsTable.name='Model Energy Bins'
	#make arrrays, columns, and table for the data energy bins
	dataLowArray=num.array(minEs)
	dataHighArray=num.array(maxEs)
	dataLowCol=pyfits.Column(name='Bin Low Edge',format='E',unit='GeV',array=dataLowArray)
	dataHighCol=pyfits.Column(name='Bin High Edge',format='E',unit='GeV',array=dataHighArray)
	dataBinsTable=pyfits.new_table([dataLowCol,dataHighCol])
	dataBinsTable.name='Data Energy Bins'
	primary=pyfits.PrimaryHDU() #make an empty primary header
	hdulist=pyfits.HDUList([primary,dataTable,dataBinsTable,fitTable,modBinsTable]) #make a list of tables for extensions
	hdulist.writeto('%s_%ibins_likeSEDout.fits' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),output_verify='ignore',clobber='yes') #write the file
	return

#check to make sure files exist, do first to avoid errors later
def fileCheck(file):
	if (not os.access(file,os.F_OK)):
		print "Error:  %s not found." %file
		return 0
	return 1

#make a residuals plot by accessing the model directly to compare point to point more exactly
def residualsPlot(sed,plot):
	ecent=sed.centers
	Emins=sed.likeIn.bins[0]
	Emaxs=sed.likeIn.bins[1]
	fluxPts=sed.data[0]
	fluxErrs=sed.data[1]
	eminus=[]
	eplus=[]
	for x,y,z in zip(ecent,Emins,Emaxs):
		eminus+=[x-(y/1000)]
		eplus+=[(z/1000)-x]
	Eminus=num.array(eminus)
	Eplus=num.array(eplus)
	modelPts=[]
	modelErrs=[]
	#get the model
	mysrc=pyLike.PointSource_cast(sed.likeIn.ubAn[sed.likeIn.source].src)
	for x in ecent: #using the same center energies as used in the data
		MeV=x*1000. #ecent is in GeV, but mysrc.spectrum assumes MeV
		arg=pyLike.dArg(MeV)
		val=mysrc.spectrum()(arg)
		modelPts+=[float(val)] 
	#allow for the possibility of the model points have errors too
	if(sed.likeIn.ubAn.covariance is None):
		modelErrs=[0]*len(modelPts)
	else:
		covArray=num.array(sed.likeIn.ubAn.covariance)
		srcCovArray=[]
		par_index_map={}
		indx=0
		for src in sed.likeIn.ubAn.sourceNames():
			parNames=pyLike.StringVector()
			sed.likeIn.ubAn[src].src.spectrum().getFreeParamNames(parNames)
			for par in parNames:
				par_index_map['::'.join((src,par))]=indx
				indx +=1
		srcPars=pyLike.StringVector()
		sed.likeIn.ubAn[sed.likeIn.source].src.spectrum().getFreeParamNames(srcPars)
		pars=['::'.join((sed.likeIn.source,x)) for x in srcPars]
		for xpar in pars:
			ix=par_index_map[xpar]
			srcCovArray.append([covArray[ix][par_index_map[ypar]] for ypar in pars])
		cov=num.array(srcCovArray)
		#the whole point here is to get the srcCovArray
		for x in ecent:
			MeV=x*1000.
			arg=pyLike.dArg(MeV)
			partials=num.array([mysrc.spectrum().derivByParam(arg,y) for y in srcPars])
			val=num.sqrt(num.dot(partials,num.dot(cov,partials))) #these should come out same as the model so convert to ph/cm^2/s/GeV as well
			modelErrs+=[float(val)]
	#calculate the residuals
	resids=[]
	residErrs=[]
	if len(sed.data[3])==0:
		for a,b,c,d in zip(fluxPts,fluxErrs,modelPts,modelErrs):
			resids+=[(a/c)-1.] #do (data-model)/model which works out to (data/model)-1
			residErrs+=[(a/c)*num.sqrt((d/c)**2+(b/a)**2)]
	else:
		newfluxPts=[]
		newfluxErrs=[]
		for m,M,f,e,g,c in zip(Emins,Emaxs,fluxPts,fluxErrs,sed.data[3],ecent):
			newfluxPts+=[f*(-g+1)*(c*1000.)**(-g)/(M**(-g+1)-m**(-g+1))] #works out to be prefactor assuming Scale parameter (E0) of c (center energy of bin)
			newfluxErrs+=[e*(-g+1)*(c*1000.)**(-g)/(M**(-g+1)-m**(-g+1))]#similar for the errors
		for a,b,c,d in zip(newfluxPts,newfluxErrs,modelPts,modelErrs):
			resids+=[(a/c)-1.] #do (data-model)/model which works out to (data/model)-1
			residErrs+=[(a/c)*num.sqrt((d/c)**2+(b/a)**2)]
	#now lets make the plots
	Resids=num.array(resids)
	RErrors=num.array(residErrs)
	energies=num.array(ecent)
	rgraph=TGraphAsymmErrors(len(energies),energies,Resids,Eminus,Eplus,RErrors,RErrors)
	gStyle.SetOptStat(0)
	gStyle.SetOptTitle(0)
	rhist=TH1F("rhist","",100,Emins[0]/1000.,Emaxs[-1]/1000.)
	rhist.GetYaxis().SetRangeUser(-1.,1.)
	rhist.SetLineStyle(2)
	rhist.SetFillColor(0)
	rhist.SetXTitle('Energy (GeV)')
	rhist.SetYTitle('(data-model)/model')
	rhist.GetXaxis().CenterTitle()
	rhist.GetYaxis().CenterTitle()
	rhist.SetMarkerStyle(20)
	rhist.SetMarkerSize(0.6)
	rcan=TCanvas('rcan','Residuals',700,350)
	rcan.SetFillColor(0)
	rcan.SetFrameBorderMode(0)
	rcan.SetBorderMode(0)
	rcan.SetTicks(1,1)
	rcan.SetLogx()
	#save as a .eps file
	rps=TPostScript('%s_%ibins_residuals.eps' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),113)
	rcan.cd()
	rhist.Draw('')
	rgraph.Draw('psame')
	rps.Close()
	if(plot==True):
		SetOwnership(rcan,False)
		SetOwnership(rhist,False)
		SetOwnership(rgraph,False)
	#save as a root file
	rfile=TFile('%s_%ibins_residuals.root' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),'RECREATE')
	rcan.Write()
	return

#run likelihood on the full ft1 to get the model and, if CoVar==True, the bowtie plot.
def fullLike(likeIn,CoVar,FTol,minE,maxE,wx=False):
	#make some energy bounds for the fit, same max and min as for the bands before but with more bins
	modEs=log_array(51,minE,maxE)
	centEs=[0.5*(e1+e2) for e1,e2 in zip(modEs[0:-1],modEs[1:])]
	#for i in range(0,len(modEs)-1):
	#	centEs+=[0.5*(modEs[i]+modEs[i+1])]
	#do the full fit
	phCorr=likeIn.phCorr
	ubAn=likeIn.ubAn
	for src in ubAn.sourceNames():
		par=ubAn.normPar(src)
		par.setValue(par.getValue()*phCorr)
		#adjusts values to account for phase cut
	try:
		ubAn.fit(tol=FTol,covar=CoVar,verbosity=0)
	except:
		try:
			print 'Trying lower tolerance of %s for full fit.' %FTol/10.
			ubAn.fit(tol=FTol/10.,covar=CoVar,verbosity=0)
		except:
			try:
				print 'Trying higher tolerance of %s for full fit.' %FTol*10.
				ubAn.fit(tol=FTol*10.,covar=CoVar,verbosity=0)
			except:
				print 'Error, no convergence.'
				return
	
	for src in ubAn.sourceNames():
		par=ubAn.normPar(src)
		par.setValue(par.getValue()/phCorr)
		par.setError(par.error()/phCorr)
		#updates the values and errors of the normalization parameters...will this adequately account for bowtie info?		
	
	#most of the following (getting model and bowtie) is taken directly from David Sanchez's pyUnfoldPlot with some minor stylistic changes
	
	#get the model
	mysrc=pyLike.PointSource_cast(ubAn[likeIn.source].src)
	#spec=[]
	spec=[float(1000.*mysrc.spectrum()(pyLike.dArg(x))) for x in centEs]
	#for x in centEs:
	#	arg=pyLike.dArg(x)
	#	val=mysrc.spectrum()(arg) #gives (I believe) dN/dE spectrum in ph/cm^2/s/MeV so need to convert to ph/cm^2/s/GeV
	#	spec+=[float(1000.*val)]
	
	if(ubAn.covariance is None):
		bt=[0]
	
	else:
		bt=[]
		covArray=num.array(ubAn.covariance)
		srcCovArray=[]
		par_index_map={}
		indx=0
		for src in ubAn.sourceNames():
			parNames=pyLike.StringVector()
			likeIn.ubAn[src].src.spectrum().getFreeParamNames(parNames)
			for par in parNames:
				par_index_map['::'.join((src,par))]=indx
				indx +=1
		srcPars=pyLike.StringVector()
		ubAn[likeIn.source].src.spectrum().getFreeParamNames(srcPars)
		pars=['::'.join((likeIn.source,x)) for x in srcPars]
		for xpar in pars:
			ix=par_index_map[xpar]
			srcCovArray.append([covArray[ix][par_index_map[ypar]] for ypar in pars])
		cov=num.array(srcCovArray)
		#the whole point here is to get the srcCovArray
		for x in centEs:
			arg=pyLike.dArg(x)
			partials=num.array([mysrc.spectrum().derivByParam(arg,y) for y in srcPars])
			val=num.sqrt(num.dot(partials,num.dot(cov,partials))) #these should come out same as the model so convert to ph/cm^2/s/GeV as well
			bt+=[float(1000.*val)]
	myfile=open('likeSED_%s_fullFitout.txt'%likeIn.source.replace(' ','_'),'w')
	print 'Full energy range model for %s:' %likeIn.source
	myfile.write('Full energy range model for %s:\n' %likeIn.source)
	print ubAn[likeIn.source]
	myfile.write('%s\n'%ubAn[likeIn.source])
	if ubAn.covariance is None:
		myfile.write('(Covariance Matrix not calculated)\n')
		print 'Flux %.1f-%.1f GeV %.1e cm^-2 s^-1' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]))
		myfile.write('Flux %.1f-%.1f GeV %.1e cm^-2 s^-1\n' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1])))
	else:
		print 'Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]),ubAn.fluxError(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]))
		myfile.write('Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1\n' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]),ubAn.fluxError(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1])))
	print "Test Statistic",ubAn.Ts(likeIn.source)
	myfile.write('Test Statistic %.2f'%ubAn.Ts(likeIn.source))
	myfile.close()
	if wx:
		ubAn.writeXml('%s_fullErange_fitmodel.xml'%likeIn.source.replace(' ','_'))
	return bt,spec

#get the model and bowtie (if covariance matrix has been calculated) from the full energy range fit if likelihood has already been run.
#this is pretty much the same as the fullLike function, just without the LIKE.fit() bit
def getFit(likeIn,minE,maxE,numBins,prtMod,ENERGIES=[],expCorrect=False,wx=False):
	#make some energy bounds for the fit, same max and min as for the bands before but with more bins
	if ENERGIES==[]:
		modEs=log_array(numBins,minE,maxE)
	else:
		modEs=ENERGIES
	centEs=[0.5*(e1+e2) for e1,e2 in zip(modEs[0:-1],modEs[1:])]
	#for i in range(0,len(modEs)-1):
	#	centEs+=[0.5*(modEs[i]+modEs[i+1])]
	#most of the following (getting model and bowtie) is taken directly from David Sanchez's pyUnfoldPlot with some minor stylistic changes
	
	#check if one needs to do exposure correction to account for not using the full phase
	phCorr=likeIn.phCorr
	ubAn=likeIn.ubAn
	if expCorrect:
		for src in ubAn.sourceNames():
			par=ubAn.normPar(src)
			par.setValue(par.getValue()/phCorr)
			par.setError(par.error()/phCorr)
		#updates the values and errors of the normalization parameters...will this adequately account for bowtie info?
	
	#get the model
	mysrc=pyLike.PointSource_cast(likeIn.ubAn[likeIn.source].src)
	spec=[float(1000.*mysrc.spectrum()(pyLike.dArg(x))) for x in centEs]
	#for x in centEs:
	#	arg=pyLike.dArg(x)
	#	val=mysrc.spectrum()(arg) #gives (I believe) dN/dE spectrum in ph/cm^2/s/MeV so need to convert to ph/cm^2/s/GeV
	#	spec+=[float(1000.*val)]
	
	if(likeIn.ubAn.covariance is None):
		bt=[0]
	
	else:
		bt=[]
		covArray=num.array(likeIn.ubAn.covariance)
		srcCovArray=[]
		par_index_map={}
		indx=0
		for src in ubAn.sourceNames():
			parNames=pyLike.StringVector()
			ubAn[src].src.spectrum().getFreeParamNames(parNames)
			for par in parNames:
				par_index_map['::'.join((src,par))]=indx
				indx +=1
		srcPars=pyLike.StringVector()
		ubAn[likeIn.source].src.spectrum().getFreeParamNames(srcPars)
		pars=['::'.join((likeIn.source,x)) for x in srcPars]
		for xpar in pars:
			ix=par_index_map[xpar]
			srcCovArray.append([covArray[ix][par_index_map[ypar]] for ypar in pars])
		cov=num.array(srcCovArray)
		#the whole point here is to get the srcCovArray
		for x in centEs:
			arg=pyLike.dArg(x)
			partials=num.array([mysrc.spectrum().derivByParam(arg,y) for y in srcPars])
			val=num.sqrt(num.dot(partials,num.dot(cov,partials))) #these should come out same as the model so convert to ph/cm^2/s/GeV as well
			bt+=[float(1000.*val)]
	if prtMod==1:
		myfile=open('likeSED_%s_fullFitout.txt'%likeIn.source.replace(' ','_'),'w')
		print 'Full energy range model for %s:' %likeIn.source
		myfile.write('Full energy range model for %s:\n' %likeIn.source)
		print ubAn[likeIn.source]
		myfile.write('%s\n'%ubAn[likeIn.source])
		if ubAn.covariance is None:
			print 'Flux %.1f-%.1f GeV %.1e cm^-2 s^-1' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]))
			myfile.write('(Covariance Matrix not calculated)\n')
			myfile.write('Flux %.1f-%.1f GeV %.1e cm^-2 s^-1\n' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1])))
		else:
			print 'Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]),ubAn.fluxError(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]))
			myfile.write('Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1\n' %(likeIn.ft1EBounds[0]/1000.,likeIn.ft1EBounds[1]/1000.,ubAn.flux(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1]),ubAn.fluxError(likeIn.source,emin=likeIn.ft1EBounds[0],emax=likeIn.ft1EBounds[1])))
		print "Test Statistic",ubAn.Ts(likeIn.source)
		myfile.write('Test Statistic %.2f'%ubAn.Ts(likeIn.source))
		myfile.close()
		if wx:
			ubAn.writeXml('%s_fullErange_fitmodel.xml'%likeIn.source.replace(' ','_'))
	return bt,spec


#calculate weighted bin centers based on the best fit model
def specWeightedCenters(likeIn):
	centers=[]
	if likeIn.custBins:
		steps=[(M-m)/100. for m,M in zip(likeIn.bins[0],likeIn.bins[1])]
		modEs=[]
		for i in range(0,likeIn.nbins):
			for j in range(0,100):
				modEs+=[likeIn.bins[0][i]+steps[i]*j]
		modEs+=[likeIn.bins[1][-1]]
		spec=getFit(likeIn,likeIn.bins[0][0],likeIn.bins[1][-1],likeIn.nbins*100+1,0,modEs)[1]
		centEs=[0.5*(e1+e2)/1000. for e1,e2 in zip(modEs[:-1],modEs[1:])]
	else:
		modEs=log_array(likeIn.nbins*100+1,likeIn.bins[0][0]/1000.,likeIn.bins[1][-1]/1000.)
		spec=getFit(likeIn,likeIn.bins[0][0],likeIn.bins[1][-1],likeIn.nbins*100+1,0)[1]
		centEs=[0.5*(e1+e2) for e1,e2 in zip(modEs[:-1],modEs[1:])]
	for i in range(0,likeIn.nbins):
		j=i*100
		jmax=(i+1)*100
		centers+=[float(num.sum(num.multiply(spec[j:jmax],centEs[j:jmax]))/num.sum(spec[j:jmax]))]
	return centers