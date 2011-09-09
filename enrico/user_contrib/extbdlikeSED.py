#Define some classes

#this class takes the very basic inputs and gathers the necessary information that one isn't likley to change
#also define class functions to easily access the other info gathering functions
class extbdlikeInput:
	#the arguments are: 
	#LIKE (an BinnedAnalysis object),
	#model (xml model for fitting the energy bands), and
	#SrcName (name of your source of interest)
	#nbins (integer number of bins)
	#model and SrcNames are assumed to be strings
	
	def __init__(self,LIKE,ft1,cc,ft2,SrcName,model='',nbins=20,phCorr=1.0):
		if(not fileCheck(ft1)):
			return
		if(not fileCheck(cc)):
			return
		if(not fileCheck(ft2)):
			return
		if(model!=''): #if given a value for model, make sure file exists
			if(not fileCheck(model)):
				return
		if(model==''): #if no value given for model, set model to that in LIKE
			model=LIKE.srcModel
		self.bdAn=LIKE
		self.sm=LIKE.binnedData.srcMaps
		self.ft1=ft1
		self.ft2=ft2
		self.cc=cc
		self.expMap=LIKE.binnedData.binnedExpMap
		self.expCube=LIKE.binnedData.expCube
		self.IRFs=LIKE.binnedData.irfs
		self.source=SrcName
		self.bandModel=model
		self.NBins=nbins
		self.obsROI=self.bdAn.binnedData.observation.roiCuts().roiCone()
		self.obsEBounds=self.bdAn.binnedData.observation.roiCuts().getEnergyCuts()
		self.pc=phCorr
	
	#print out the parameters
	def Print(self):
		print self.bdAn
		print 'Spacecraft file:',self.ft2
		print 'Event file:',self.ft1
		print 'Model for energy bands:', self.bandModel
		print 'Source of Interest:', self.source
	
	#makes bins for plotting taking max energy into account, makes NBins bands for full rnge and then accesses plotEbounds function
	def plotBins(self,MaxE):
		ebounds=log_array(self.NBins+1,self.obsEBounds[0],self.obsEBounds[1])
		self.bins=plotEbounds(ebounds,MaxE)
		self.nbins=len(self.bins[0])
		self.makeFiles()
		self.custBins=0
		
	#function for those who want to define their own plot bins in a manner different from the plotBins function
	def customBins(self,MIN,MAX):#here MIN and MAX are the bin low and high edges (respectively) as you've defined them
		self.bins=[MIN,MAX]
		self.nbins=len(MIN)
		self.makeFiles()
		self.custBins=1
	
	#runs both makeFT1Bands and makeExpMaps
	def makeFiles(self):
		#makeFT1Bands(self,evclsmin,evclsmax)
		makeEFile(self)
		makeCCube(self)
		makeSMap(self)
	
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
		self.SR=getSR(self)
		self.FitBounds=[Emin,Emax]
		self.Fit=getFit(self,Emin,Emax,51,1,expCorrect=expCorrect,wx=writeXML)
			
	#run the full fit on LIKE, objective to get model and bowtie plot info
	def fullFit(self,ftol=1e-3,CoVar=False,toltype='ABS',Emin=0,Emax=0,writeXML=False): #must say CoVar=True if one wants bowtie info
		if(toltype!='ABS' and toltype!='REL'):
			print 'Invalid toltype parameter, valid choices are "REL" for relative tolerance and "ABS" for absolute tolerance.'
			return
		if(toltype=='ABS'):
			self.bdAn.setFitTolType(1)
		else:
			self.bdAn.setFitTolType(0)
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
	
class extbdlikeSED:
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
	
	#for the next two functions, leave ecent as an optional input in the event that the centers are determined via a different method than
	#in the getECent function
	
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
			self.likeIn.bdAn()
			residualsPlot(self,plot)
		except:
			print 'Full energy range BinnedAnalysis object has been deleted (probably with likeInput.delFullLike() function), can not make residuals plot.'
			return

from GtApp import GtApp
import numpy as num
import pyfits
import os
from BinnedAnalysis import *
from ROOT import TCanvas, TH1F, gStyle, TLegend, TFile, TArrow, TLine, TPostScript, TF1, SetOwnership, TGraph, TGraphAsymmErrors
from UpperLimits import UpperLimits
import pyLikelihood as pyLike
from math import sin,cos,log10,pi,ceil,floor
print 'This is extbdlikeSED version 12.'
print 'Adapted for STs v9r23p1 and above.'
#function to get sr extension of extended source for use in fitBands and plotSED functions
def getSR(likeIn):
	pref=getParamIndxnowarn(likeIn.bdAn,likeIn.source,'Prefactor')
	if pref!=-1:
		E0=getParamIndxnowarn(likeIn.bdAn,likeIn.source,'Scale')
		Gamma=likeIn.bdAn[getParamIndxnowarn(likeIn.bdAn,likeIn.source,'Index')]
		Gam=(Gamma.value() if Gamma.getScale()<0 else -1.*Gamma.value())
		PREF=likeIn.bdAn.flux(likeIn.source,likeIn.obsEBounds[0],likeIn.obsEBounds[1])*E0.getValue**(-Gam)*(Gam-1.)/(likeIn.obsEBounds[0]**(1-Gam)-likeIn.obsEBounds[1]**(1-Gam))
		sr=likeIn.bdAn[pref].value()*likeIn.bdAn[pref].parameter.getScale()/PREF
	else:
		Integral=getParamIndx(likeIn.bdAn,likeIn.source,'Integral')
		Upper=getParamIndx(likeIn.bdAn,likeIn.source,'UpperLimit')
		Lower=getParamIndx(likeIn.bdAn,likeIn.source,'LowerLimit')
		INTEGRAL=likeIn.bdAn.flux(likeIn.source,likeIn.bdAn[Lower].value()*likeIn.bdAn[Lower].parameter.getScale(),likeIn.bdAn[Upper].value()*likeIn.bdAn[Upper].parameter.getScale())
		#print INTEGRAL
		#print likeIn.bdAn[Integral].value()*likeIn.bdAn[Integral].parameter.getScale()
		sr=INTEGRAL/(likeIn.bdAn[Integral].value()*likeIn.bdAn[Integral].parameter.getScale())
	return sr
#this function creates an array with npts-1 logarithmically spaced bins
#borrowed from macro Jim sent to do profile likelihood
def log_array(npts, xmin, xmax):
	xstep = num.log(xmax/xmin)/(npts - 1)
	return xmin*num.exp(num.arange(npts, dtype=num.float)*xstep)

#this function creates the different energy band ft1 files from the original ft1
#def makeFT1Bands(likeIn,evmin,evmax):
	#gtselect=GtApp('gtselect') #note, checks for existence of ft1 bands, if found won't remake 
	#i=0
	#while i<likeIn.nbins:
		#out='%s_%ibins_band%i.fits' %(likeIn.source.replace(' ','_'),likeIn.NBins,i) #ft1 has name of source and number of bins to be specific about what was done
		#if(not os.access(out,os.F_OK)):
			#gtselect.run(infile=likeIn.ft1,outfile=out,ra=likeIn.obsROI[0],dec=likeIn.obsROI[1],rad=likeIn.obsROI[2],tmin=0,tmax=0,emin=likeIn.bins[0][i],emax=likeIn.bins[1][i],evclsmin=evmin,evclsmax=evmax,zmax=0,chatter=0) #zmax=0 makes no zenith selection
		#else:
			#print '  -%s' %out,'already exists, skipping gtselect for this energy band.-'
		#i+=1
	#return

def makeEFile(likeIn):
	gtbindef=GtApp('gtbindef')
	mybins=likeIn.bins
	Energies=[]
	#construct the set of energy bins we want making sure that there are 10/decade
	for i in range(0,len(mybins[0])):
		mymin=mybins[0][i]
		mymax=mybins[1][i]
		Nene=int(ceil(10*(log10(mymax)-log10(mymin))))
		Energies+=[l for l in log_array(Nene+1,mymin,mymax)]#have to be tricky as log_array returns an array which can't be added to a list
		if i!=len(mybins[0])-1:#need to avoid doubling up on values
			Energies=Energies[:-1]
	#now we write an ascii file with the bin definitions
	myname='%s_binDefs.ascii'%likeIn.source.replace(' ','_')
	myfile=open(myname,'w')#set to overwrite
	for i in range(0,len(Energies)-1):
		if i!=len(Energies)-2:
			myfile.write('%.4e\t%.4e\n'%(Energies[i],Energies[i+1]))
		else:
			myfile.write('%.4e\t%.4e'%(Energies[i],Energies[i+1]))
	myfile.close()
	gtbindef.run(bintype="E",binfile=myname,outfile=myname.replace('.ascii','.fits'),energyunits="MeV",chatter=0)#note, by default gtbindef is also set to overwrite existing files
	likeIn.ebf=myname.replace('.ascii','.fits')
	return

#this function will calcluate the exposure maps for each energy band, one exposure cube (same as that used for the full fit) is sufficient
def makeCCube(likeIn):
	CC=pyfits.open(likeIn.cc)
	gtbin=GtApp('gtbin')
	Src=likeIn.source
	NBins=likeIn.NBins
	ebins=likeIn.bins
	ft2=likeIn.ft2
	obsROI=likeIn.obsROI
	Proj=CC[0].header['CTYPE1'][-3:]
	likeIn.proj=Proj
	coord=('GAL' if CC[0].header['CTYPE1'][0]=='G' else 'CEL')
	likeIn.coord=coord
	coordX=CC[0].header['CRVAL1']
	coordY=CC[0].header['CRVAL2']
	likeIn.cX=coordX
	likeIn.cY=coordY
	out='%s_%ibins_cc.fits' %(Src.replace(' ','_'),NBins)
	likeIn.binsz=abs(CC[0].header['CDELT1'])
	if(not os.access(out,os.F_OK)):
		gtbin.run(evfile=likeIn.ft1,scfile=ft2,outfile=out,algorithm='CCUBE',ebinalg='FILE',ebinfile=likeIn.ebf,coordsys=coord,xref=coordX,yref=coordY,nxpix=CC[0].header['NAXIS1'],nypix=CC[0].header['NAXIS2'],binsz=abs(CC[0].header['CDELT1']),proj=Proj,axisrot=0,chatter=0)
	else:
		print '%s already exists, skipping gtbin.' %out
	CC.close()
	return

#make source maps for each energy band
def makeSMap(likeIn):
	#gtexpcube2=GtApp('gtexpcube2')
	gtsrcmaps=GtApp('gtsrcmaps')
	Src=likeIn.source
	NBins=likeIn.NBins
	IRFs=likeIn.IRFs
	ft2=likeIn.ft2
	ft1=likeIn.ft1
	expCube=likeIn.expCube
	bandModel=likeIn.bandModel
	mybins=likeIn.bins
	out='%s_%ibins_sm.fits' %(Src.replace(' ','_'),NBins)
	cc='%s_%ibins_cc.fits' %(Src.replace(' ','_'),NBins)
	em='%s_%ibins_%s_bdem.fits' %(Src.replace(' ','_'),NBins,IRFs)
	ROIrad=getRad(likeIn.ft1)#get the ROI info, namely the extraction radius
	mypix=(ROIrad+10)*2./likeIn.binsz#need to make the exposure map this big to avoid an error with gtsrcmaps
	if not os.access(em,os.F_OK):
		os.system('gtexpcube2 ebinalg=FILE infile=%s cmap=none outfile=%s irfs=%s nxpix=%i nypix=%i binsz=%f coordsys=%s xref=%f yref=%f axisrot=0 proj=%s ebinfile=%s chatter=0'%(expCube,em,IRFs,mypix,mypix,likeIn.binsz,likeIn.coord,likeIn.cX,likeIn.cY,likeIn.proj,likeIn.ebf))#for some reason the GtApp version wouldn't run with ebinalg="FILE"
		#gtexpcube2.run(infile=expCube,cmap='none',outfile=em,irfs=IRFs,nxpix=mypix,nypix=mypix,binsz=likeIn.binsz,coordsys=likeIn.coord,xref=likeIn.cX,yref=likeIn.cY,axisrot=0,proj=likeIn.proj,ebinalg='FILE',ebinfile=likeIn.ebf,chatter=4)
		#gtexpcube2.run(infile=expCube,cmap=cc,outfile=em,irfs=IRFs,ebinalg='LOG',emin=mybins[0][0],emax=mybins[1][-1],chatter=4)
	else:
		print '%s already exists, skipping gtexpcube2.'%em
	if not os.access(out,os.F_OK):
		gtsrcmaps.run(scfile=ft2,expcube=expCube,cmap=cc,srcmdl=bandModel,bexpmap=em,outfile=out,irfs=IRFs,chatter=0)
	else:
		print '%s already exists, skipping gtsrcmaps.' %out
	return

#creates lower and upper bounds for energy bins, takes into account maximum energy found wihtin 95% containment of the source
def plotEbounds(energies,maxE):
	i=len(energies)-1
	j=0
	while i>=0: #start at the top, check to see which energy bin contains maxE, corresponds to highest energy bin for plotting
		if(energies[i]>=maxE and energies[i-1]<=maxE):
			j=i
			i=-1
		i-=1
	mins=energies[0:j]
	maxs=energies[1:j+1]
	return mins,maxs

#this function takes in the lower and upper limits on the energy bands and produces a weighted center for each bin, returns in units of GeV
#def centEnergy(MIN,MAX):
	#cent=[]
	#for m,M in zip(MIN,MAX): #get the energy center by finding the mid point of the log10 of the E bounds
		#cent+=[10**(0.5*(num.log10(m)+num.log10(M)))] #this particular way of getting the bin center was taken from on old IDL macro Max sent me
	#GeVcent=[]
	#for x in cent:
		#GeVcent+=[float(x/1000.)]
	
	#return GeVcent

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
	flux=likeIn.bdAn.flux
	nbins=likeIn.nbins
	pc=likeIn.pc
	SR=likeIn.SR
	mybins=likeIn.bins
	sm='%s_%ibins_sm.fits' %(Src.replace(' ','_'),NBins)
	em='%s_%ibins_%s_bdem.fits' %(Src.replace(' ','_'),NBins,IRFs)
	band_obs=BinnedObs(sm,expCube,em,IRFs)
	for i in range(0,nbins): #first, set up for a likelihood run
		print '  -Runnng Likelihood for band%i-' %i
		#band_obs=BinnedObs(sm,expCube,em,IRFs)
		band_like=BinnedAnalysis(band_obs,bandModel,opt)
		emin=mybins[0][i]
		emax=mybins[1][i]
		band_like.setEnergyRange(emin+1,emax-1)#should set the analysis to only do between these two energies, +/- 1 to make sure we get the proper range, avoids rounding issues
		band_like.setFitTolType(ttype)
		stype=band_like.model.srcs[Src].spectrum().genericName()
		for src in band_like.sourceNames():
			par=band_like.normPar(src)
			par.setValue(par.getValue()*pc)
		#then set the scale factor to the center of the energy band, make sure it's frozen, and get index for prefactor while you're at it
		if stype=='PowerLaw':
			scale=getParamIndx(band_like,Src,'Scale') #this is where you have to use PowerLaw, for PowerLaw2 these parameters don't exist and this will cause problems
			pref=getParamIndx(band_like,Src,'Prefactor')
			band_like.freeze(scale)
			band_like[scale].setBounds(20,5e5)
			band_like[scale].setScale(1)
			band_like[scale]=1000.*ecent[i] #put center energies in units of GeV but xml files use MeV
			#multiplier=band_like[pref].getScale() #need to get the scale of the prefactor so the values will not be too large
			try:
				logFlux=log10(flux(Src,emin=emin,emax=emax)/(emax-emin))
			except:
				logFlux=-14
			newScale=max(int(floor(logFlux)),-14)
			band_like[pref].setScale(10**newScale)
			multiplier=10**newScale
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
					band_like.freeze(HIGH) #just in case, make sure these aren't fit values
					band_like.freeze(LOW)
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
				band_like.fit(tol=1,verbosity=0,optimizer=opt)
				fail=0			
			except:
				try:
					band_like.fit(tol=1*10,verbosity=0,optimizer=opt)
					fail=0
				except:
					try:
						band_like.fit(tol=1./10,verbosity=0,optimizer=opt)
						fail=0
					except:
						print "Fit with optimizer %s with tolerances ~1 to look for negative or zero TS sources failed, if error bars are unrealistically small you may need to redo the fit for energy band %i manually" %(opt,i)
						fail=1
			if fail==0:
				for src in band_like.sourceNames():
					if src!=Src:
						par=band_like.normPar(src)
						if band_like[src].type=='PointSource' and band_like.Ts(src)<=0 and par.isFree()==True:
							band_like.deleteSource(src)
							print "   -Removing %s from the model" %src
			scale=getParamIndx(band_like,Src,'Scale') 
			pref=getParamIndx(band_like,Src,'Prefactor')
			#do the actual fit
			try:
				band_like.fit(tol=ftol,verbosity=0)
			except:
				try:
					print 'Trying lower tolerance of %s for band%i.' %(ftol/10,i)
					band_like.fit(tol=ftol/10,verbosity=0)
				except:
					try:
						print 'Trying higher tolerance of %s for band%i.' %(ftol*10,i)
						band_like.fit(tol=ftol*10,verbosity=0)
					except:
						print 'No convergence for band%i, skipping.' %i
						pts+=[0]
						errs+=[0]
						tsPts+=[0]
						pass
			#get the prefactor, error, and ts values
			val=band_like[pref].value()
			err=band_like[pref].error()
			TS=band_like.Ts(Src)
			if(TS<tslim or (i==(nbins-1) and lastbinUL==True)): #calculate 95% upperlimit if source TS<tslim, 25 by default (corresponds to ~5sigma)
				try:
					band_like.freeze(pref)
					ul=UpperLimits(band_like)
					UL=ul[Src].compute(emin=emin,emax=emax)
					val=UL[1]
					err=0
					print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
				except:
					try:
						print '   Tyring higher tolerance of %s for band %i to get good starting point for upper limit calculations.' %(ftol*10,i)
						band_like[pref].setFree(1)
						band_like.fit(tol=ftol*10,verbosity=0)
						TS=band_like.Ts(Src)
						if(TS<tslim or (i==(nbins-1) and lastbinUL==True)):
							band_like.freeze(pref)
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
							band_like.fit(tol=ftol/10,verbosity=0)
							TS=band_like.Ts(Src)
							if(TS<tslim or (i==(nbins-1) and lastbinUL==True)):
								band_like.freeze(pref)
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
			pts+=[val*multiplier*SR/pc]
			errs+=[err*multiplier*SR/pc]
		if stype=='PowerLaw2':
			Upper=getParamIndx(band_like,Src,'UpperLimit')
			Lower=getParamIndx(band_like,Src,'LowerLimit')
			Integral=getParamIndx(band_like,Src,'Integral')
			Index=getParamIndx(band_like,Src,'Index')
			band_like.freeze(Upper)
			band_like.freeze(Lower)
			band_like[Upper].setBounds(20,5e5)
			band_like[Lower].setBounds(20,5e5)
			band_like[Lower].setScale(1)
			band_like[Upper].setScale(1)
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
					band_like.freeze(HIGH) #just in case, make sure these aren't fit values
					band_like.freeze(LOW)
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
				band_like.fit(tol=1,verbosity=0,optimizer=opt)
				fail=0			
			except:
				try:
					band_like.fit(tol=1*10,verbosity=0,optimizer=opt)
					fail=0
				except:
					try:
						band_like.fit(tol=1./10,verbosity=0,optimizer=opt)
						fail=0
					except:
						print "Fit with optimizer %s with tolerances ~1 to look for negative or zero TS sources failed, if error bars are unrealistically small you may need to redo the fit for energy band %i manually" %(opt,i)
						fail=1
			if fail==0:
				for src in band_like.sourceNames():
					if src!=Src:
						par=band_like.normPar(src)
						if band_like[src].type=='PointSource' and band_like.Ts(src)<=0 and par.isFree()==True:
							band_like.deleteSource(src)
							print "   -Removing %s from the model" %src
			Upper=getParamIndx(band_like,Src,'UpperLimit')
			Lower=getParamIndx(band_like,Src,'LowerLimit')
			Integral=getParamIndx(band_like,Src,'Integral')
			Index=getParamIndx(band_like,Src,'Index')
			try:
				band_like.fit(tol=ftol,verbosity=0)
			except:
				try:
					print 'Trying lower tolerance of %s for band%i.' %(ftol/10,i)
					band_like.fit(tol=ftol/10,verbosity=0)
				except:
					try:
						print 'Trying higher tolerance of %s for band%i.' %(ftol*10,i)
						band_like.fit(tol=ftol*10,verbosity=0)
					except:
						print 'No convergence for band%i, skipping.' %i
						pts+=[0]
						errs+=[0]
						tsPts+=[0]
						gamma+=[0]
						pass
			val=band_like[Integral].value()
			err=band_like[Integral].error()
			TS=band_like.Ts(Src)
			gam=band_like[Index].value()*indxMult*-1.
			if(TS<tslim  or (i==(nbins-1) and lastbinUL==True)): #calculate 95% upperlimit if source TS<9, i.e. less than 3 sigma detection in each energy band
				try:
					band_like.freeze(Integral)
					ul=UpperLimits(band_like)
					UL=ul[Src].compute(emin=emin,emax=emax)
					val=UL[1]
					err=0
					#need to redo the fit with band_like set to Upper Limit value to get correct spectral index for that value
					band_like[Integral]=val
					band_like.freeze(Integral)
					band_like.fit(tol=ftol,verbosity=0)
					gam=band_like[Index].value()*indxMult*-1.
					print '    NOTE: Band%i,' %i,'with center energy',ecent[i],'GeV, quoting 95% upper limit on flux.'
				except:
					try:
						print '   Tyring higher tolerance of %s for band %i to get good starting point for upper limit calculations.' %(ftol*10,i)
						band_like[Integral].setFree(1)
						band_like.fit(tol=ftol*10,verbosity=0)
						TS=band_like.Ts(Src)
						if(TS<tslim  or (i==(nbins-1) and lastbinUL==True)):
							band_like.freeze(Integral)
							ul=UpperLimits(band_like)
							UL=ul[Src].compute(emin=emin,emax=emax)
							val=UL[1]
							err=0
							band_like[Integral]=val
							band_like.freeze(Integral)
							band_like.fit(tol=ftol*10,verbosity=0)
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
							band_like.fit(tol=ftol/10,verbosity=0)
							TS=band_like.Ts(Src)
							if(TS<tslim  or (i==(nbins-1) and lastbinUL==True)):
								band_like.freeze(Integral)
								ul=UpperLimits(band_like)
								UL=ul[Src].compute(emin=emin,emax=emax)
								val=UL[1]
								err=0
								band_like[Integral]=val
								band_like.freeze(Integral)
								band_like.fit(tol=ftol/10,verbosity=0)
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
			pts+=[val*multiplier*SR/pc]
			errs+=[err*multiplier*SR/pc]
			gamma+=[gam]
			if wx:
				band_like.writeXml('%s_%ibins_band%i_fitmodel.xml'%(Src.replace(' ','_'),NBins,i))
			del band_like
			#del band_obs
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
	
def getParamIndxnowarn(fit,name,NAME):
	ID=-1
	spec=fit[name].funcs['Spectrum']
	for indx, parName in zip(spec._parIds, spec.paramNames):
		if(parName==NAME):
			ID = indx
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
	modEs=num.array(modelEcent)
	modFs=num.array(sed.Fit)
	modesq=[y*x**2*0.001602 for x,y in zip(modelEcent,sed.Fit)]
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
	GeVBands=[float(x/1000.) for x in minEs]
	GeVBands+=[float(maxEs[-1]/1000.)]
	Bands=num.array(GeVBands) #make an array that TH1F will recognize for binning purposes
	#if using PowerLaw 2 model, need to make a few adjustments to the data points and errors
	if len(sed.data[3])==0:
		GeVflxPts=[float(x*1000.) for x in fluxPts if x!=0]#convert to units of cm^-2 s^-1 GeV^-1 for points and errors
		GeVfluxErrs=[float(x*1000.) for x,X in zip(fluxErrs,fluxPts) if X!=0]
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
	flxArrows=[]
	newecent=[]
	for e,f in zip(ecent,fluxPts):
		if f!=0:
			newecent+=[e]
	energies=num.array(newecent)
	esqf=[]
	esqfErr=[]
	if(len(sed.data[3])==0):
		for x,y,z in zip(energies,GeVfluxPts,GeVfluxErrs):
			if y!=0:
				esqf+=[y*x**2*0.001602]
				esqfErr+=[z*x**2*0.001602]
	else:
		for x,y,z,g,e,E in zip(energies,fluxPts,fluxErrs,gammas,minEs,maxEs):
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
			flxArrows+=[TArrow(z,x,z,x*0.50,0.025,'|->')]	
	
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
		hold=[]
		for x,z in zip(sed.Fit,modelEcent):
			hold+=[x*z**2*0.001602]
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
	writeSpecFits(sed,modBins,modelEcent,[x/1000. for x in newminEs],[x/1000. for x in newmaxEs],GeVfluxPts,GeVfluxErrs,ecent,EsqFluxes,EsqFErrors) #write the output to a sed.Fits file
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
	hdulist.writeto('%s_%ibins_extbdlikeSEDout.fits' %(sed.likeIn.source.replace(' ','_'),sed.likeIn.NBins),output_verify='ignore',clobber='yes') #write the file
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
	mysrc=pyLike.DiffuseSource_downcastAsDiffuse(sed.likeIn.bdAn[sed.likeIn.source].src)
	for x in ecent: #using the same center energies as used in the data
		MeV=x*1000. #ecent is in GeV, but mysrc.spectrum assumes MeV
		arg=pyLike.dArg(MeV)
		val=mysrc.spectrum()(arg)
		modelPts+=[float(val*sed.likeIn.SR)] 
	#allow for the possibility of the model points have errors too
	if(sed.likeIn.bdAn.covariance is None):
		modelErrs=[0]*len(modelPts)
	else:
		covArray=num.array(sed.likeIn.bdAn.covariance)
		srcCovArray=[]
		par_index_map={}
		indx=0
		for src in sed.likeIn.bdAn.sourceNames():
			parNames=pyLike.StringVector()
			sed.likeIn.bdAn[src].src.spectrum().getFreeParamNames(parNames)
			for par in parNames:
				par_index_map['::'.join((src,par))]=indx
				indx +=1
		srcPars=pyLike.StringVector()
		sed.likeIn.bdAn[sed.likeIn.source].src.spectrum().getFreeParamNames(srcPars)
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
			modelErrs+=[float(val*sed.likeIn.SR)]
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
	centEs=[]
	for i in range(0,len(modEs)-1):
		centEs+=[0.5*(modEs[i]+modEs[i+1])]
	#correct for difference in exposure
	for src in likeIn.bdAn.sourceNames():
		par=likeIn.bdAn.normPar(src)
		par.setValue(par.getValue()*likeIn.pc)
	
	#do the full fit
	try:
		likeIn.bdAn.fit(tol=FTol,covar=CoVar,verbosity=0)
	except:
		try:
			print 'Trying lower tolerance of %s for full fit.' %FTol/10.
			likeIn.bdAn.fit(tol=FTol/10.,covar=CoVar,verbosity=0)
		except:
			try:
				print 'Trying higher tolerance of %s for full fit.' %FTol*10.
				likeIn.bdAn.fit(tol=FTol*10.,covar=CoVar,verbosity=0)
			except:
				print 'Error, no convergence.'
				return
			
	
	#most of the following (getting model and bowtie) is taken directly from David Sanchez's pyUnfoldPlot with some minor stylistic changes
	
	likeIn.SR=getSR(likeIn)
	
	#get the model
	mysrc=pyLike.DiffuseSource_downcastAsDiffuse(likeIn.bdAn[likeIn.source].src)
	spec=[]
	for x in centEs:
		arg=pyLike.dArg(x)
		val=mysrc.spectrum()(arg) #gives (I believe) dN/dE spectrum in ph/cm^2/s/MeV so need to convert to ph/cm^2/s/GeV
		spec+=[float(1000.*val*likeIn.SR/likeIn.pc)] #correct to full exposure value
	
	if(likeIn.bdAn.covariance is None):
		bt=[0]
	
	else:
		bt=[]
		covArray=num.array(likeIn.bdAn.covariance)
		srcCovArray=[]
		par_index_map={}
		indx=0
		for src in likeIn.bdAn.sourceNames():
			parNames=pyLike.StringVector()
			likeIn.bdAn[src].src.spectrum().getFreeParamNames(parNames)
			for par in parNames:
				par_index_map['::'.join((src,par))]=indx
				indx +=1
		srcPars=pyLike.StringVector()
		likeIn.bdAn[likeIn.source].src.spectrum().getFreeParamNames(srcPars)
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
			bt+=[float(1000.*val*likeIn.SR/likeIn.pc)] #correct to full exposure value
	myfile=open('extbdlikeSED_%s_fullFitout.txt'%likeIn.source.replace(' ','_'),'w')
	print 'Full energy range model for %s:' %likeIn.source
	myfile.write('Full energy range model for %s:\n' %likeIn.source)
	print likeIn.bdAn[likeIn.source]
	myfile.write('%s\n'%likeIn.bdAn[likeIn.source])
	if likeIn.bdAn.covariance is None:
		myfile.write('(Covariance Matrix not calculated)\n')
		print 'Flux %.1f-%.1f GeV %.1e cm^-2 s^-1' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]))
		myfile.write('Flux %.1f-%.1f GeV %.1e cm^-2 s^-1\n' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1])))
	else:
		print 'Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]),likeIn.bdAn.fluxError(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]))
		myfile.write('Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1\n' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]),likeIn.bdAn.fluxError(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1])))
	print 'Test Statistic',likeIn.bdAn.Ts(likeIn.source)
	myfile.write('Test Statistic %.2f'%likeIn.bdAn.Ts(likeIn.source))
	myfile.close()
	if wx:
		likeIn.bdAn.writeXml('%s_fullErange_fitmodel.xml'%likeIn.source.replace(' ','_'))
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
	
	if expCorrect:
		for src in likeIn.bdAn.sourceNames():
			par=bdAn.normPar(src)
			par.setValue(par.getValue()/likeIn.pc)
			par.setError(par.error()/likeIn.pc)
	
	#most of the following (getting model and bowtie) is taken directly from David Sanchez's pyUnfoldPlot with some minor stylistic changes
	
	#get the model
	mysrc=pyLike.DiffuseSource_downcastAsDiffuse(likeIn.bdAn[likeIn.source].src)
	spec=[]
	for x in centEs:
		arg=pyLike.dArg(x)
		val=mysrc.spectrum()(arg) #gives (I believe) dN/dE spectrum in ph/cm^2/s/MeV so need to convert to ph/cm^2/s/GeV
		spec+=[float(1000.*val*likeIn.SR)]
	
	if(likeIn.bdAn.covariance is None):
		bt=[0]
	
	else:
		bt=[]
		covArray=num.array(likeIn.bdAn.covariance)
		srcCovArray=[]
		par_index_map={}
		indx=0
		for src in likeIn.bdAn.sourceNames():
			parNames=pyLike.StringVector()
			likeIn.bdAn[src].src.spectrum().getFreeParamNames(parNames)
			for par in parNames:
				par_index_map['::'.join((src,par))]=indx
				indx +=1
		srcPars=pyLike.StringVector()
		likeIn.bdAn[likeIn.source].src.spectrum().getFreeParamNames(srcPars)
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
			bt+=[float(1000.*val*likeIn.SR)]
	if prtMod==1:
		myfile=open('extbdlikeSED_%s_fullFitout.txt'%likeIn.source.replace(' ','_'),'w')
		print 'Full energy range model for %s:' %likeIn.source
		myfile.write('Full energy range model for %s:\n' %likeIn.source)
		print likeIn.bdAn[likeIn.source]
		myfile.write('%s\n'%likeIn.bdAn[likeIn.source])
		if likeIn.bdAn.covariance is None:
			myfile.write('(Covariance Matrix not calculated)\n')
			print 'Flux %.1f-%.1f GeV %.1e cm^-2 s^-1' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]))
			myfile.write('Flux %.1f-%.1f GeV %.1e cm^-2 s^-1\n' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1])))
		else:
			print 'Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]),likeIn.bdAn.fluxError(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]))
			myfile.write('Flux %.1f-%.1f GeV %.1e +/- %.1e cm^-2 s^-1\n' %(likeIn.obsEBounds[0]/1000.,likeIn.obsEBounds[1]/1000.,likeIn.bdAn.flux(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1]),likeIn.bdAn.fluxError(likeIn.source,emin=likeIn.obsEBounds[0],emax=likeIn.obsEBounds[1])))
		print 'Test Statistic',likeIn.bdAn.Ts(likeIn.source)
		myfile.write('Test Statistic %.2f',likeIn.bdAn.Ts(likeIn.source))
		myfile.close()
		if wx:
			likeIn.bdAn.writeXml('%s_fullErange_fitmodel.xml'%likeIn.source.replace(' ','_'))
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

def getRad(ft1):
	file=pyfits.open(ft1)
	num=file[1].header['NDSKEYS']
	header=file[1].header
	right='POS(RA,DEC)'
	i=1
	keynum=0
	while i<=num:  #this step is necessary since it is not clear that the POS key word will have the same number always
		word='DSTYP%i' %i
		test=file[1].header[word]
		if(test==right):
			keynum=i
			i=num
		i+=1
	if(keynum==0):  #DSKEYS start numbering at 1, if this value hasn't been updated, KEYword doesn't exist
		print 'Error: No position keyword found in fits header (assuming position is RA and DEC.  Exiting...'
		exit()
	keyword='DSVAL%i' %keynum
	ra,dec,rad=header[keyword].strip('CIRCLE()').split(',') #gets rid of the circle and parenthesis part and splits around the comma
	file.close()
	return float(rad)