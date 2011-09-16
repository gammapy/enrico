import gtfunction
import fitfunction
import pyfits
import string
import numpy
import os
from math import log10

def SubstracFits(Map1,Map2):
	print "Substracting : ",Map1," to ",Map2

	arr1 = pyfits.getdata(Map1)
	arr2 = pyfits.getdata(Map2)

	outfile = "Substract_Model_cmap.fits"
	pyfits.writeto(outfile,arr1-arr2)

	outfile2 = "Residual_Model_cmap.fits"

	pyfits.writeto(outfile2,(arr1-arr2)/arr2)

def PrintResult(Fit,Current_Obs):
	print '# *********************************************************************'
	print "# *** Model Result ***\n"
	
	print Fit.model
	print 	
	print "Source Name\tNpred\tTS"
	for src in Fit.model.srcNames :
		print src,"\t%2.3f\t%2.3f"%(Fit.NpredValue(src),Fit.Ts(src))
	print 
	print '# *********************************************************************'
	print 

	print "Values and Errors [Minos errors] for "+Current_Obs.srcname
	print "TS : ",Fit.Ts(Current_Obs.srcname)
	stype = Fit.model.srcs[Current_Obs.srcname].spectrum().genericName()
	if stype == 'PowerLaw2' :
	
		Flux = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Integral').value()
		ErrFlux = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Integral').error()
		Scale = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Integral').getScale()

		try :
			Interror = Fit.minosError(Current_Obs.srcname, 'Integral')
			print " Integral:  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ] %2.0e"%(Flux,ErrFlux,Interror[0],Interror[1],Scale)
		except :
			print " Integral:  %2.2f +/-  %2.2f  %2.0e"%(Flux,ErrFlux,Scale)

	if stype == 'PowerLaw' :
		Flux = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Prefactor').value()
		ErrFlux = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Prefactor').error()
		Scale = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Prefactor').getScale()
		try :
			Interror = Fit.minosError(Current_Obs.srcname, 'Prefactor')
			print " Prefactor:  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ] %2.0e"%(Flux,ErrFlux,Interror[0],Interror[1],Scale)
		except :
			print " Prefactor:  %2.2f +/-  %2.2f %2.0e"%(Flux,ErrFlux,Scale)

	if stype == 'PowerLaw2' or stype == 'PowerLaw' :
		Gamma = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Index').value()
		ErrGamma = Fit[Current_Obs.srcname].funcs['Spectrum'].getParam('Index').error()
		try :
			Gamerror = Fit.minosError(Current_Obs.srcname, 'Index')
			print " Index:  %2.2f +/-  %2.2f [ %2.2f, + %2.2f ]"%(Gamma,ErrGamma,Gamerror[0],Gamerror[1])
		except :
			print " Index:  %2.2f +/-  %2.2f"%(Gamma,ErrGamma)


def RemoveWeakSources(Fit,SourceName = None):
    print '# *********************************************************************'
    print "# *** Remove all the weak sources ***"
    NoWeakSrcLeft = False
    while not(NoWeakSrcLeft) :
	NoWeakSrcLeft = True
	for src in Fit.model.srcNames :
		if Fit.Ts(src)< 1 and not(src==SourceName) : #TODO
			print "delete source : ",src	
			NoWeakSrcLeft = False
			Fit.deleteSource(src)
	if not(NoWeakSrcLeft) :
		print "# *** Re-optimize ***"
		Fit.optimize(0)
	print 
    return Fit

def GetFlux(Fit):
    print "Source Flux : "
    for src in Fit.model.srcNames :
	try :
		print src+" : %2.2e +/-  %2.2e"%(Fit.flux(src),Fit.fluxError(src))
	except :
		pass
    print

def GetCovar(srcname,Fit):
	ptsrc = pyLike.PointSource_cast(Fit[srcname].src)
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
        Fit[srcname].src.spectrum().getFreeParamNames(srcpars)
        pars = ["::".join((srcname, x)) for x in srcpars]
        for xpar in pars:
           	ix = par_index_map[xpar]
           	my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
        print "The covariance matrix is :\n",numpy.array(my_covar)
	print 


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


def ChangeModel(inFit,Em0,Em1,TSlimit=50) :
	E0 = int(pow(10,(log10(Em1)+log10(Em0))/2))

	Fit = inFit
	for name in Fit.model.srcNames :


	   if Fit.model.srcs[name].spectrum().genericName() == 'PowerLaw' :
		IdPref = getParamIndx(Fit,name,'Prefactor')
		IdEScale  = getParamIndx(Fit,name,'Scale')
		Flux    = Fit[IdPref].value()
		Scale = Fit[IdPref].getScale()
		Escale  = Fit[IdEScale].value()

		IdGamma  = getParamIndx(Fit,name,'Index')
		Gamma = Fit[IdGamma].value()

		NewFlux  = Flux*pow(E0/Escale,-Gamma)*Scale
		NormFlux = pow(10,log10(NewFlux)-int(log10(NewFlux)))
		NewScale =  pow(10,int(log10(NewFlux)))

		Fit[IdEScale].setBounds(0.0,4e5)
		Fit[IdPref].setBounds(0,Flux*1000)	

		Fit[IdPref].setScale(NewScale)
		Fit[IdPref] = NormFlux
		Fit[IdPref].setBounds(NormFlux*0.05,NormFlux*50)
		Fit[IdEScale] = E0
		Fit[IdEScale].setBounds(E0*0.05,E0*50)

	   if Fit.model.srcs[name].spectrum().genericName() == 'PowerLaw2' :
		IdInt = getParamIndx(Fit,name,'Integral')
		IdEmin = getParamIndx(Fit,name,'LowerLimit')
		IdEmax = getParamIndx(Fit,name,'UpperLimit')

		IdGamma  = getParamIndx(Fit,name,'Index')
		Gamma = Fit[IdGamma].value()


		Flux    = Fit[IdInt].value()
		Scale = Fit[IdInt].getScale()

		Emin    = Fit[IdEmin].value()
		Emax    = Fit[IdEmax].value()

		D = pow(Em1,Gamma+1)-pow(Em0,Gamma+1)
		N = pow(Emax,Gamma+1)-pow(Emin,Gamma+1)
		print D," ",N


		NewFlux = Flux*D/N*Scale
		NormFlux = pow(10,log10(NewFlux)-int(log10(NewFlux)))
		NewScale =  pow(10,int(log10(NewFlux)))

		Fit[IdInt].setBounds(0,Flux*1000)	

		Fit[IdInt].setScale(NewScale)
		Fit[IdInt] = NormFlux
		Fit[IdInt].setBounds(NormFlux*0.05,NormFlux*50)
		Fit[IdEmin] = Em0
		Fit[IdEmax] = Em1

	return Fit


def Analysis(folder,Configuration,tag="",convtyp='-1'):
	Obs = gtfunction.Observation(folder,Configuration,convtyp,tag=tag)

	print
	print '# *********************************************************************'
	print '# * 0 - SUMMARY '+tag
	print '# *********************************************************************'
	Obs.printSum()

	runfit  = fitfunction.MakeFit(Obs,Configuration)

	if Configuration['enricobehavior']['FitsGeneration'] == 'yes' :
		runfit.PreparFit()
	return runfit, Obs

def PrepareEbin(Fit,runfit):
	NEbin = int(runfit.Configuration['enricobehavior']['NumEnergyBins'])
	NewConfig = runfit.Configuration
	NewConfig['enricobehavior']['NumEnergyBins'] = '0'
	NewConfig['out'] = runfit.Configuration['out']+'/Ebin'
	NewConfig['enricobehavior']['ResultPlots'] = 'no'
	NewConfig['enricobehavior']['FitsGeneration'] = 'yes'
	NewConfig['enricobehavior']['TSlimit'] = NewConfig['enricobehavior']['TSEnergyBins']

	tag = runfit.Configuration['file']['tag']

	lEmax = log10(float(runfit.Configuration['energy']['emax']))
	lEmin = log10(float(runfit.Configuration['energy']['emin']))
	print "Preparing submission of fit into energy bins"
	print "Emin = ",float(runfit.Configuration['energy']['emin'])," Emax = ",float(runfit.Configuration['energy']['emax'])," Nbins = ",NEbin
	ener = numpy.logspace(lEmin,lEmax,NEbin+1)
	os.system("mkdir -p "+runfit.Configuration['out']+'/Ebin')
	paramsfile = []

	for ibin in xrange(NEbin):
		E = int(pow(10,(log10(ener[ibin+1])+log10(ener[ibin]))/2))
		print "Submition # ",ibin," at energy ",E
		ChangeModel(Fit,ener[ibin],ener[ibin+1])

		NewModel = NewConfig['out']+"/"+runfit.Configuration['target']['name']+"_"+str(E)+".xml"
		Fit.writeXml(NewModel)
		NewConfig['energy']['emin'] = str(ener[ibin])
		NewConfig['energy']['emax'] = str(ener[ibin+1])
		NewConfig['file']['tag'] = tag+'_Ebin_'+str(ibin) 

		paramsfile.append(runfit.Configuration['out']+'/'+runfit.Configuration['target']['name']+"_"+str(E)+".conf")
		NewConfig.write(open(paramsfile[ibin], 'w'))

	return paramsfile

