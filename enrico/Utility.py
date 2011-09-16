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


#def ReadParam(filename) :
#	FileAna = open(filename,"r")
#	lines =  FileAna.readlines()
#	FileAna.close()
#	dictionary = {}
#	i=0
#	for line in lines:
#		dictionary.update({string.split(lines[i])[0]:string.split(lines[i])[1]})		
#		i+=1
#	return dictionary

#def SaveParam(filename,dictionary) :
#	FileAna = open(filename,"w")

#	for key in dictionary.iterkeys():
#    		FileAna.write(key+' '+str(dictionary[key])+'\n')
#	FileAna.close()

def ChangeModel(inFit,E0,TSlimit=50) :
	Fit = inFit
	for name in Fit.model.srcNames :
	   if Fit.model.srcs[name].spectrum().genericName() == 'PowerLaw' :
		IdPref = getParamIndx(Fit,name,'Prefactor')
		IdGamma  = getParamIndx(Fit,name,'Index')
		IdScale  = getParamIndx(Fit,name,'Scale')

		Flux    = Fit[IdPref].value()
		ErrFlux = Fit[IdPref].error()
		Scale = Fit[IdPref].getScale()
		Escale  = Fit[IdScale].value()
		Gamma = Fit[IdGamma].value()

		NewFlux  = Flux*pow(E0/Escale,-Gamma)*Scale
		NormFlux = pow(10,log10(NewFlux)-int(log10(NewFlux)))
		NewScale =  pow(10,int(log10(NewFlux)))

		Fit[IdScale].setBounds(0.0,4e5)
		Fit[IdPref].setBounds(0,Flux*1000)	

		Fit[IdPref].setScale(NewScale)
		Fit[IdPref] = NormFlux
		Fit[IdPref].setBounds(NormFlux*0.05,NormFlux*50)
		Fit[IdScale] = E0
		Fit[IdScale].setBounds(E0*0.05,E0*50)
	return Fit


def Analysis(folder,Configuration,tag="",convtyp='-1'):
	Obs = gtfunction.Observation(folder,Configuration,convtyp,tag=tag)
#	Option.update(DicParams)

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
	NEbin = int(runfit.Configuration['enricobehavior']['TSlimit'])
	NewConfig = Config
	NewConfig['enricobehavior']['NumEnergyBins'] = '0'

	lEmax = log10(float(runfit.Configuration['energy']['emax']))
	lEmin = log10(float(runfit.Configuration['energy']['emin']))

	ener = numpy.logspace(lEmin,lEmax,NEbin+1)
	os.system("mkdir -p "+runfit.Configuration['out']+'/Ebin')
	paramsfile = []

	for ibin in xrange(NEbin):
		E = int(pow(10,(log10(ener[ibin+1])+log10(ener[ibin]))/2))
		ChangeModel(Fit,E)
		NewModel = runfit.Configuration['out']+"/Ebin/"+runfit.Configuration['target']['name']+"_"+str(E)+".xml"
		Fit.writeXml(NewModel)
		NewConfig['energy']['emin'] = str(ener[ibin])
		NewConfig['energy']['emax'] = str(ener[ibin+1])
		NewConfig['file']['tag'] = runfit.Option['tag']+'_Ebin'+str(ibin) 
		NewConfig['out'] = runfit.Configuration['out']+'/Ebin'
		NewConfig['enricobehavior']['ResultPlots'] = 'no'
		NewConfig['enricobehavior']['FitsGeneration'] = 'yes'
		NewConfig['enricobehavior']['TSlimit'] = NewConfig['enricobehavior']['TSEnergyBins']


	NewConfig.write(open(paramsfile[ibin], 'w'))
#		paramsfile.append(runfit.Option.get('folder')+"/Ebin/"+runfit.Option.get('srcname')+"_"+str(E)+".conf")
#		SaveParam(paramsfile[ibin],NewConfig)
	return paramsfile

