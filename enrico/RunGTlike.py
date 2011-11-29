
import os,sys,string
import SummedLikelihood
#import fitfunction 
#import gtfunction 
import Utility
from enrico.config import get_config

def run(infile) :
	Configuration = get_config(infile)

	folder = Configuration['out']
	os.system('mkdir -p '+folder)


	TSlimit = float(Configuration['UpperLimit']['TSlimit'])
	SPindex = float(Configuration['Spectrum']['FreeSpectralIndex'])
	SummedLike = Configuration['Spectrum']['SummedLike']

	if SummedLike == 'yes' :
		### create an instance of Observation for the FRONT events
		runfitfront, Obs_FRONT = Utility.Analysis(folder,Configuration,tag="FRONT",convtyp=0)
	
		### create an instance of Observation for the BACK events
		runfitback, Obs_BACK = Utility.Analysis(folder,Configuration,tag="BACK",convtyp=1)

		FitB = runfitback.CreateFit()
		FitF = runfitfront.CreateFit()
	
		Fit = SummedLikelihood.SummedLikelihood()
		Fit.addComponent(FitB)
		Fit.addComponent(FitF)
		runfit = runfitback

	else :
	
		### create an instance of Observation 
		runfit,Obs = Utility.Analysis(folder,Configuration,tag="",convtyp=Configuration['analysis']['convtype'])
		Fit = runfit.CreateFit()

	Result = runfit.PerformFit(Fit)
	Utility.DumpResult(Result,Configuration)

	if Configuration['Spectrum']['ResultPlots'] == 'yes' :
		runfit.PlotSED(Fit)
		outXml = folder+"/"+runfit.Observation.srcname+"_"+Configuration['file']['tag']+"_out.xml"
		if SummedLike == 'yes' :
			runfitback.ModelMap(outXml)
			runfitfront.ModelMap(outXml)
		else :
			runfit.ModelMap(outXml)

	if int(Configuration['Ebin']['NumEnergyBins']) > 0 :
		configfiles = Utility.PrepareEbin(Fit,runfit)
	
		for conf in configfiles:
			os.system('enrico_submit '+conf)


if __name__=='__main__':

	try :
		infile = sys.argv[1]
	except :
		print('FATAL: Config file not found.')
		sys.exit(1)

	run(infile)

