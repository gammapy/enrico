
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
	SPindex = float(Configuration['enricobehavior']['FreeSpectralIndex'])
	SummedLike = Configuration['enricobehavior']['SummedLike']

	if SummedLike == 'yes' :
		### create an instance of Observation for the FRONT events
		runfitfront, Obs_FRONT = Utility.Analysis(folder,Configuration,tag="FRONT",convtyp='Front')
	
		### create an instance of Observation for the BACK events
		runfitback, Obs_BACK = Utility.Analysis(folder,Configuration,tag="BACK",convtyp='Back')

		FitB = runfitback.CreateFit()
		FitF = runfitfront.CreateFit()
	
		Fit = SummedLikelihood.SummedLikelihood()
		Fit.addComponent(FitB)
		Fit.addComponent(FitF)
		runfit = runfitback

	else :
	
		### create an instance of Observation 
		runfit,Obs = Utility.Analysis(folder,Configuration,tag="",convtyp="All")
		Fit = runfit.CreateFit()

	runfit.PerformFit(Fit)


	if Configuration['enricobehavior']['ResultPlots'] == 'yes' :
		runfit.PlotSED(Fit)
		if SummedLike == 'yes' :
			runfitback.ModelMap(folder+"/"+runfit.Observation.srcname+"_out.xml")
			runfitfront.ModelMap(folder+"/"+runfit.Observation.srcname+"_out.xml")
		else :
			runfit.ModelMap(folder+"/"+runfit.Observation.srcname+"_out.xml")

	if int(Configuration['enricobehavior']['NumEnergyBins']) > 0 :
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

