import Utility
import os
import string
import numpy,array
import ROOT
import Style
from enrico.RunGTlike import run
import enrico.pyPlot as pyPlot
import enrico.environ as environ
from enrico.config import get_config
from enrico.submit import call

def PrepareLC(infile) :

	Configuration = get_config(infile)

	Tag = Configuration['file']['tag']
	Nbin = Configuration['LightCurve']['NLCbin']
	tmin = Configuration['time']['tmin']
	tmax = Configuration['time']['tmax']

	Configuration['UpperLimit']['TSlimit'] = Configuration['LightCurve']['TSLightCurve']
	Configuration['out'] +='/LightCurve'

	Configuration['Spectrum']['ResultPlots'] = 'no' 
	Configuration['Ebin']['NumEnergyBins'] = 0
	Configuration['UpperLimit']['envelope'] ='no' 

	AllConfigFile = []
	for i in xrange(Nbin):
		dt = (tmax-tmin)/Nbin
		Configuration['time']['tmin'] = tmin + i*dt
		Configuration['time']['tmax'] = tmin + (i+1)*dt
		Configuration['file']['tag'] = Tag+'_LC_'+str(i)
		filename = Configuration['out']+"/Config_"+str(Configuration['time']['tmin'])+"_"+str(Configuration['time']['tmax'])

		if Configuration['LightCurve']['Prepare'] :
			print "writting the config file ",filename
			Configuration.write(open(filename, 'w'))

		AllConfigFile.append(filename)
	return AllConfigFile


def MakeLC(infile) :
	ROOT.gROOT.SetBatch(ROOT.kTRUE) 
	Style.RootStyle()

	Doplot = True;

	enricodir = environ.DIRS[7][1]

	Configuration = get_config(infile)

	folder = Configuration['out']
	os.system('mkdir -p '+folder+'/LightCurve')

	AllConfigFile = PrepareLC(infile)

	Nbin = Configuration['LightCurve']['NLCbin']

	Time = numpy.array(Nbin*[0.])
	TimeErr = numpy.array(Nbin*[0.])
	Flux = numpy.array(Nbin*[0.])
	FluxErr = numpy.array(Nbin*[0.])

	Npred = numpy.array(Nbin*[0.])
	TS = numpy.array(Nbin*[0.])
	Configuration['Spectrum']['FitsGeneration'] == Configuration['LightCurve']['FitsGeneration'] 

	if Configuration['LightCurve']['FitsGeneration'] == 'yes' or Configuration['LightCurve']['Re-Fit'] == 'yes':
		for i in xrange(Nbin):

			if Configuration['LightCurve']['Submit'] == 'yes':
				Doplot = False
				cmd = "enrico_fit "+AllConfigFile[i]

				scriptname=folder+"/LightCurve/LC_Script_"+str(i)+".sh"
				JobLog = folder+"/LightCurve/LC_Job_"+str(i)+".log"
				JobName = Configuration['target']['name']+"_"+str(i)+".log"

				call(cmd,enricodir,scriptname,JobLog,JobName)
#				os.system("(source "+enricodir+"/init.sh );(sleep 1);(enrico_submit "+AllConfigFile[i]+")")
			else :
				run(AllConfigFile[i])
#				os.system("(source "+enricodir+"/init.sh );(sleep 1);(enrico_fit "+AllConfigFile[i]+")")

	if Configuration['LightCurve']['Plot'] == 'yes' and Doplot:
		for i in xrange(Nbin):
			print "Reading "+AllConfigFile[i]
			CurConfig = get_config(AllConfigFile[i])
			ResultDic = Utility.ReadResult(CurConfig)
			
			Time[i] = (ResultDic.get("tmax")+ResultDic.get("tmin"))/2.
			TimeErr[i] =(ResultDic.get("tmax")-ResultDic.get("tmin"))/2.
			if ResultDic.has_key('Ulvalue') :
				Flux[i] = ResultDic.get("Ulvalue")
				FluxErr[i] = 0
			else :
				Flux[i] = ResultDic.get("Flux")
				FluxErr[i] = ResultDic.get("dFlux")

			Npred[i] = ResultDic.get("Npred")
			TS[i] = ResultDic.get("TS")

		if Configuration['LightCurve']['DiagnosticPlots'] == 'yes':
			gTHNpred,TgrNpred = pyPlot.PlotNpred(Npred,Flux,FluxErr)
			CanvNpred = ROOT.TCanvas()
			gTHNpred.Draw()
			TgrNpred.Draw('zP')
			CanvNpred.Print(folder+'/LightCurve/'+ Configuration['target']['name']+'_Npred.eps')
			CanvNpred.Print(folder+'/LightCurve/'+ Configuration['target']['name']+'_Npred.C')

			gTHTS,TgrTS = pyPlot.PlotTS(Time,TimeErr,TS)
			CanvTS = ROOT.TCanvas()
			gTHTS.Draw()
			TgrTS.Draw('zP')
			CanvTS.Print(folder+'/LightCurve/'+ Configuration['target']['name']+'_TS.eps')
			CanvTS.Print(folder+'/LightCurve/'+ Configuration['target']['name']+'_TS.C')



		gTHLC,TgrLC,ArrowLC = pyPlot.PlotLC(Time,TimeErr,Flux,FluxErr)
		CanvLC = ROOT.TCanvas()
		gTHLC.Draw()
		TgrLC.Draw('zP')

		for i in xrange(len(ArrowLC)):
			ArrowLC[i].Draw()

		CanvLC.Print(folder+'/LightCurve/'+ Configuration['target']['name']+'_LC.eps')
		CanvLC.Print(folder+'/LightCurve/'+ Configuration['target']['name']+'_LC.C')

