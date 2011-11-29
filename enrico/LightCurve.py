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

		if Configuration['LightCurve']['Prepare'] == 'yes':
			Configuration.write(open(filename, 'w'))

		AllConfigFile.append(filename)
	return AllConfigFile


def WriteToAscii(Time,TimeErr,Flux,FluxErr,TS,Npred,filename):
	flc = open(filename,'w')
	flc.write("Time (MET) Delta_Time Flux(ph cm-2 s-1) Delta_Flux TS Npred\n")
	for i in xrange(len(Time)):
		flc.write(str(Time[i])+"\t"+str(TimeErr[i])+"\t"+str(Flux[i])+"\t"+str(FluxErr[i])+"\t"+str(TS[i])+"\t"+str(Npred[i])+"\n")
	flc.close()


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

	Time = []# numpy.array(Nbin*[0.])
	TimeErr = []#numpy.array(Nbin*[0.])
	Flux = []#numpy.array(Nbin*[0.])
	FluxErr = []#numpy.array(Nbin*[0.])

	Npred = []#numpy.array(Nbin*[0.])
	TS = []#numpy.array(Nbin*[0.])

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
		print "Reading output files"
		for i in xrange(Nbin):
			CurConfig = get_config(AllConfigFile[i])
			try :
				ResultDic = Utility.ReadResult(CurConfig)
			except :
				print "WARNING : fail reading the configuration file : ",AllConfigFile[i]
				print "Job Number : ",i
				print "Please have a look at this job log file"
				continue
			
			Time.append((ResultDic.get("tmax")+ResultDic.get("tmin"))/2.)
			TimeErr.append((ResultDic.get("tmax")-ResultDic.get("tmin"))/2.)
			if ResultDic.has_key('Ulvalue') :
				Flux.append(ResultDic.get("Ulvalue"))
				FluxErr.append(0)
			else :
				Flux.append(ResultDic.get("Flux"))
				FluxErr.append(ResultDic.get("dFlux"))

			Npred.append(ResultDic.get("Npred"))
			TS.append(ResultDic.get("TS"))

		TS = numpy.array(TS)
		Npred = numpy.array(Npred)
		Time = numpy.array(Time)
		TimeErr = numpy.array(TimeErr)
		Flux = numpy.array(Flux)
		FluxErr = numpy.array(FluxErr)

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

		lcfilename = folder+'/LightCurve/'+ Configuration['target']['name']+"_results.dat"
		print "Write to Ascii file : ",lcfilename
		WriteToAscii(Time,TimeErr,Flux,FluxErr,TS,Npred,lcfilename)
