# Author: David Sanchez dsanchez@llr.in2p3.fr
# Begun: 2008-09-15
#developed for FGST

#defintion of CLASS Run
class Params :
	"define some usefull parametere for the plot"
	def __init__(self,srcname,Emin=100,Emax=3e5,extend=False,PlotName="LAT_SED",LineColor=1,AreaColor=2,plotRes=False):
		self.Emin = Emin
		self.Emax = Emax
		self.N = 500
		self.srcname = srcname
		self.extend = extend
		self.PlotName = PlotName
		self.LineColor = LineColor
		self.AreaColor = AreaColor
		self.plotRes = plotRes

#defintion of CLASS Result
class Result :
	"internal class designed to get the results"
	Prefac_UL = 0.
	def __init__(self,like,Params)	:
		self.like = like
		self.ra = like[Params.srcname].funcs['Position'].getParam('RA').value()
		self.dec = like[Params.srcname].funcs['Position'].getParam('DEC').value()
		self.Model = like[Params.srcname].funcs['Spectrum'].genericName()
		try :
			self.TS=like.Ts(Params.srcname)
		except RuntimeError:
			self.TS=-1
        	self.ptsrc = pyLike.PointSource_cast(like[Params.srcname].src)
        	par_index_map = {}
        	indx = 0
        	for src in like.sourceNames():
        	    parNames = pyLike.StringVector()
        	    like[src].src.spectrum().getFreeParamNames(parNames)
        	    for par in parNames:
         	       par_index_map["::".join((src, par))] = indx
         	       indx += 1
       	 	#
       	 	# Build the source-specific covariance matrix.
        	#
        	if like.covariance is None:
            		raise RuntimeError("Covariance matrix has not been computed.")
        	covar = num.array(like.covariance)
        	if len(covar) != len(par_index_map):
            		raise RuntimeError("Covariance matrix size does not match the " +
                               "number of free parameters.")
        	my_covar = []
        	srcpars = pyLike.StringVector()
        	like[Params.srcname].src.spectrum().getFreeParamNames(srcpars)
        	pars = ["::".join((Params.srcname, x)) for x in srcpars]
        	for xpar in pars:
            		ix = par_index_map[xpar]
            		my_covar.append([covar[ix][par_index_map[ypar]] for ypar in pars])
        	self.covar = numpy.array(my_covar)
		self.srcpars = srcpars	


#GLAST import
from UnbinnedAnalysis import * 
#import pyIrfLoader as irf_loader
import pyLikelihood as pyLike
#import IntegralUpperLimit

#Others import
import string
from math import *
from array import array
import ROOT
import numpy
import pyfits
import os,sys
import copy
import bisect

#defintion of Usefull functions




########################################################################################################
########################################################################################################
########################################################################################################
def MakeError(result,Params):
	like = result.like
        estep = num.log(Params.Emax/Params.Emin)/(Params.N - 1)
        energies = Params.Emin*numpy.exp(estep*num.arange(numpy.float(Params.N)))
	err = numpy.array(Params.N*[0.])
	j=0
	for ene in energies:   
		arg = pyLike.dArg(ene)
		partials = numpy.array(len(result.srcpars)*[0.])
		for i in xrange(len(result.srcpars)):
			x = result.srcpars[i]
			partials[i] = result.ptsrc.spectrum().derivByParam(arg, x)
        	err[j] =  numpy.sqrt(numpy.dot(partials, numpy.dot(result.covar, partials)))
		j+=1
	return energies**2*err*1.602e-6

def dNde(result,energy):
        arg = pyLike.dArg(energy)
        return result.ptsrc.spectrum()(arg)

def MakeFlux(result,Params):
	"Compute differential Flux distribution and corresponding energy and return a numpy array"
	E = numpy.logspace(log10(Params.Emin),log10(Params.Emax),Params.N)
	Flux = numpy.array(Params.N*[0.])
	for i in xrange(Params.N):
		Flux[i] = dNde(result,E[i])
	return E,Flux


def MakeFluxInBin(result,Params,Emin_bin,Emax_bin):
	"Compute differential Flux distribution and corresponding energy and return a numpy array"
	E=numpy.logspace(log10(Emin_bin),log10(Emax_bin),Params.N)
	Flux = numpy.array(Params.N*[0.])
	for i in xrange(Params.N):
		Flux[i] = dNde(result,E[i])
	return E,Flux

def MakeSED(result,Params):
	"Compute Spectral energy distribution and corresponding energy and return a numpy array"
	E=numpy.logspace(log10(Params.Emin),log10(Params.Emax),Params.N)
	nuFnu = numpy.array(Params.N*[0.])
	for i in xrange(Params.N):
		nuFnu[i] = E[i]**2*dNde(result,E[i])*1.602e-6
	return E,nuFnu


def MakeSEDInBin(result,Params,Emin_bin,Emax_bin):
	"Compute Spectral energy distribution and corresponding energy and return a numpy array"
	E=numpy.logspace(log10(Emin_bin),log10(Emax_bin),Params.N)
	nuFnu = numpy.array(Params.N*[0.])
	for i in xrange(Params.N):
		nuFnu[i] = E[i]**2*dNde(result,E[i])*1.602e-6
	return E,nuFnu	




def CountsPlot(Result,Parameter):
	
	ROOT.gROOT.SetBatch(ROOT.kTRUE) 
	imName = "tmp.fits"
	Result.like.writeCountsSpectra(imName)
	image = pyfits.open(imName)

	#loop on the source names to find the good one
	j = 0
	for ID in image[1].data.names:
		if ID==Parameter.srcname :
			indice = j
		j+=1
	emax = image[3].data.field(0)
	emin = image[3].data.field(1)
	E = array('f',((emax+emin)/2.))
	err_E = array('f',((-emax+emin)/2.))
	src = array('f',(image[1].data.field(indice)))
	Nbin = len(src)

	obs = array('f',(image[1].data.field(0)))
	obs_err = array('f',numpy.sqrt(image[1].data.field(0)))
	summ = 0
	for i in xrange(len(image[1].data.names)-1):
		summ = summ +image[1].data.field(i+1)
	other = other = array('f',summ - image[1].data.field(indice))
	summ = array('f', summ)

	residual = array('f', Nbin*[0.])
	Dres = array('f', Nbin*[0.])
	zero =  array('f', Nbin*[0.])
	for i in xrange(Nbin):
	   try :
		residual[i] = (obs[i]-summ[i])/summ[i]
	   except :
		residual[i] = 0.
	   try :
		Dres[i] = abs(residual[i]*obs_err[i]/obs[i])
	   except :
		Dres[i] = 0.
	cmap = ROOT.TCanvas("cmap")
	cmap.SetLogy()
	cmap.SetLogx()
	ghh = ROOT.TH2F("ghh","",80,E[0]-err_E[0],E[-1]-err_E[-1],100,min(obs),max(obs)*2);
	ghh.SetStats(000)
	ghh.SetXTitle("E (MeV) ")
	ghh.SetYTitle("Counts")
	ghh.Draw()	
	tgrobs = ROOT.TGraphErrors(Nbin,E,obs,err_E,obs_err)
	tgrobs.SetLineColor(1)
	tgrobs.SetMarkerColor(1)
	tgrobs.SetMarkerStyle(1)
	tgrobs.Draw("P")

	tgrother = ROOT.TGraph(Nbin,E,other)
	tgrother.SetLineWidth(2)
	tgrother.SetLineStyle(2)
	tgrother.Draw("L")

	tgr = ROOT.TGraph(Nbin,E,src)
	tgr.SetLineColor(1)
	tgr.SetLineWidth(2)
	tgr.Draw("L")

	tgrsum = ROOT.TGraph(Nbin,E,summ)
	tgrsum.SetLineStyle(3)
	tgrsum.SetLineWidth(2)
	tgrsum.Draw("L")


	cres = ROOT.TCanvas("cres")
	cres.SetLogx()
	gh = ROOT.TH2F("gh","",80,E[0]-err_E[0],E[-1]-err_E[-1],100,-1,1);
	gh.SetStats(000)
	gh.SetXTitle("E (MeV) ")
	gh.SetYTitle("Residauls")
	gh.Draw()
	tgres = ROOT.TGraphErrors(Nbin,E,residual,err_E,Dres)
	tgres.Draw("P")
	tg0 = ROOT.TGraph(Nbin,E,zero)
	tg0.SetLineStyle(2)
	tg0.Draw("L")

	# Save the plots in different formats
	cmap.Print(Parameter.PlotName+"_CountsPlot.eps")
	cmap.Print(Parameter.PlotName+"_CountsPlot.C")
	cmap.Print(Parameter.PlotName+"_CountsPlot.png")
	cres.Print(Parameter.PlotName+"_ResPlot.eps")
	cres.Print(Parameter.PlotName+"_ResPlot.C")
	cres.Print(Parameter.PlotName+"_ResPlot.png")
	os.system("rm "+imName)
	image.close()

def RootStyle():
	ROOT.gROOT.SetStyle("Plain");
	ROOT.gStyle.SetPalette(1);


# Canvas
	ROOT.gStyle.SetCanvasColor(10);

# Frame
	ROOT.gStyle.SetFrameBorderMode(0);
	ROOT.gStyle.SetFrameFillColor(0);

# Pad
	ROOT.gStyle.SetPadBorderMode(0);
	ROOT.gStyle.SetPadColor(0);
	ROOT.gStyle.SetPadTopMargin(0.07);
	ROOT.gStyle.SetPadLeftMargin(0.13);
	ROOT.gStyle.SetPadRightMargin(0.11);
	ROOT.gStyle.SetPadBottomMargin(0.1);
	ROOT.gStyle.SetPadTickX(1);#make ticks be on all 4 sides.
	ROOT.gStyle.SetPadTickY(1);

# histogram
	ROOT.gStyle.SetHistFillStyle(0);
	ROOT.gStyle.SetOptTitle(0);

# histogram title
	ROOT.gStyle.SetTitleSize(0.22);
	ROOT.gStyle.SetTitleFontSize(2);
	ROOT.gStyle.SetTitleFont(42);
	ROOT.gStyle.SetTitleFont(62,"xyz");
	ROOT.gStyle.SetTitleYOffset(1.0);
	ROOT.gStyle.SetTitleXOffset(1.0);
	ROOT.gStyle.SetTitleXSize(0.04);
	ROOT.gStyle.SetTitleYSize(0.04);
	ROOT.gStyle.SetTitleX(.15);
	ROOT.gStyle.SetTitleY(.98);
	ROOT.gStyle.SetTitleW(.70);
	ROOT.gStyle.SetTitleH(.05);

# statistics box
	ROOT.gStyle.SetStatFont(42);
	ROOT.gStyle.SetStatX(.91);
	ROOT.gStyle.SetStatY(.90);
	ROOT.gStyle.SetStatW(.15);
	ROOT.gStyle.SetStatH(.15);

# axis labels
	ROOT.gStyle.SetLabelFont(42,"xyz");
	ROOT.gStyle.SetLabelSize(0.035,"xyz");
	ROOT.gStyle.SetGridColor(16);

	ROOT.gStyle.SetLegendBorderSize(0);


########################################################################################################
########################################################################################################
########################################################################################################


def Tgraph(like,Parameter):
	ROOT.gROOT.SetBatch(ROOT.kTRUE) 
	RootStyle()
#	
	Res = Result(like,Parameter)

	E,SED = MakeSED(Res,Parameter)
	Err = MakeError(Res,Parameter)	

        try :
		CountsPlot(Res,Parameter)
	except :
		pass

## Save all in ascii file
#log(E)  log (E**2*dN/dE)   log(E**2*dN/dE_err)  is_dot (0,1) is_upper (0,1) 
	save_file = open(Parameter.PlotName+'.dat','w')
	save_file.write("log(E)  log (E**2*dN/dE)   log(E**2*dN/dE_err)   \n")
	for i in xrange(Parameter.N):
		save_file.write("%12.4e  %12.4e  %12.4e \n" % (E[i], SED[i], Err[i]))

	Fluxp = SED+Err
	Fluxm = SED-Err
	ErrorFlux = array('f',(2*Parameter.N+1)*[0])
	ErrorE = array('f',(2*Parameter.N+1)*[0])

	for i in xrange(Parameter.N):
		ErrorFlux[i] = Fluxp[i]
		ErrorE[i] = E[i]
	for i in xrange(Parameter.N):
		ErrorFlux[Parameter.N+i] = Fluxm[Parameter.N-i-1]
		ErrorE[Parameter.N+i] = E[Parameter.N-i-1]
	ErrorFlux[-1] = Fluxp[0]
	ErrorE[-1] = E[0]


	c_plot = ROOT.TCanvas("c_plot")
	c_plot.SetLogx()
	c_plot.SetLogy()

	gh = ROOT.TH2F("gh","",10000,E[0]*0.8,E[-1]*1.5,100,min(SED[0]-Err[0],SED[-1]-Err[-1])*0.2,max(SED[0]+Err[0],SED[-1]+Err[-1])*3)
	gh.SetStats(000)
	gh.SetTitle("Fermi SED")
	gh.SetXTitle("E [MeV]")
	gh.SetYTitle("E^{2}dN/dE [ erg cm^{-2} s^{-1} ] ")
	gh.Draw()

	tgr = ROOT.TGraph(Parameter.N,E,SED)
	tgr.SetLineWidth(2)
	tgr.SetLineColor(Parameter.LineColor)

	tgerr = ROOT.TGraph(2*Parameter.N+1,ErrorE,ErrorFlux)
	tgerr.SetLineStyle(1)
	tgerr.SetLineWidth(1)
	if (Parameter.AreaColor==0):
		tgerr.SetLineColor(2)
		tgerr.Draw("L")
	else :
		tgerr.SetFillColor(Parameter.AreaColor)
		tgerr.Draw("FL")
		tgr.Draw('L')

	c_plot.Print(Parameter.PlotName+'.C')
	c_plot.Print(Parameter.PlotName+'.eps')
	c_plot.Print(Parameter.PlotName+'.png')

