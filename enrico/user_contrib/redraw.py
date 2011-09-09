#import pylab
from ROOT import TFile,TGraph,TGraphAsymmErrors,SetOwnership,TH1F,gStyle,TLegend,TCanvas,TArrow,TLine
import pyfits
import numpy as num
import os

print 'This is redraw.py version 02'
print 'This macro is designed to use output files from the likeSED family of macros'

#function to redraw just the E^2dN/dE plot, useful if scripts were run in background
#and Plot function called with plot=False
def redrawE2dNdE(filename,flxmod=1.,energymod=1.,covar=False,ytitle='E^{2}dN/dE (erg cm^{-2} s^{-1})',xtitle='Energy (GeV)'):
	gStyle.SetOptTitle(0)
	gStyle.SetOptStat(0)
	gStyle.SetTextFont(132)
	gStyle.SetTitleFont(132,'xyz')
	gStyle.SetLabelFont(132,'xyz')
	flxGraph=getflxGraph(filename,flxmod,energymod)
	flxGraph.SetMarkerStyle(21)
	flxGraph.SetMarkerSize(0.7)
	myfile=pyfits.open(filename)
	if covar:
		modGraph,top,bot=getmodGraph(filename,flxmod,energymod,covar)
		b=myfile[3].data.field("Center Energy")[-1]
		f=myfile[3].data.field('Center Energy')[0]
		fdn=(myfile[3].data.field('model')[0]-myfile[3].data.field('bowtie')[0])*(myfile[3].data.field('center energy')[0]**2*1.602e-3)
		bdn=(myfile[3].data.field('model')[-1]-myfile[3].data.field('bowtie')[-1])*(myfile[3].data.field('center energy')[-1]**2*1.602e-3)
		fup=(myfile[3].data.field('model')[0]+myfile[3].data.field('bowtie')[0])*(myfile[3].data.field('center energy')[0]**2*1.602e-3)
		bup=(myfile[3].data.field('model')[-1]+myfile[3].data.field('bowtie')[-1])*(myfile[3].data.field('center energy')[-1]**2*1.602e-3)
		front=TLine(f*energymod,fdn*flxmod,f*energymod,fup*flxmod)
		back=TLine(b*energymod,bdn*flxmod,b*energymod,bup*flxmod)
		back.SetLineColor(2)
		front.SetLineColor(2)
	else:
		modGraph=getmodGraph(filename,flxmod,energymod,covar)
	#get min and max for dummy hist x-axis
	xmin=myfile[2].data.field('bin low edge')[0]*0.8*energymod
	xmax=myfile[2].data.field('bin high edge')[-1]*1.2*energymod
	#get min and max for dummy hist y-axis, doesn't check model info but can adjust as necessary
	ymax=max([(V+E)*flxmod for V,E in zip(myfile[1].data.field('E^2dN/dE'),myfile[1].data.field('E^2dN/dE_Error'))])
	ymin=min([(V-E)*flxmod for V,E in zip(myfile[1].data.field('E^2dN/dE'),myfile[1].data.field('E^2dN/dE_Error'))])
	ymin=(ymax*1e-3 if ymin<=0 else ymin)
	dummy=TH1F('dummy','',1000,xmin,xmax)
	dummy.GetYaxis().SetRangeUser(ymin*.5,ymax*2)
	dummy.GetYaxis().SetTitle(ytitle)
	dummy.GetYaxis().CenterTitle()
	dummy.GetXaxis().SetTitle(xtitle)
	dummy.GetXaxis().CenterTitle()
	dummy.GetXaxis().SetTitleSize(0.05)
	dummy.GetXaxis().SetTitleOffset(0.9)
	dummy.GetYaxis().SetTitleSize(0.05)
	flxarrows=[]
	for f,e,c in zip(myfile[1].data.field('E^2dN/dE'),myfile[1].data.field('E^2dN/dE_Error'),myfile[1].data.field('Center Energy')):
		if e==0:
			flxarrows+=[TArrow(c*energymod,f*flxmod,c*energymod,f*flxmod*0.5,0.025,'|->')]
	flxCanvas=TCanvas("flxCanvas","E^{2}dN/dE Spectrum",700,500)
	flxCanvas.SetFillColor(0)
	flxCanvas.SetFrameBorderMode(0)
	flxCanvas.SetBorderMode(0)
	flxCanvas.SetTicks(1,1)
	flxCanvas.SetLogx()
	flxCanvas.SetLogy()
	dummy.Draw()
	flxGraph.Draw('psame')
	modGraph.Draw('csame')
	if covar:
		top.Draw('csame')
		bot.Draw('csame')
		front.Draw()
		back.Draw()
	for x in flxarrows:
		x.Draw()
	leg=TLegend(0.15,0.15,0.55,0.3)
	leg.SetBorderSize(1)
	leg.SetFillStyle(0)
	leg.AddEntry(flxGraph,'Energy Band Fits','p')
	leg.AddEntry(modGraph,'Maximum Likelihood Model','l')
	leg.Draw()
	SetOwnership(dummy,False)
	SetOwnership(flxGraph,False)
	SetOwnership(modGraph,False)
	SetOwnership(flxCanvas,False)
	SetOwnership(leg,False)
	if covar:
		SetOwnership(top,False)
		SetOwnership(bot,False)
		SetOwnership(front,False)
		SetOwnership(back,False)
	for x in flxarrows:
		SetOwnership(x,False)
	myfile.close()
	return

#function to simply redraw the three canvases from the (bd)likeSED macro
def redrawSED(filename):
	if fileCheck(filename)==0:
		return
	tfile=TFile(filename)
	flx=tfile.Get('flxCanvas')
	cnt=tfile.Get('cntCanvas')
	Ts=tfile.Get('tsCanvas')
	Ts.Draw()
	cnt.Draw()
	flx.Draw()
	SetOwnership(Ts,False)
	SetOwnership(cnt,False)
	SetOwnership(flx,False)
	return

#function to return a TGraphAsymmErrors corresponding to the E^2dN/dE points from the (bd)likeSED macro with potential to change units
def getflxGraph(filename,flxmod=1.,energymod=1.): #flxmod and energymod are the necessary multipliers to go from default units to units of your choice
	if fileCheck(filename)==0:
		return
	#access the data file and columns necessary
	sfile=pyfits.open(filename)
	pts=sfile[1].data.field('E^2dN/dE')
	errs=sfile[1].data.field('E^2dN/dE_Error')
	cents=sfile[1].data.field('Center Energy')
	low=sfile[2].data.field('Bin Low Edge')
	high=sfile[2].data.field('Bin High Edge')
	#need to create some empty lists, need to refill even if flxmod=1, for instance, since they need to be cast as floats (before the num.array step) for some reason
	ene=[]
	plus=[]
	minus=[]
	flx=[]
	dflx=[]
	#fill the energy specific lists and create arrays for the TGraph
	for c,l,h in zip(cents,low,high):
		ene+=[float(energymod*c)]
		plus+=[float(energymod*(h-c))]
		minus+=[float(energymod*(c-l))]
	Ene=num.array(ene)
	Eplus=num.array(plus)
	Eminus=num.array(minus)
	#Fill the flux specific lists and create arrays for the TGraph
	for f,d in zip(pts,errs):
		flx+=[float(flxmod*f)]
		dflx+=[float(flxmod*d)]
	Flx=num.array(flx)
	dFlx=num.array(dflx)
	#Make the TGraph
	flxGraph=TGraphAsymmErrors(len(Ene),Ene,Flx,Eminus,Eplus,dFlx,dFlx)
	return flxGraph

#function to return a TGraph for the best fit model and, if computed, bowtie plots
def getmodGraph(filename,flxmod=1.,energymod=1.,covar=False): #similar arguments as for getflxGraph function, to get bowtie graphs need covar=True
	if fileCheck(filename)==0:
		return
	#access the data file and columns necessary
	sfile=pyfits.open(filename)
	modf=sfile[3].data.field('Model')
	cents=sfile[3].data.field('Center Energy')
	#make empty lists to fill for later use in TGraph
	model=[]
	ene=[]
	#Fill lists for model and energies
	for m,c in zip(modf,cents):
		model+=[float(flxmod*(m*c**2*1.602e-3))] #Note that the model is in differential flux, need to multiply by E^2 and convert to erg
		ene+=[float(energymod*c)]
	#make arrays for TGraph
	Model=num.array(model)
	Ene=num.array(ene)
	#make model TGraph
	modGraph=TGraph(len(Ene),Ene,Model)
	#Create more list and TGraphs if covar=True
	if covar==True:
		try:
			bt=sfile[3].data.field('BowTie')
			mminus=[]
			mplus=[]
			for b,c,m in zip(bt,cents,modf):
				mplus+=[float(flxmod*((m+b)*c**2*1.602e-3))]
				mminus+=[float(flxmod*((m-b)*c**2*1.602e-3))]
			Mplus=num.array(mplus)
			Mminus=num.array(mminus)
			plusGraph=TGraph(len(Ene),Ene,Mplus)
			minusGraph=TGraph(len(Ene),Ene,Mminus)
			plusGraph.SetLineColor(2)
			minusGraph.SetLineColor(2)
			return modGraph,plusGraph,minusGraph
		except: #just in case someone wants a bowtie plot when that information was not calcluated previously
			print 'Model error information was not computed when likeSED macro was run.  Either use covar=False or rerun likeSED macro with CoVar=True for fullFit(...) function.'
			return
		
	return modGraph

#check to make sure files exist, do first to avoid errors later
def fileCheck(file):
	if (not os.access(file,os.F_OK)):
		print "Error:  %s not found." %file
		return 0
	return 1

#######################USAGE NOTES##############################
#
#This is meant to be a supplement to the likeSED and bdlikeSED python macros I've written to aid in replotting
#the output spectra and in using the spectra on larger band SEDs.
#
#There are 3 functions of interest:
#redrawSED merely accesses the root file, which is automatically made when calling the Plot function of a bdlikeSED or likeSED object,
#to redraw the default canvases.
#
#getflxGraph accesses the fits file, also made automatically when calling the Plot function, and makes a TGraphAsymmErrors which can be used to plot on a full band SED.
#
#getmodGraph accesses the same fits file getflxGraph does but makes a TGraph for the full energy range model (and bowtie plot) for use on a full band SED.
#
#To use any of these functions, the first step is: from redraw import *
#
#The redraw function:
#There is one argument for this function, a string which corresponds to the filename of the output root file (srcname_nbins_histos.root, where one replaces
#srcname with the name of your source as it appears in the xml file sans any spaces and n with the number of bins over the full energy range you gave to 
#the (bd)likeInput object).
#
#The getflxGraph function:
#This has three arguments, in order they are 1)a string which corresponds to the name of the default output fits file from (bd)likeSED
#for likeSED (unbinned) this is srcname_nbins_likeSEDout.fits and for bdlikeSED (binned) this is srcname_nbins_bdlikeSEDout.fits, with
#the subsitutions for 'srcname' and 'n' as described above for the redraw function; 2)a factor to multiply the flux by in the case that one
#desires to use different units, the default units are erg cm^-2 s^-1; & 3)a factor to multiply the energy by in the case that one desires to
#use different units, the default units are GeV.
#To use the function without changing units one simply does: mygraph=getflxGraph('myfile.fits')
#To change only the flux units: mygraph=getflxGraph('myfile.fits',1.602e3)
#To change both units: mygraph=getflxGraph('myfile.fits',1.602e3,1000)
#To change only the energy units: mygraph=getflxGraph('myfile.fits',1.,1000) or mygraph=getflxGraph('myfile.fits',energymod=1000)
#Note that you always have to set this equal to something as it returns a TGraphAsymmErrors object, use option 'p' to get points (i.e. to draw on
#a plot with other points use mygraph.Draw('psame')).
#Note that if there are upper limits to some of the points this does not create the TArrows for those, they will appear as points with vertical error
#bars of length zero.
#
#The getmodGraph function:
#This has four agruments, the first three are exactly the same as described above for the getflxGraph function but returns simply a TGraph.  The fourth argument denotes whether
#or not you want a bowtie plot as well.  The default for this is False, set to True if you want a bowtie plot.  Note that to get a bowtie plot you had 
#to have run the fullFit function of the (bd)likeInput object with the option CoVar=True initially or else that information was not calculated.
#To use the function without the bowtie plot follow the same prescription as outlined above for the getflxGraph function.
#To use the function to get a bowtie plot with no flux or energy unit changing: mymodel,bowtieplus,bowtieminus=getmodGraph('myfile.fits',covar=True)
#To use the function to get a bowtie plot and change units: mymodel,bowtieplus,bowtieminus=getmodGraph('myfile.fits',1.602e3,1000,True)
#The argument names are the same as outlined above if one wants to only change one set of units.
#Note that, when using covar=True, this returns three(3) TGraphs.  To draw them as curves us option 'c' (i.e. to draw on a plot with other data
#use mymodel.Draw('csame');bowtieplus.Draw('csame');bowtieminus.Draw('csame')).
#By default the two bowtie TGraphs are set to have red lines.
#################################################################