class srcList:
	#arguments are:
	#sources (string, filename of LAT source list fits file in catalog format)
	#ft1 (string, filename of event file for which the xml will be used, only used to extract ROI info)
	#out (string, name of output xml file, defaults to mymodel.xml)
	def __init__(self,sources,ft1,out='mymodel.xml'):
		if not fileCheck(sources): #check that file exists
			print "Error:  %s not found." %sources
			return
		if fileCheck(out):
			print 'Warning: %s already exists, file will be overwritten if you proceed with makeModel.' %out
		self.srcs=sources
		self.out=out
		self.roi=getPos(ft1)
	
	#define a quick print function to make sure everything looks irght
	def Print(self):
		print 'Source list file: ',self.srcs
		print 'Output file name: ',self.out
		print 'Selecting %s degrees around (ra,dec)=(%s,%s)' %(self.roi[2],self.ra,self.dec)
	
	#make the xml file
	#GDfile is galactic diffuse model and ISOfile is optional Isotropic template file
	#radLim is an angular radius limit beyond which all the spectral parameters of all sources will be fixed
	#tslim, for 2FGL only, sources with TS below this limit will have spectral parameters fixed
	#signif, for 1FGL only, sources with signiv_avg below this value will have spectral parameters fixed
	#psForce, if true forces extended sources to be case as point sources
	def makeModel(self,GDfile='$(GLAST_EXT)/extFiles/v0r9/galdiffuse/gll_iem_v02.fit',GDname='GAL_v02',ISOfile=None,ISOname='Extragalactic Diffuse',extDir='',radLim=-1,signif=4,psForce=False): #the version number may differ, but default syntax from ModelEditor didn't work
		#quick check to see if 1FGL or 2FGL fits file
		self.radLim=(self.roi[2] if radLim<=0 else radLim)
		self.psF=psForce
		self.extD=extDir
		mycheck=pyfits.open(self.srcs)
		try:
			test=mycheck[1].data.field('Cutoff')
			mycheck.close()
			print 'Creating file and adding sources for 2FGL'
			use=2
		except:
			mycheck.close()
			print 'Creating file and adding sources for 1FGL'
			use=1
		if use==1:
			addSrcs1(self,GDfile,GDname,ISOfile,ISOname,signif)
		else:
			addSrcs2(self,GDfile,GDname,ISOfile,ISOname,signif)
	
import pyfits
import os
from xml.dom.minidom import parseString as pS
from numpy import floor,log10,cos,sin,arccos,pi,array,log
acos=arccos
#import ROOT #note that this is only done to turn tab completion on for functions and filenames
print "This is make2FGLxml version 04."
#print "NOTE: You must have run gtselect on the event file you use as input."
d2r=pi/180.

#function to cycle through the source list and add point source entries
def addSrcs2(sL,GD,GDn,ISO,ISOn,signif):
	model=open(sL.out,'w') #open file in write mode, overwrites other files of same name
	file=pyfits.open(sL.srcs) #open source list file and access necessary fields, requires LAT source catalog definitions and names
	data=file[1].data
	extendedinfo=file[4].data
	extName=extendedinfo.field('Source_Name')	
	extFile=extendedinfo.field('Spatial_Filename')
	name=data.field('Source_Name')
	ExtName=data.field('Extended_Source_Name')
	ra=data.field('RAJ2000')
	dec=data.field('DEJ2000')
	flux=data.field('Flux_Density')
	pivot=data.field('Pivot_Energy')
	index=data.field('Spectral_Index')
	cutoff=data.field('Cutoff')
	spectype=data.field('SpectrumType')
	beta=data.field('beta')
	#extended=data.field('Extended')
	sigma=data.field('Signif_Avg')
	model.write('<?xml version="1.0" ?>\n')
	model.write('<source_library title="source library">\n')
	model.write('\n<!-- Point Sources -->\n')
	step=(sL.roi[2]+5)/5 #divide ROI radius plus 5 degrees into 5 steps for ordering of sources, get next highest int for ease
	i=1
	radii=[]
	ptSrcNum=0
	extSrcNum=0
	while i<6:
		if i*step<=sL.roi[2]+5:
			radii+=[step*i]
		else:
			radii+=[sL.roi[2]+5] #just in case of rounding errors
		i+=1
	for x in radii:
		if x==sL.roi[2]+5.:
			model.write('\n<!-- Sources between [%s,%s] degrees of ROI center -->\n' %(x-step,x))
		else:
			model.write('\n<!-- Sources between [%s,%s) degrees of ROI center -->\n' %(x-step,x))
		for n,f,i,r,d,p,c,t,b,S,EN in zip(name,flux,index,ra,dec,pivot,cutoff,spectype,beta,sigma,ExtName):
			dist=angsep(sL.roi[0],sL.roi[1],r,d) #check that source is within ROI radius + 5 degress of ROI center
			if r==sL.roi[0] and d==sL.roi[1]:
				dist=0.0
			if (dist<x and dist>=x-step) or (x==sL.roi[2]+5. and dist==x):
				sn=''
				for N in n.split(' '):
					sn+=N
				if EN!='' and not sL.psF:
					extSrcNum+=1
					if EN[0].isdigit():#use the extended source name to set these apart even further
						Name='<source name="_%s" type="DiffuseSource">\n' %EN
					else:
						Name='<source name="%s" type="DiffuseSource">\n' %EN
				else:
					ptSrcNum+=1
					if sn[0].isdigit():
						Name='<source name="_%s" type="PointSource">\n' %sn
					else:
						Name='<source name="%s" type="PointSource">\n' %sn
				if EN!='VelaX' and EN!='MSH15-52':
					if t=='PowerLaw':
						spec=PLspec(sL,f,i,p,dist,S,signif)
					elif t=='PowerLaw2':#no value for flux from 100 MeV to 100 GeV in fits file0
						if i!=1.:#so calculate it by integrating PowerLaw spectral model
							F=f*p**i/(-i+1)*(1.e5**(-i+1)-1.e2**(-1+1))
						else:
							F=f*p*log(1.e3)
						spec=PL2spec(sL,F,i,dist,S,signif)
					elif t=='LogParabola':
						spec=LPspec(sL,f,i,p,b,dist,S,signif)
					else:
						spec=COspec(sL,f,i,p,c,dist,S,signif)
				else:
					if EN=='VelaX':
						spec=VXspec(sL,i,dist)
					else:
						spec=MSHspec(sL,i,dist)
				if EN!='' and not sL.psF:
					eFile=None
					for exN,exf in zip(extName,extFile):
						if exN.replace(' ','')==EN:
							if len(sL.extD)>0:
								if sL.extD[-1]=='/':
									eFile=sL.extD+exf
								else:
									eFile=sL.extD+'/'+exf
							else:
								eFile=exf
							break
					if eFile==None:
						skydir='\t<spatialModel file="%s" type="SpatialMap">\n'%sL.extD
						print 'Extended source %s (%s) in ROI, could not find name of matching template file.'%(EN,sn)
					else:
						skydir='\t<spatialModel file="%s" type="SpatialMap">\n'%eFile
						print 'Extended source %s (%s) in ROI, verify that %s is the correct path.'%(EN,sn,eFile)
					skydir+='\t\t<parameter free="0" max="1000" min="0.001" name="Prefactor" scale="1" value="1"/>\n'
					skydir+='\t</spatialModel>\n'
				else:
					skydir='\t<spatialModel type="SkyDirFunction">\n'
					skydir+='\t\t<parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="%s"/>\n' %r
					skydir+='\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="%s"/>\n' %d
					skydir+='\t</spatialModel>\n'
				skydir+='</source>'
				(src,)=(Name+spec+skydir,)
				ptsrc=pS(src).getElementsByTagName('source')[0]
				ptsrc.writexml(model)
				model.write('\n')
	file.close() #close file
	if not sL.psF:
		print 'Added %i point sources and %i extended sources'%(ptSrcNum,extSrcNum)
		if extSrcNum>0:
			print 'If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument psForce=True'
	else:
		print 'Added %i point sources, note that any extended sources in ROI were modeled as point sources becaue psForce option was set to True'%ptSrcNum
	#add galactic diffuse with PL spectrum, fix index to zero for general use, those who want it to be free can unfreeze parameter manually
	model.write('\n<!-- Diffuse Sources -->\n')
	Name='\n<source name="%s" type="DiffuseSource">\n' %GDn
	spec='\t<spectrum type="PowerLaw">\n'
	spec+='\t\t<parameter free="1" max="10" min="0" name="Prefactor" scale="1" value="1"/>\n'
	spec+='\t\t<parameter free="0" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n'
	spec+='\t\t<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n'
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel file="%s" type="MapCubeFunction">\n' %GD
	skydir+='\t\t<parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	galdiff=pS(src).getElementsByTagName('source')[0]
	galdiff.writexml(model)
	model.write('\n')
	if ISO==None: #null means no file so assume user wants an isotropic power law component
		Name='<source name="%s" type="DiffuseSource">\n' %ISOn
		spec='\t<spectrum type="PowerLaw">\n'
		spec+='\t\t<parameter free="1" max="1e3" min="1e-3" name="Prefactor" scale="1e-7" value="1"/>\n'
		spec+='\t\t<parameter free="1" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>\n'
		spec+='\t\t<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n'
		spec+='\t</spectrum>\n'
	else: #if a file is given assume it is an isotropic template
		Name='<source name="%s" type="DiffuseSource">\n' %ISOn
		spec='\t<spectrum type="FileFunction" file="%s">\n' %ISO
		spec+='\t\t<parameter free="1" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n'
		spec+='\t</spectrum>\n'
	#both of the above options use the same spatial model
	skydir='\t<spatialModel type="ConstantValue">\n'
	skydir+='\t\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	iso=pS(src).getElementsByTagName('source')[0]
	iso.writexml(model)	
	model.write('\n</source_library>')
	model.close()
	return

def PLspec(sL,f,i,p,dist,sig,siglim):
	fscale=int(floor(log10(f)))
	spec='\t<spectrum type="PowerLaw">\n'
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	elif(dist>sL.radLim):
		spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	elif sig<siglim:
		spec+='\t<!-- Source Signif=%s is less than specified limit of %s -->\n' %(sig,siglim)
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	else:
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="1" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="Scale" scale="1.0" value="%f"/>\n' %p
	spec+='\t</spectrum>\n'
	return spec

def PL2spec(sL,F,i,dist,sig,siglim):
	fscale=int(floor(log10(F)))
	spec='\t<spectrum type="PowerLaw2">\n'
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
		spec+='\t\t<parameter free="0" max="5" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	elif(dist>sL.radLim):
		spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
		spec+='\t\t<parameter free="0" max="5" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	elif sig<siglim:
		spec+='\t<!-- Source Signif=%s is less than specified limit of %s -->\n' %(sig,siglim)
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
		spec+='\t\t<parameter free="0" max="5" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	else:
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Integral" scale="1e%i" value="%s"/>\n'%(fscale,F/10**fscale)
		spec+='\t\t<parameter free="1" max="5" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="LowerLimit" scale="1" value="1e2"/>\n'
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="UpperLimit" scale="1" value="1e5"/>\n'
	spec+='\t</spectrum>\n'
	return spec

def VXspec(sL,i,dist):
	spec='\t<spectrum type="PowerLaw2">\n'
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
	spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e-7" value="4.73"/>\n'
	spec+='\t\t<parameter free="0" max="5" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="LowerLimit" scale="1" value="1e2"/>\n'
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="UpperLimit" scale="1" value="2e5"/>\n'
	spec+='\t</spectrum>\n'
	return spec

def MSHspec(sL,i,dist):
	spec='\t<spectrum type="PowerLaw2">\n'
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
	spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Integral" scale="1e-8" value="0.291"/>\n'
	spec+='\t\t<parameter free="0" max="5" min="0" name="Index" scale="-1.0" value="%s"/>\n' %i
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="LowerLimit" scale="1" value="1e3"/>\n'
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="UpperLimit" scale="1" value="1e5"/>\n'
	spec+='\t</spectrum>\n'
	return spec


def COspec(sL,f,i,p,c,dist,sig,siglim):
	fscale=int(floor(log10(f)))
	spec='\t<spectrum type="PLSuperExpCutoff">\n'
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="0" max="1e5" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	elif(dist>sL.radLim):
		spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="0" max="1e5" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	elif sig<siglim:
		spec+='\t<!-- Source Signif=%s is less than specified limit of %s -->\n' %(sig,siglim)
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="0" max="1e5" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	else:
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="1" max="5.0" min="0.0" name="Index1" scale="-1.0" value="%s"/>\n' %i
		if c<=1e5:
			spec+='\t\t<parameter free="1" max="1e5" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%c
		else:
			spec+='\t\t<parameter free="0" max="%.2e" min="1e2" name="Cutoff" scale="1.0" value="%f"/>\n'%(2.*c,c)
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="Scale" scale="1.0" value="%f"/>\n'%p
	spec+='\t\t<parameter free="0" max="5" min="0" name="Index2" scale="1.0" value="1"/>\n'
	spec+='\t</spectrum>\n'
	return spec

def LPspec(sL,f,i,p,b,dist,sig,siglim):
	fscale=int(floor(log10(f)))
	spec='\t<spectrum type="LogParabola">\n'
	spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
	if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
		spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	elif(dist>sL.radLim):
		spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	elif sig<siglim:
		spec+='\t<!-- Source Signif=%s is less than specified limit of %s -->\n' %(sig,siglim)
		spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="0" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	else:
		spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="norm" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
		spec+='\t\t<parameter free="1" max="5.0" min="0.0" name="alpha" scale="1.0" value="%s"/>\n' %i
		spec+='\t\t<parameter free="1" max="10.0" min="0.0" name="beta" scale="1.0" value="%s"/>\n'%b
	spec+='\t\t<parameter free="0" max="5e5" min="30" name="Eb" scale="1.0" value="%s"/>\n'%p
	spec+='\t</spectrum>\n'
	return spec

def addSrcs1(sL,GD,GDn,ISO,ISOn,signif):
	model=open(sL.out,'w') #open file in write mode, overwrites other files of same name
	file=pyfits.open(sL.srcs) #open source list file and access necessary fields, requires LAT source catalog definitions and names
	data=file[1].data
	try:
		name=data.field('Source_Name')
	except:
		name=data.field('NickName')
	ra=data.field('RA')
	dec=data.field('DEC')
	flux=data.field('Flux_Density')
	pivot=data.field('Pivot_Energy')
	index=data.field('Spectral_Index')
	sigma=data.field('Signif_Avg')
	model.write('<?xml version="1.0" ?>\n')
	model.write('<source_library title="source library">\n')
	model.write('\n<!-- Point Sources -->\n')
	step=(sL.roi[2]+5)/5 #divide ROI radius plus 5 degrees into 5 steps for ordering of sources, get next highest int for ease
	i=1
	radii=[]
	while i<6:
		if i*step<=sL.roi[2]+5:
			radii+=[step*i]
		else:
			radii+=[sL.roi[2]+5] #just in case of rounding errors
		i+=1
	for x in radii:
		if x==sL.roi[2]+5.:
			model.write('\n<!-- Sources between [%s,%s] degrees of ROI center -->\n' %(x-step,x))
		else:
			model.write('\n<!-- Sources between [%s,%s) degrees of ROI center -->\n' %(x-step,x))
		for n,f,i,r,d,p,s in zip(name,flux,index,ra,dec,pivot,sigma):
			dist=angsep(sL.roi[0],sL.roi[1],r,d) #check that source is within ROI radius + 5 degress of ROI center
			if (dist<x and dist>=x-step) or (x==sL.roi[2]+5. and dist==x):
				fscale=int(floor(log10(f)))
				sn=''
				for N in n.split(' '):
					sn+=N
				if sn[0].isdigit():
						Name='<source name="_%s" type="PointSource">\n' %sn
				else:
					Name='<source name="%s" type="PointSource">\n' %sn
				spec='\t<spectrum type="PowerLaw">\n'
				spec+='\t<!-- Source is %s degrees away from ROI center -->\n' %dist
				if(dist>sL.roi[2]): #if beyond ROI, shouldn't attempt to fit parameters
					spec+='\t<!-- Source is outside ROI, all parameters should remain fixed -->\n'
					spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
					spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
				elif(dist>sL.radLim):
					spec+='\t<!-- Source is outside specified radius limit of %s -->\n'%sL.radLim
					spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
					spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
				elif s<signif:
					spec+='\t<!-- Source significance %s is less than specified limit of %s -->\n' %(s,signif)
					spec+='\t\t<parameter free="0" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
					spec+='\t\t<parameter free="0" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
				else:
					spec+='\t\t<parameter free="1" max="1e4" min="1e-4" name="Prefactor" scale="1e%i" value="%s"/>\n' %(fscale,f/10**fscale)
					spec+='\t\t<parameter free="1" max="5.0" min="0.0" name="Index" scale="-1.0" value="%s"/>\n' %i
				spec+='\t\t<parameter free="0" max="5e5" min="30" name="Scale" scale="1.0" value="%f"/>\n' %p
				spec+='\t</spectrum>\n'
				skydir='\t<spatialModel type="SkyDirFunction">\n'
				skydir+='\t\t<parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="%s"/>\n' %r
				skydir+='\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="%s"/>\n' %d
				skydir+='\t</spatialModel>\n'
				skydir+='</source>'
				(src,)=(Name+spec+skydir,)
				ptsrc=pS(src).getElementsByTagName('source')[0]
				ptsrc.writexml(model)
				model.write('\n')
	file.close() #close file
	#add galactic diffuse with PL spectrum, fix index to zero for general use, those who want it to be free can unfreeze parameter manually
	model.write('\n<!-- Diffuse Sources -->\n')
	Name='\n<source name="%s" type="DiffuseSource">\n' %GDn
	spec='\t<spectrum type="PowerLaw">\n'
	spec+='\t\t<parameter free="1" max="10" min="0" name="Prefactor" scale="1" value="1"/>\n'
	spec+='\t\t<parameter free="0" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n'
	spec+='\t\t<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n'
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel file="%s" type="MapCubeFunction">\n' %GD
	skydir+='\t\t<parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	galdiff=pS(src).getElementsByTagName('source')[0]
	galdiff.writexml(model)
	model.write('\n')
	if ISO==None: #null means no file so assume user wants an isotropic power law component
		Name='<source name="%s" type="DiffuseSource">\n' %ISOn
		spec='\t<spectrum type="PowerLaw">\n'
		spec+='\t\t<parameter free="1" max="1e3" min="1e-3" name="Prefactor" scale="1e-7" value="1"/>\n'
		spec+='\t\t<parameter free="1" max="-1.0" min="-3.5" name="Index" scale="1.0" value="-2.1"/>\n'
		spec+='\t\t<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n'
		spec+='\t</spectrum>\n'
	else: #if a file is given assume it is an isotropic template
		Name='<source name="%s" type="DiffuseSource">\n' %ISOn
		spec='\t<spectrum type="FileFunction" file="%s">\n' %ISO
		spec+='\t\t<parameter free="1" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n'
		spec+='\t</spectrum>\n'
	#both of the above options use the same spatial model
	skydir='\t<spatialModel type="ConstantValue">\n'
	skydir+='\t\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	iso=pS(src).getElementsByTagName('source')[0]
	iso.writexml(model)	
	model.write('\n</source_library>')
	model.close()
	return

#this function searches the header of the ft1 to find the Position keyword and extract the ra and dec values
def getPos(ft1):
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
	try:
		ra,dec,rad=header[keyword].strip('CIRCLE()').split(',') #gets rid of the circle and parenthesis part and splits around the comma
		float(ra)
	except:
		ra,dec,rad=header[keyword].strip('circle()').split(',')
	file.close()
	return float(ra),float(dec),float(rad)
	
#calculates the angular separation between two points on the sky
def angsep(ra1,dec1,ra2,dec2):
	ra1*=d2r
	dec1*=d2r
	ra2*=d2r
	dec2*=d2r
	diffCosine=cos(dec1)*cos(dec2)*cos(ra1-ra2)+sin(dec1)*sin(dec2)
	dC='%.10f'%diffCosine#when the source is right at the center of the roi python sometimes adds extraneous digits at the end of the value i.e. instead of 1.0
	#it returns 1.0000000000000024, which throws an error with the acos function
	return acos(float(dC))/d2r #returns values between 0 and pi radians

#Check if a given file exists or not
def fileCheck(file):
	if (not os.access(file,os.F_OK)):
		return 0
	return 1