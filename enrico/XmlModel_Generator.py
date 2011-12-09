# Author David Sanchez dsanchez@llr.in2p3.fr
# script to look for near source of given RA and DEC in a catalog
# begun Oct 2009

import sys, getopt, subprocess, xml.dom.minidom, os, time, math
from xml.sax import ContentHandler, make_parser
from math import *
import sys
import pyfits,xml.dom.minidom,numpy

import enrico.environ as env
from enrico.config import get_config


def fluxScale(flux_value):
    return 10**math.floor(math.log10(flux_value)+0.5)

def meanEnergy(emin, emax, index_value):
    x=emax/emin;
    if index_value==-2.0:
        eflux = emax*math.log(x)/(x-1)
    elif index_value==-1.0:
        eflux = emin*(x-1)/math.log(x)
    else:
        eflux = emin*(index_value+1)/(index_value+2)*\
                (x**(index_value+2)-1)/(x**(index_value+1)-1)
    return eflux

def addParameter(el, name, free, value, scale, min, max):
    doc = el.ownerDocument
    param = doc.createElement('parameter')
    param.setAttribute('name',name)
    param.setAttribute('free','%d'%free)
    param.setAttribute('scale','%g'%scale)
    param.setAttribute('value','%g'%value)
    param.setAttribute('max','%g'%max)
    param.setAttribute('min','%g'%min)
    el.appendChild(param)

def addDiffusePL(lib, file, free=1, value=1.0, scale=1.0, max=10.0, min=1.0,
               name = 'EG_v02'):
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('file',file)
    spec.setAttribute('type','FileFunction')
    addParameter(spec, 'Normalization', free, 1, 1, 0.001, 1000)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type','ConstantValue')
    addParameter(spatial,'Value',0,1,1.0,0.0,10.0)
    src.appendChild(spatial)
    lib.appendChild(src)

def addGalprop(lib, file, free=1, value=1.0, scale=1.0, max=10.0, min=.010,
               name = 'GAL_v02'):
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','ConstantValue')
    addParameter(spec, 'Value', free, value, scale, min, max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('file',file)
    spatial.setAttribute('type','MapCubeFunction')
    addParameter(spatial, 'Normalization', 0, 1, 1, 0.001, 1000)
    src.appendChild(spatial)
    lib.appendChild(src)


def addDiffusePL2(lib, file, name ,emin=200, emax=3e5,
                   flux_free=1, flux_value=1.6e-8, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index_free=1, index_value=-2.0,
                   index_min=-5.0, index_max=-0.5):

    elim_min = 30;
    elim_max = 300000;
    if emin<elim_min:
        elim_min = emin
    if emax>elim_max:
        elim_max = emax
    if flux_scale == 0:
        flux_scale=fluxScale(flux_value)
    flux_value /= flux_scale
 
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','DiffuseSource')   
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','PowerLaw2')
    addParameter(spec,'Integral',
                 flux_free,flux_value,flux_scale,flux_min,flux_max)
    addParameter(spec,'Index',index_free,index_value,1.0,index_min,index_max)
    addParameter(spec,'LowerLimit',0,emin,1.0,elim_min,elim_max)
    addParameter(spec,'UpperLimit',0,emax,1.0,elim_min,elim_max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('file',file)
    spatial.setAttribute('type','SpatialMap')
    addParameter(spatial, 'Prefactor', 0, 1, 1, 0.001, 1000)
    src.appendChild(spatial)
    lib.appendChild(src)

#def addPSPowerLaw1(lib, name, ra, dec, emin=200, emax=100000, eflux=0,
#                   flux_free=1, flux_value=1e-9, flux_scale=0,
#                   flux_max=1000.0, flux_min=1e-5,
#                   index_free=1, index_value=-2.0,
#                   index_min=-5.0, index_max=-0.5):
def addPSPowerLaw1(lib, name, ra, dec, eflux=0,
                   flux_free=1, flux_value=1e-9, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index_free=1, index_value=-2.0,
                   index_min=-5.0, index_max=-0.5):
    elim_min = 30;
    elim_max = 300000;
#    if emin<elim_min:
#        elim_min = emin
#    if emax>elim_max:
#        elim_max = emax 
#    if eflux==0:
#        eflux =2e5# meanEnergy(emin,emax,index_value)
#        flux_value *= (eflux/100.0)**index_value
#    else :
#	flux_value *= 0.1*(eflux/100.0)**index_value
    if flux_scale == 0:
        flux_scale=fluxScale(flux_value)
    flux_value /= flux_scale        
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','PointSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','PowerLaw')
    addParameter(spec,'Prefactor',
                 flux_free,flux_value,flux_scale,flux_min,flux_max)
    addParameter(spec,'Index',index_free,index_value,1.0,index_min,index_max)
    addParameter(spec,'Scale',0,eflux,1.0,elim_min,elim_max)
#    addParameter(spec,'Scale',0,2e5,1.0,elim_min,elim_max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type','SkyDirFunction')
    addParameter(spatial,'RA',0,ra,1.0,-360.0,360.0)
    addParameter(spatial,'DEC',0,dec,1.0,-90.0,90.0)
    src.appendChild(spatial)
    lib.appendChild(src)

def addPSPowerLaw2(lib, name, ra, dec, emin=200, emax=3e5,
                   flux_free=1, flux_value=1.6e-8, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index_free=1, index_value=-2.0,
                   index_min=-5.0, index_max=-0.5):
    elim_min = 30;
    elim_max = 300000;
    if emin<elim_min:
        elim_min = emin
    if emax>elim_max:
        elim_max = emax
    if flux_scale == 0:
        flux_scale=fluxScale(flux_value)
    flux_value /= flux_scale
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','PointSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','PowerLaw2')
    addParameter(spec,'Integral',
                 flux_free,flux_value,flux_scale,flux_min,flux_max)
    addParameter(spec,'Index',index_free,index_value,1.0,index_min,index_max)
    addParameter(spec,'LowerLimit',0,emin,1.0,elim_min,elim_max)
    addParameter(spec,'UpperLimit',0,emax,1.0,elim_min,elim_max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type','SkyDirFunction')
    addParameter(spatial,'RA',0,ra,1.0,-360.0,360.0)
    addParameter(spatial,'DEC',0,dec,1.0,-90.0,90.0)
    src.appendChild(spatial)
    lib.appendChild(src)

def addPSLogparabola(lib, name, ra, dec,  enorm=300,
                   norm_free=1, norm_value=1e-9, norm_scale=0,
                   norm_max=1000.0, norm_min=1e-5,
                   alpha_free=1, alpha_value=1.0,
                   alpha_min=.5, alpha_max=5.,
                   beta_free=1, beta_value=1.0,
                   beta_min=0.5, beta_max=5.0):
    elim_min = 30;
    elim_max = 300000;

    if enorm==0:
        enorm =2e5# meanEnergy(emin,emax,index_value)
        norm_value *= (enorm/100.0)**alpha_value
    if norm_scale == 0:
        norm_scale=fluxScale(norm_value)
    norm_value /= norm_scale        
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','PointSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','LogParabola')
    addParameter(spec,'norm',
                 norm_free,norm_value,norm_scale,norm_min,norm_max)
    addParameter(spec,'alpha',alpha_free,alpha_value,1.0,alpha_min,alpha_max)
    addParameter(spec,'Eb',0,enorm,1.0,elim_min,elim_max)
    addParameter(spec,'beta',beta_free,beta_value,1.0,beta_min,beta_max)
#    addParameter(spec,'Scale',0,2e5,1.0,elim_min,elim_max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type','SkyDirFunction')
    addParameter(spatial,'RA',0,ra,1.0,-360.0,360.0)
    addParameter(spatial,'DEC',0,dec,1.0,-90.0,90.0)
    src.appendChild(spatial)
    lib.appendChild(src)


def addPSBrokenPowerLaw2(lib, name, ra, dec, emin=200, emax=100000,
                         ebreak_free=0, ebreak=0, ebreak_min=0, ebreak_max=0,
                         flux_free=1, flux_value=1.6, flux_scale=1e-6,
                         flux_max=1000.0, flux_min=1e-5,
                         index_lo_free=1, index_lo_value=-2.0,
                         index_lo_min=-5.0, index_lo_max=-1.0,
                         index_hi_free=1, index_hi_value=-2.0,
                         index_hi_min=-5.0, index_hi_max=-1.0):
    elim_min = 30;
    elim_max = 300000;
    if emin<elim_min:
        elim_min = emin
    if emax>elim_max:
        elim_max = emax 
    if ebreak_min == 0:
        ebreak_min = emin
    if ebreak_max == 0:
        ebreak_max = emax
    if ebreak == 0:
        ebreak = math.sqrt(ebreak_min*ebreak_max)
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','PointSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','BrokenPowerLaw2')
    addParameter(spec,'Integral',
                 flux_free,flux_value,flux_scale,flux_min,flux_max)
    addParameter(spec,'Index1',
                 index_lo_free,index_lo_value,1.0,index_lo_min,index_lo_max)
    addParameter(spec,'Index2',
                 index_hi_free,index_hi_value,1.0,index_hi_min,index_hi_max)
    addParameter(spec,'BreakValue',
                 ebreak_free,ebreak,1.0,ebreak_min,ebreak_max)
    addParameter(spec,'LowerLimit',0,emin,1.0,elim_min,elim_max)
    addParameter(spec,'UpperLimit',0,emax,1.0,elim_min,elim_max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type','SkyDirFunction')
    addParameter(spatial,'RA',0,ra,1.0,-360.0,360.0)
    addParameter(spatial,'DEC',0,dec,1.0,-90.0,90.0)
    src.appendChild(spatial)
    lib.appendChild(src)


def addPSPLSuperExpCutoff(lib, name, ra, dec, eflux=0,
                   flux_free=1, flux_value=1e-9, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index1_free=1, index1_value=-2.0,
                   index1_min=-5.0, index1_max=-0.5,
                   cutoff_free=1, cutoff_value=1e4,
                   cutoff_min=200, cutoff_max=3e5,
                   index2_free=0, index2_value=-1.0,
                   index2_min=-5.0, index2_max=-0.5):
    elim_min = 30;
    elim_max = 300000;
#    if emin<elim_min:
#        elim_min = emin
#    if emax>elim_max:
#        elim_max = emax 
#    if eflux==0:
#        eflux =2e5# meanEnergy(emin,emax,index_value)
#        flux_value *= (eflux/100.0)**index_value
#    else :
#	flux_value *= 0.1*(eflux/100.0)**index_value
    if flux_scale == 0:
        flux_scale=fluxScale(flux_value)
    flux_value /= flux_scale        
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name',name)
    src.setAttribute('type','PointSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type','PLSuperExpCutoff')
    addParameter(spec,'Prefactor',
                 flux_free,flux_value,flux_scale,flux_min,flux_max)
    addParameter(spec,'Index1',index1_free,index1_value,1.0,index1_min,index1_max)
    addParameter(spec,'Scale',0,eflux,1.0,elim_min,elim_max)
    addParameter(spec,'Cutoff',cutoff_free,cutoff_value,1.0,cutoff_min,cutoff_max)
    addParameter(spec,'Index2',index2_free,index2_value,1.0,index2_min,index2_max)

    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type','SkyDirFunction')
    addParameter(spatial,'RA',0,ra,1.0,-360.0,360.0)
    addParameter(spatial,'DEC',0,dec,1.0,-90.0,90.0)
    src.appendChild(spatial)
    lib.appendChild(src)
 

def GetlistFromFits(Configuration,catalog,verbosity=1):

	srcname=Configuration['target']['name']
	ra_src=Configuration['target']['ra']
	dec_src=Configuration['target']['dec']
	emin = Configuration['energy']['emin']

	roi  = Configuration['space']['rad']
	max_radius  = Configuration['model']['max_radius']
	min_significance = Configuration['model']['min_significance']
	model = Configuration['target']['spectrum']

	cfile = pyfits.open(catalog)
	data=cfile[1].data


	names = data.field('Source_Name')
	ra    = data.field('RAJ2000')
	dec   = data.field('DEJ2000')
	flux  = data.field('Flux_Density')
	pivot = data.field('Pivot_Energy')
	index = data.field('Spectral_Index')
	cutoff = data.field('Cutoff')
	spectype = data.field('SpectrumType')
	beta  = data.field('beta')
	sigma = data.field('Signif_Avg')

	listSource = [{'name' :srcname, 'ra' : ra_src, 'dec' : dec_src, 'flux' : 1e-9, 'index' : -2, 'scale' : emin, 'cutoff' : 1e4, 'beta' : 0.1, 'IsFree' : 1 , 'SpectrumType' : model}]

	Nfree = 1
	for i in xrange(len(names)):

		r = calcAngSepDeg(float(ra[i]), float(dec[i]), ra_src, dec_src)
		if  r<max_radius and r>.1 and  sigma[i]>min_significance :
			Nfree += 1 
			listSource.append({'name' :names[i], 'ra' : ra[i], 'dec' : dec[i], 'flux' : flux[i], 'index' : -index[i], 'scale' : pivot[i], 'cutoff' : cutoff[i], 'beta' : beta[i], 'IsFree' : 1 , 'SpectrumType' : spectype[i]})
		else :
			if  r<roi and r>.1  and  sigma[i]>min_significance :
				listSource.append({'name' :names[i], 'ra' : ra[i], 'dec' : dec[i], 'flux' : flux[i], 'index' : -index[i], 'scale' : pivot[i], 'cutoff' : cutoff[i], 'beta' : beta[i], 'IsFree' : 0 , 'SpectrumType' : spectype[i]})


	print "Add ",len(listSource)," sources in the ROI of ",roi," degrees"
	print Nfree," sources have free parameters inside ",max_radius," degrees"
	return listSource


def calcAngSepDeg(ra0, dec0, ra1, dec1):
        '''Return the angular separation between two objects. Use the
        special case of the Vincenty formula that is accurate for all
        distances'''
        C = pi/180
        d0 = C*dec0
        d1 = C*dec1
        r12 = C*(ra0-ra1)
        cd0 = cos(d0)
        sd0 = sin(d0)
        cd1 = cos(d1)
        sd1 = sin(d1)
        cr12 = cos(r12)
        sr12 = sin(r12)
        num = sqrt((cd0*sr12)**2 + (cd1*sd0-sd1*cd0*cr12)**2)
        den = sd0*sd1+cd0*cd1*cr12
        return atan2(num,den)/C

def IsIn( Name,ListOfSource):
	for Src in ListOfSource:
		if Src == Name :
			return True
	return False



def WriteXml(lib,doc,srclist,Configuration):

	emin = Configuration['energy']['emin']
	emax = Configuration['energy']['emax']
	

	#test if the user profide diffuse files
	if Configuration['model']['diffuse_gal_dir'] =="" :
		Gal_dir = env.DIFFUSE_DIR
	else :
		Gal_dir = Configuration['model']['diffuse_gal_dir']

	if Configuration['model']['diffuse_iso_dir'] =="" :
		Iso_dir = env.DIFFUSE_DIR
	else :
		Iso_dir = Configuration['model']['diffuse_iso_dir']


	if Configuration['model']['diffuse_gal'] =="" :
		Gal = Gal_dir+"/"+env.DIFFUSE_GAL
	else :
		Gal = Gal_dir+"/"+Configuration['model']['diffuse_gal']

	if Configuration['model']['diffuse_gal_dir'] =="" :
		Iso = Iso_dir+"/"+env.DIFFUSE_ISO_SOURCE
	else :
		Iso =  Iso_dir+"/"+Configuration['model']['diffuse_iso']


	i=0

	addDiffusePL(lib, Iso, free=1, value=1.0, scale=1.0, max=10.0, min=1.0, name = 'iso_p7v6source.txt')

	addGalprop(lib, Gal, free=1, value=1.0, scale=1.0, max=10.0, min=.010, name = 'gal_2yearp7v6_v0')

#	listSource = [{'name' :srcname, 'ra' : ra_src, 'dec' : dec_src, 'flux' : 1e-9, 'index' : -2, 'scale' : emin, 'cutoff' : 1e4, 'beta' : 0.1, 'IsFree' : 1 , 'SpectrumType' : model}]

	for i in xrange(len(srclist)):
		name = srclist[i].get('name')
		ra = srclist[i].get('ra')
		dec = srclist[i].get('dec')
		free = srclist[i].get('IsFree')
		spectype = srclist[i].get('SpectrumType')
		if spectype == "PowerLaw" :
			addPSPowerLaw1(lib, name, ra, dec,
                                          eflux=srclist[i].get('scale'),
                                          flux_free=free,flux_value= srclist[i].get('flux'),
                                          index_free=free,index_value=srclist[i].get('index'))
		if spectype == "PowerLaw2" :
			addPSPowerLaw2(lib, name, ra, dec, 
                                        emin=emin, emax=emax,
                                        flux_free=free,flux_value= srclist[i].get('flux'),
                                        index_free=free,index_value=srclist[i].get('index'))
		if spectype == "LogParabola" :
			addPSLogparabola(lib, name, ra, dec, enorm=srclist[i].get('scale'),
                                          norm_free=free, norm_value=srclist[i].get('flux'),
                                          alpha_free=free, alpha_value=abs(srclist[i].get('index')), 
                                          beta_free=free, beta_value=srclist[i].get('beta'))
		if spectype == "PLExpCutoff" :
			addPSPLSuperExpCutoff(lib, name, ra, dec,
                                          eflux=srclist[i].get('scale'),
                                          flux_free=free,flux_value= srclist[i].get('flux'),
                                          index1_free=free,index1_value=srclist[i].get('index'),
                                          cutoff_free=free, cutoff_value=srclist[i].get('cutoff'))
#		addPSPowerLaw1(lib, name, ra, dec, eflux=srclist[i].get('scale'),flux_free=free,flux_value= srclist[i].get('flux'),index_free=free,index_value=srclist[i].get('index'))

	folder = Configuration['out']
	os.system('mkdir -p '+folder)

	output=Configuration['file']['xml']
	print "write the Xml file in ",output
	open(output,'w').write(doc.toprettyxml('  '))

def CreateLib():
	domimpl = xml.dom.minidom.getDOMImplementation()
	doc = domimpl.createDocument(None, "source_library", None)
	lib = doc.documentElement
	lib.setAttribute("title", "source library")
	lib.appendChild(doc.createComment('Source library created by %s at %s'%(sys.argv[0],time.asctime())))
	return lib,doc


def Xml_to_Reg(Filename,listSource,Prog = None) :

	fds9 = open(Filename+".reg","w")
	if not(Prog==None) :
		fds9.write('# File generated by '+Prog+'\n')
	fds9.write('# Author David Sanchez\n')
	fds9.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
	fds9.write("fk5\n")

	for src in listSource:
		ra = src.get('ra')
		dec = src.get('dec')
		name = src.get('name')
 		if src.get('IsFree') :
			color = 'blue'
		else :
			color = 'green'	  
		fds9.write("circle("+str(ra)+","+str(dec)+",0.2)  # color="+color+" text={"+str(name)+"}\n") 
	fds9.close()


if __name__=='__main__':

	try :
		infile = sys.argv[1]
	except :
		print('FATAL: Config file not found.')
		sys.exit(1)

	try :
		Catalog=sys.argv[2]
	except :
		print('FATAL: Catalog file not found.')
		sys.exit(1)

	Configuration = get_config(infile)

	folder = Configuration['out']
	os.system('mkdir -p '+folder)

	lib = CreateLib()

	srclist =GetlistFromFits(Configuration,Catalog)
#LogParabola
	#Write donw the XML file
	WriteXml(lib,doc,srclist,Configuration)
	Xml_to_Reg(folder+"/Roi_model",srclist,Prog = sys.argv[0])
