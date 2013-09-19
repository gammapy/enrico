"""Central place for the XML generation"""
import os,sys
import xml.dom.minidom
import numpy as np
import pyfits
import utils

import enrico.environ as env
from enrico.config import get_config

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def addParameter(el, name, free, value, scale, min, max):
    """Add a parameter to a source"""
    doc = el.ownerDocument
    param = doc.createElement('parameter')
    param.setAttribute('name', name)
    param.setAttribute('free', '%d' % free)
    param.setAttribute('scale', '%g' % scale)
    param.setAttribute('value', '%g' % value)
    param.setAttribute('max', '%g' % max)
    param.setAttribute('min', '%g' % min)
    el.appendChild(param)

def addDiffusePL(lib, file, free=1, value=1.0, max=10.0, min=1.0,
               name='iso_p7v6source'):
    """Add the diffuse extragalactic diffuse"""
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('file', file)
    spec.setAttribute('type', 'FileFunction')
    addParameter(spec, 'Normalization', free, 1, 1, 0.001, 1000)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('type', 'ConstantValue')
    addParameter(spatial, 'Value', 0, 1, 1.0, 0.0, 10.0)
    src.appendChild(spatial)
    lib.appendChild(src)

def addGalprop(lib, file, free=1, value=1.0, scale=1.0, max=10.0, min=.010,
               name='gal_2yearp7v6_v0'):
    """Add the diffuse galactic diffuse"""
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type', 'ConstantValue')
    addParameter(spec, 'Value', free, value, scale, min, max)
    src.appendChild(spec)
    spatial = doc.createElement('spatialModel')
    spatial.setAttribute('file', file)
    spatial.setAttribute('type', 'MapCubeFunction')
    addParameter(spatial, 'Normalization', 0, 1, 1, 0.001, 1000)
    src.appendChild(spatial)
    lib.appendChild(src)

def addPSPowerLaw1(lib, name, ra, dec, eflux=0,
                   flux_free=1, flux_value=1e-9, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index_free=1, index_value=-2.0,
                   index_min=-5.0, index_max=-0.5,extendedName=""):
    """Add a source with a POWERLAW1 model"""
    elim_min = 30
    elim_max = 300000
    if flux_scale == 0:
        flux_scale = utils.fluxScale(flux_value)
    flux_value /= flux_scale
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    if extendedName=="":
      src.setAttribute('type', 'PointSource')
    else:
      src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type', 'PowerLaw')
    addParameter(spec, 'Prefactor',
                 flux_free, flux_value, flux_scale, flux_min, flux_max)
    addParameter(spec, 'Index', index_free, index_value, 1.0,
                 index_min, index_max)
    addParameter(spec, 'Scale', 0, eflux, 1.0, elim_min, elim_max)
    src.appendChild(spec)
    spatial = AddSpatial(doc,ra,dec,extendedName)
    src.appendChild(spatial)
    lib.appendChild(src)


def addPSPowerLaw2(lib, name, ra, dec, emin=200, emax=3e5,
                   flux_free=1, flux_value=1.6e-8, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index_free=1, index_value=-2.0,
                   index_min=-5.0, index_max=-0.5,extendedName=""):
    """Add a source with a POWERLAW2 model"""
    elim_min = 30
    elim_max = 300000
    if emin < elim_min:
        elim_min = emin
    if emax > elim_max:
        elim_max = emax
    if flux_scale == 0:
        flux_scale = utils.fluxScale(flux_value)
    flux_value /= flux_scale
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    if extendedName=="":
      src.setAttribute('type', 'PointSource')
    else:
      src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type', 'PowerLaw2')
    addParameter(spec, 'Integral',
                 flux_free, flux_value, flux_scale, flux_min, flux_max)
    addParameter(spec, 'Index', index_free, index_value, 1.0,
                 index_min, index_max)
    addParameter(spec, 'LowerLimit', 0, emin, 1.0, elim_min, elim_max)
    addParameter(spec, 'UpperLimit', 0, emax, 1.0, elim_min, elim_max)
    src.appendChild(spec)
    spatial = AddSpatial(doc,ra,dec,extendedName)
    src.appendChild(spatial)
    lib.appendChild(src)


def addPSLogparabola(lib, name, ra, dec, enorm=300,
                   norm_free=1, norm_value=1e-9, norm_scale=0,
                   norm_max=1000.0, norm_min=1e-5,
                   alpha_free=1, alpha_value=1.0,
                   alpha_min=.5, alpha_max=5.,
                   beta_free=1, beta_value=1.0,
                   beta_min=0.0005, beta_max=5.0,extendedName=""):
    """Add a source with a LOGPARABOLA model"""
    elim_min = 30
    elim_max = 300000

    if enorm == 0:
        enorm = 2e5  # meanEnergy(emin,emax,index_value)
        norm_value *= (enorm / 100.0) ** alpha_value
    if norm_scale == 0:
        norm_scale = utils.fluxScale(norm_value)
    norm_value /= norm_scale
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    if extendedName=="":
      src.setAttribute('type', 'PointSource')
    else:
      src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type', 'LogParabola')
    addParameter(spec, 'norm',
                 norm_free, norm_value, norm_scale, norm_min, norm_max)
    addParameter(spec, 'alpha', alpha_free, alpha_value, 1.0,
                 alpha_min, alpha_max)
    addParameter(spec, 'Eb', 0, enorm, 1.0, elim_min, elim_max)
    addParameter(spec, 'beta', beta_free, beta_value, 1.0, beta_min, beta_max)
    src.appendChild(spec)
    spatial = AddSpatial(doc,ra,dec,extendedName)
    src.appendChild(spatial)
    lib.appendChild(src)


def addPSBrokenPowerLaw2(lib, name, ra, dec, emin=200, emax=100000,
                         ebreak_free=0, ebreak=0, ebreak_min=0, ebreak_max=0,
                         flux_free=1, flux_value=1.6, flux_scale=1e-6,
                         flux_max=1000.0, flux_min=1e-5,
                         index_lo_free=1, index_lo_value=-2.0,
                         index_lo_min=-5.0, index_lo_max=-1.0,
                         index_hi_free=1, index_hi_value=-2.0,
                         index_hi_min=-5.0, index_hi_max=-1.0,extendedName=""):
    """Add a source with a BROKENPOWERLAW2 model"""
    elim_min = 30
    elim_max = 300000
    if emin < elim_min:
        elim_min = emin
    if emax > elim_max:
        elim_max = emax
    if ebreak_min == 0:
        ebreak_min = emin
    if ebreak_max == 0:
        ebreak_max = emax
    if ebreak == 0:
        ebreak = np.sqrt(ebreak_min * ebreak_max)
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    if extendedName=="":
      src.setAttribute('type', 'PointSource')
    else:
      src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type', 'BrokenPowerLaw2')
    addParameter(spec, 'Integral',
                 flux_free, flux_value, flux_scale, flux_min, flux_max)
    addParameter(spec, 'Index1',
                 index_lo_free, index_lo_value, 1.0,
                 index_lo_min, index_lo_max)
    addParameter(spec, 'Index2',
                 index_hi_free, index_hi_value, 1.0,
                 index_hi_min, index_hi_max)
    addParameter(spec, 'BreakValue',
                 ebreak_free, ebreak, 1.0, ebreak_min, ebreak_max)
    addParameter(spec, 'LowerLimit', 0, emin, 1.0, elim_min, elim_max)
    addParameter(spec, 'UpperLimit', 0, emax, 1.0, elim_min, elim_max)
    src.appendChild(spec)
    spatial = AddSpatial(doc,ra,dec,extendedName)
    src.appendChild(spatial)
    lib.appendChild(src)


def addPSPLSuperExpCutoff(lib, name, ra, dec, eflux=0,
                   flux_free=1, flux_value=1e-9, flux_scale=0,
                   flux_max=1000.0, flux_min=1e-5,
                   index1_free=1, index1_value=-2.0,
                   index1_min=-5.0, index1_max=-0.5,
                   cutoff_free=1, cutoff_value=1e4,
                   cutoff_min=200, cutoff_max=3e5,
                   index2_free=0, index2_value=1.0,
                   index2_min=0.0, index2_max=5.0,extendedName=""):
    """Add a source with a SUPEREXPCUTOFF model"""
    elim_min = 30
    elim_max = 300000
    if flux_scale == 0:
        flux_scale = utils.fluxScale(flux_value)
    flux_value /= flux_scale
    doc = lib.ownerDocument
    src = doc.createElement('source')
    src.setAttribute('name', name)
    if extendedName=="":
      src.setAttribute('type', 'PointSource')
    else:
      src.setAttribute('type', 'DiffuseSource')
    spec = doc.createElement('spectrum')
    spec.setAttribute('type', 'PLSuperExpCutoff')
    addParameter(spec, 'Prefactor',
                 flux_free, flux_value, flux_scale, flux_min, flux_max)
    addParameter(spec, 'Index1', index1_free, index1_value, 1.0,
                 index1_min, index1_max)
    addParameter(spec, 'Scale', 0, eflux, 1.0, elim_min, elim_max)
    addParameter(spec, 'Cutoff', cutoff_free, cutoff_value, 1.0,
                 cutoff_min, cutoff_max)
    addParameter(spec, 'Index2', index2_free, index2_value, 1.0,
                 index2_min, index2_max)

    src.appendChild(spec)
    spatial = AddSpatial(doc,ra,dec,extendedName)
    src.appendChild(spatial)
    lib.appendChild(src)

def AddSpatial(doc,ra,dec,extendedName=""):
    spatial = doc.createElement('spatialModel')
    if extendedName=="":
      spatial.setAttribute('type', 'SkyDirFunction')
      addParameter(spatial, 'RA', 0, ra, 1.0, -360.0, 360.0)
      addParameter(spatial, 'DEC', 0, dec, 1.0, -90.0, 90.0)
    else : 
      from environ import CATALOG_TEMPLATE_DIR
      from os.path import join
      filename = extendedName+'.fits'
      filename = filename.replace(' ','')
      spatialModel = join(CATALOG_TEMPLATE_DIR, filename)
      spatial.setAttribute('type', 'SpatialMap')
      spatial.setAttribute('file', spatialModel)
      addParameter(spatial, 'Prefactor', 1, 1, 1.0, 0.0001,1000)

    return spatial

def GetlistFromFits(config, catalog):
    """Read the config and catalog file and generate the list of sources to include"""
    #Get the informations for the config file
    srcname = config['target']['name']
    ra_src = config['target']['ra']
    dec_src = config['target']['dec']
    ra_space = config['space']['xref']
    dec_space = config['space']['yref']
    emin = config['energy']['emin']
    roi = config['space']['rad']
    max_radius = config['model']['max_radius']
    min_significance = config['model']['min_significance']
    model = config['target']['spectrum']

    if model == "Generic":
        log.warning("Generic model found. Will turn it to PowerLaw")
        model = "PowerLaw"

    #read the catalog file
    cfile = pyfits.open(catalog)
    data = cfile[1].data
    names = data.field('Source_Name')
    ra = data.field('RAJ2000')
    dec = data.field('DEJ2000')
    flux = data.field('Flux_Density')
    pivot = data.field('Pivot_Energy')
    index = data.field('Spectral_Index')
    try :  # valid for the 2FGH, not for the 1FHL
      cutoff = data.field('Cutoff')
      beta = data.field('beta')
      spectype = data.field('SpectrumType')
    except :
      cutoff = np.zeros(names.size)
      beta = np.zeros(names.size)
      spectype = np.array(names.size*["PowerLaw"])
      pivot *= 1e3 ## energy in the 1FHL are in GeV
      flux *= 1e3
    try :
      extendedName = data.field('Extended_Source_Name')
    except:
      logging.warning("Cannot find th extended source list: please check the xml")
      extendedName = np.array(names.size*[""])

    sigma = data.field('Signif_Avg')

    sources = []
    Nfree = 0
    Nextended = 0
    #loop over all the sources of the catalog
    for i in xrange(len(names)):
        #distance from the center of the maps
        rspace = utils.calcAngSepDeg(float(ra[i]), float(dec[i]), ra_space, dec_space)
        #distance for the target
        rsrc = utils.calcAngSepDeg(float(ra[i]), float(dec[i]), ra_src, dec_src)

        # if the source has a separation less than 0.1deg to the target and has
        # the same model type as the one we want to use, insert as our target
        # with our given coordinates
        if rsrc < .1 and sigma[i] > min_significance and spectype[i] == model and extendedName[i]=="":
            Nfree += 1
            sources.insert(0,{'name': srcname, 'ra': ra_src, 'dec': dec_src,
                            'flux': flux[i], 'index': -index[i], 'scale': pivot[i],
                            'cutoff': cutoff[i], 'beta': beta[i], 'IsFree': 1,
                            'SpectrumType': spectype[i], 'ExtendedName': extendedName[i]})

        elif  rsrc < max_radius and rsrc > .1 and  sigma[i] > min_significance:
            # if the source is close to the target : add it as a free source
            Nfree += 1
            sources.append({'name': names[i], 'ra': ra[i], 'dec': dec[i],
                            'flux': flux[i], 'index': -index[i], 'scale': pivot[i],
                            'cutoff': cutoff[i], 'beta': beta[i], 'IsFree': 1,
                            'SpectrumType': spectype[i], 'ExtendedName': extendedName[i]})
            if not(extendedName[i]==""):
                print "Adding extended source ",extendedName[i]," 2FGL name is ",names[i]
                Nextended=+1
        else:
            # if the source is inside the ROI: add it as a frozen source
            if  rspace < roi and rsrc > .1  and  sigma[i] > min_significance:
                sources.append({'name': names[i], 'ra': ra[i], 'dec': dec[i],
                                'flux': flux[i], 'index': -index[i], 'scale': pivot[i],
                                'cutoff': cutoff[i], 'beta': beta[i], 'IsFree': 0,
                                'SpectrumType': spectype[i],'ExtendedName': extendedName[i]})
                if not(extendedName[i]==""):
                  print "Adding extended source ",extendedName[i]," 2FGL name is ",names[i]
                  Nextended=+1


    # if the target has not been added from catalog, add it now
    if sources[0]['name']!=srcname:
        Nfree += 1
        #add the target to the list of sources in first position
        sources.insert(0,{'name':srcname, 'ra': ra_src, 'dec': dec_src,
                       'flux': 1e-9, 'index':-2, 'scale': emin,
                       'cutoff': 1e4, 'beta': 0.1, 'IsFree': 1,
                       'SpectrumType': model,'ExtendedName': ""})

    print "Add ", len(sources), " sources in the ROI of ", roi, " degrees"
    print Nfree, " sources have free parameters inside ", max_radius, " degrees"
    print Nextended, " source(s) is (are) extended"
    return sources

def IsIn(name, sources):
    for source in sources:
        if source == name:
            return True
    return False

def WriteXml(lib, doc, srclist, config):
    """Fill and write the library of sources into an XML file"""
    emin = config['energy']['emin']
    emax = config['energy']['emax']

    Galname = "GalDiffModel"
    Isoname = "IsoDiffModel"

    #test if the user provides diffuse files. if not  use the default one
    if config['model']['diffuse_gal_dir'] == "":
        Gal_dir = env.DIFFUSE_DIR
    else:
        Gal_dir = config['model']['diffuse_gal_dir']

    if config['model']['diffuse_iso_dir'] == "":
        Iso_dir = env.DIFFUSE_DIR
    else:
        Iso_dir = config['model']['diffuse_iso_dir']

    if config['model']['diffuse_gal'] == "":
        Gal = Gal_dir + "/" + env.DIFFUSE_GAL
    else:
        Gal = Gal_dir + "/" + config['model']['diffuse_gal']

    if config['model']['diffuse_iso'] == "":
        Iso = Iso_dir + "/" + env.DIFFUSE_ISO_SOURCE

        if config['analysis']['irfs'] == "P7CLEAN_V6":# and config['analysis']['convtype'] == -1:
            Iso = Iso_dir + "/" + env.DIFFUSE_ISO_CLEAN

    else:
        Iso = Iso_dir + "/" + config['model']['diffuse_iso']

    #add diffuse sources
    addDiffusePL(lib, Iso, free=1, value=1.0,
                 max=10.0, min=1.0, name=Isoname)
    addGalprop(lib, Gal, free=1, value=1.0, scale=1.0,
               max=10.0, min=.010, name=Galname)

    # loop over the list of sources and add it to the library
    for i in xrange(len(srclist)):
        name = srclist[i].get('name')
        ra = srclist[i].get('ra')
        dec = srclist[i].get('dec')
        free = srclist[i].get('IsFree')
        spectype = srclist[i].get('SpectrumType')
        extendedName = srclist[i].get('ExtendedName')
        # Check the spectrum model
        if spectype == "PowerLaw":
            addPSPowerLaw1(lib, name, ra, dec,
                              eflux=srclist[i].get('scale'),
                              flux_free=free, flux_value=srclist[i].get('flux'),
                              index_free=free, index_value=srclist[i].get('index'),extendedName=extendedName)
        if spectype == "PowerLaw2":
            addPSPowerLaw2(lib, name, ra, dec,
                            emin=emin, emax=emax,
                            flux_free=free, flux_value=srclist[i].get('flux'),
                            index_free=free, index_value=srclist[i].get('index'),extendedName=extendedName)
        if spectype == "LogParabola":
            addPSLogparabola(lib, name, ra, dec, enorm=srclist[i].get('scale'),
                              norm_free=free, norm_value=srclist[i].get('flux'),
                              alpha_free=free, alpha_value=abs(srclist[i].get('index')),
                              beta_free=free, beta_value=srclist[i].get('beta'),extendedName=extendedName)
        if spectype == "PLExpCutoff":
            addPSPLSuperExpCutoff(lib, name, ra, dec,
                              eflux=srclist[i].get('scale'),
                              flux_free=free, flux_value=srclist[i].get('flux'),
                              index1_free=free, index1_value=srclist[i].get('index'),
                              cutoff_free=free, cutoff_value=srclist[i].get('cutoff'),extendedName=extendedName)

    folder = config['out']
    os.system('mkdir -p ' + folder)

    output = config['file']['xml']
    print "write the Xml file in ", output
    open(output, 'w').write(doc.toprettyxml('  '))#save it


def CreateLib():
    """@todo: document me"""
    import sys
    import time
    domimpl = xml.dom.minidom.getDOMImplementation()
    doc = domimpl.createDocument(None, "source_library", None)
    lib = doc.documentElement
    lib.setAttribute("title", "source library")
    lib.appendChild(doc.createComment('Source library created by %s at %s' %
                                      (sys.argv[0], time.asctime())))
    return lib, doc


def Xml_to_Reg(Filename, listSource, Prog=None):
    """Convert model from xml format to ds9 format"""
    fds9 = open(Filename + ".reg", "w")
    if Prog != None:
        fds9.write('# File generated by ' + Prog + '\n')
    fds9.write('global color=green dashlist=8 3 width=1 '
               'font="helvetica 10 normal roman" select=1 highlite=1 '
               'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    fds9.write("fk5\n")

    for src in listSource:
        ra = src.get('ra')
        dec = src.get('dec')
        name = src.get('name')
        if src.get('IsFree'):
            color = 'blue'
        else:
            color = 'green'
        fds9.write("circle(" + str(ra) + "," + str(dec) + ",0.2)  # color=" +
                   color + " text={" + str(name) + "}\n")
    fds9.close()


def XmlMaker(config):
  folder = config['out']
  os.system('mkdir -p ' + folder)

# test if the user provide a catalog or not.
#if not use the default one
  if config['environ']['FERMI_CATALOG_DIR'] == '':
    catalogDir = env.CATALOG_DIR
    print "use the default location of the catalog"
  else:
    catalogDir = config['environ']['FERMI_CATALOG_DIR']

  if config['environ']['FERMI_CATALOG'] == '':
    catalog = catalogDir + "/" + env.CATALOG
    print "use the default catalog"
  else:
    catalog = catalogDir + "/" + config['environ']['FERMI_CATALOG']

  print "Use the catalog : ", catalog

  lib, doc = CreateLib()
  srclist = GetlistFromFits(config, catalog)
  WriteXml(lib, doc, srclist, config)
  Xml_to_Reg(folder + "/Roi_model",
                            srclist, Prog=sys.argv[0])
