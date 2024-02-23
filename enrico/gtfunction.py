"""
gtfunction.py written by David Sanchez : david.sanchez@lapp.in2p3.fr
Collection of function to run the ST tools.
The description of each ST tool is given in the dedicated NASA website
The Observation class contains all the variables needed  to run the ST like file path, energy, time etc, ...
begun October 2010
"""
import os
import sys
import glob
import shutil
from enrico import Loggin
from time import sleep
from random import random
from math import sqrt, log10
#from gt_apps import evtbin, maketime, diffResps, expCube, expMap, srcMaps, model_map, filter
from gt_apps import evtbin, maketime, diffResps, expCube, expMap, srcMaps, model_map, filter
from GtApp import GtApp
from enrico import utils
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as aascii
def run_retry(macro,tries=5,compress=False):
    """
    The Fermi LAT sometimes fail with annoying Runtime errors,
    try to catch them and re-run the macro that failed. 
    Do that up to 5 times with random waiting times.
    """
    mes = Loggin.Message()

    # Try to write the temporary output to a temporary file and then move it, 
    # to avoid broken files left all over the place if the macro is interrupted
    try :
        orig_name = str(macro['outfile']).replace(".gz","")
    except:
        orig_name = str(macro['infile']).replace(".gz","") #for gtexposure which has no outfile

    try:
        macro['outfile'] = macro['outfile']+".tmpout"
    except:
        is_out_in_tmp = False
    else:
        is_out_in_tmp = True

    for retry in range(tries):
        try:
            macro.run()
        except RuntimeError:
            # Wait between 10 and 20 seconds to try again
            mes.warning("An error ocurred, retrying ...")
            sleep(10.*random()+10.)
            continue
        else:
            if is_out_in_tmp:
                shutil.move(macro['outfile'],orig_name)

            # Compress the output if needed and the file exists
            if os.path.isfile(orig_name):
                if not utils.is_gz_file(orig_name):
                    if compress:
                        cmd = "gzip -f "+orig_name 
                        print('Compressing file: '+ cmd)
                        os.system(cmd)

            return(macro)
    mes.error("An error ocurred and could not be recovered. Exiting!")
    sys.exit(1)    

class Observation:
    # init function of the Observation class.
    # folder : folder where the produced fits files will be stored.
    # configuration is the config of enrico (contains the variable)

    def __init__(self,folder,Configuration,tag=""):
        self.Configuration = Configuration
        self.tag = tag
        self.folder = folder
        self.LoadConfiguration()

    def LoadConfiguration(self):
        #Read the configuration object and init all the variable
        filetag = self.Configuration['file']['tag']
        self.inttag  = "_"+filetag
        if not(self.tag==""):
            self.inttag+="_"+self.tag

        self.srcname   = self.Configuration['target']['name']
        self.modelname = self.Configuration['target']['spectrum']
        self.ft1       = self.Configuration['file']['event']
        self.ft2       = self.Configuration['file']['spacecraft']
        self.xmlfile   = self.Configuration['file']['xml']

        self.SimXmlfile = self.Configuration['ObservationSimulation']['infile'] 
        self.srcList    = self.Configuration['ObservationSimulation']['srclist'] 

        #Fits files and optional gz compression
        if self.Configuration['file']['compress_fits'] == "yes":
            self.gzflag=".gz"
        else:
            self.gzflag=""

        self.eventcoarse   = self.folder+'/'+self.srcname+"_"+filetag+"_EvtCoarse.fits"+self.gzflag
        self.eventfile     = self.folder+'/'+self.srcname+self.inttag+"_Evt.fits"+self.gzflag
        self.mktimefile    = self.folder+'/'+self.srcname+self.inttag+"_MkTime.fits"+self.gzflag
        self.Cubename      = self.folder+'/'+self.srcname+self.inttag+"_ltCube.fits"+self.gzflag
        self.Mapname       = self.folder+'/'+self.srcname+self.inttag+"_ExpMap.fits"+self.gzflag
        self.BinnedMapfile = self.folder+'/'+self.srcname+self.inttag+"_BinnedMap.fits"+self.gzflag
        self.cmapfile      = self.folder+'/'+self.srcname+self.inttag+"_CountMap.fits"+self.gzflag
        self.lcfile        = self.folder+'/'+self.srcname+self.inttag+"_applc.fits"+self.gzflag
        self.ccube         = self.folder+'/'+self.srcname+self.inttag+"_CCUBE.fits"+self.gzflag
        self.srcMap        = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_srcMap.fits"+self.gzflag
        self.ModelMapFile  = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_ModelMap.fits"+self.gzflag
        self.BinDef        = self.folder+'/'+self.srcname+self.inttag+"_BinDef.fits"+self.gzflag
        self.Probfile      = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_prob.fits"+self.gzflag
        self.psf           = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_psf.fits"+self.gzflag
        self.rel_diff_file = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_ResidualMap.fits"+self.gzflag
        self.abs_diff_file = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_SubtractMap.fits"+self.gzflag
        self.drmfile       = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_eDRM.fits"+self.gzflag
        self.effbkgfile    = self.folder+'/'+self.srcname+self.inttag+"_"+self.modelname+"_effbkgfile.fits"+self.gzflag
        self.alphabkgfile  = self.folder+'/'+self.srcname+"_"+self.modelname+"_alphabkgfile.fits"+self.gzflag
        self.wtsmapfile    = self.folder+'/'+self.srcname+"_"+self.modelname+"_wtsmapfile.fits"+self.gzflag
        self.gtifitsfile   = self.folder+'/'+self.srcname+"_"+filetag+"_GTI.fits"

        #Variables
        print("DAVID :",self.Configuration['time']['type'])
        if ('MJD' in self.Configuration['time']['type']):
            get_met = lambda t: utils.MJD_to_met(float(t))
        elif ('JD' in self.Configuration['time']['type']):
            get_met = lambda t: utils.JD_to_met(float(t))
        else:
            get_met = lambda t: float(t)
        
        #use energy dispersion corrections? This will extend Emin and Emax
        self.use_edisp = bool(self.Configuration['analysis']['EnergyDispersion']=='yes' and self.Configuration["analysis"]["likelihood"]=="binned")
        
        self.t1        = get_met(self.Configuration['time']['tmin'])
        print("DAVID :", self.Configuration['time']['tmin'])
        print("DAVID :", get_met(self.Configuration['time']['tmin']))
        self.t2        = get_met(self.Configuration['time']['tmax'])
        self.Emin      = float(self.Configuration['energy']['emin'])
        self.Emax      = float(self.Configuration['energy']['emax'])
        if (self.use_edisp):
            self.Emin_ext  = 10**(log10(self.Emin)-0.3)
            self.Emax_ext  = 10**(log10(self.Emax)+0.3)
        else:
            self.Emin_ext  = self.Emin
            self.Emax_ext  = self.Emax

        self.ra        = float(self.Configuration['space']['xref'])
        self.dec       = float(self.Configuration['space']['yref'])
        self.roi       = float(self.Configuration['space']['rad'])
        self.irfs,_    = utils.GetIRFS(self.Configuration['event']['evclass'],self.Configuration['event']['evtype'])
        #self.irfs      = self.irfs
        self.likelihood = self.Configuration['analysis']['likelihood']

        #Apply cuts in event selections? 
        # (roicuts should not be applied twice, it makes ST crash)
        self.roicuts   = bool(self.Configuration['analysis']['evtroicuts']=='yes')
        self.timecuts  = bool(self.Configuration['analysis']['evttimecuts']=='yes')


        #diffuse Response
        self.diffrspflag = self.folder+'/'+self.srcname+self.inttag+"_diffrsp.flag"

        #Maps binning
        self.binsz     = self.Configuration['space']['binsz']
        self.npix      = int(2*self.roi/self.binsz)
        self.npixCntMp = int(sqrt(2.)*self.roi/self.binsz)

        #tool options
        self.clobber = self.Configuration['clobber']

    def run_retry_compress(self,macro,tries=5):
        compress = self.Configuration['file']['compress_fits'] == "yes"
        run_retry(macro,tries,compress)

    def printSum(self):
        """Print a summary of the value stored in the class"""
        print("Source = ",self.srcname)
        print("RA\t=\t",self.ra," degrees")
        print("Dec\t=\t",self.dec," degrees")
        print("Start\t=\t",self.t1,"  MET (s)")
        print("Stop\t=\t",self.t2,"  MET (s)")
        print("ROI\t=\t",self.roi," degrees")
        print("E min\t=\t",self.Emin," MeV")
        print("E max\t=\t",self.Emax," MeV")
        print("E min ext\t=\t",self.Emin_ext," MeV")
        print("E max ext\t=\t",self.Emax_ext," MeV")
        print("IRFs\t=\t",self.irfs)
        print("evclass\t=\t",self.Configuration['event']['evclass'])
        print("evtype\t=\t",self.Configuration['event']['evtype'])
        if  self.irfs == 'CALDB':
            print(("Corresponding IRFs\t=\t",\
            utils.GetIRFS(self.Configuration['event']['evclass'],\
            self.Configuration['event']['evtype'])))

    def Gtbin(self):
        """Run gtbin with the CMAP option. A square count map is produced enclosed inside the roi"""
        if (self.clobber=="no" and os.path.isfile(self.cmapfile)):
            #print("File exists and clobber is False")
            return(0)
        evtbin['evfile'] = self.mktimefile
        evtbin['scfile'] = self.ft2
        evtbin['outfile'] = self.cmapfile
        evtbin['algorithm'] = "CMAP"
        evtbin['nxpix'] = self.npixCntMp
        evtbin['nypix'] = self.npixCntMp
        evtbin['binsz'] = self.binsz
        evtbin['coordsys'] = self.Configuration['space']['coordsys']
        evtbin["emin"] = self.Emin_ext
        evtbin["emax"] = self.Emax_ext
        evtbin['xref'] = self.ra
        evtbin['yref'] = self.dec
        evtbin['axisrot'] = 0
        evtbin['proj'] = self.Configuration['space']['proj']
        evtbin['clobber'] = self.clobber
        #evtbin.run()
        self.run_retry_compress(evtbin)
    
    def GtBinDef(self,filename):
        if (self.clobber=="no" and os.path.isfile(self.BinDef)):
            #print("File exists and clobber is False")
            return(0)
        bindef = GtApp('gtbindef', 'Likelihood')
        bindef['bintype'] = 'T'
        bindef['binfile'] = filename
        bindef['outfile'] = self.BinDef
        #bindef.run()
        self.run_retry_compress(bindef)

    def GtExposure(self):
        exposure = GtApp('gtexposure', 'Likelihood')
        exposure['infile'] = self.lcfile
        exposure['scfile'] = self.ft2
        exposure['target'] = self.srcname
        #if  self.irfs != 'CALDB':
        #exposure['evtype']= self.Configuration['event']['evtype']
        exposure['irfs'] = self.irfs
        exposure['srcmdl'] = "none"
        exposure['specin'] = -self.Configuration['AppLC']['index']
        exposure['clobber'] = self.clobber
        #exposure.run()
        self.run_retry_compress(exposure)

    def GtLCbin(self,dt=60):
        """Run gtbin with the LC option. the default dt is 60 sec and data can be rebinned after.
          Can also take a file as input to define the time bins"""
        if (self.clobber=="no" and os.path.isfile(self.lcfile)):
            #print("File exists and clobber is False")
            return(0)
        evtbin['evfile'] = self.mktimefile
        evtbin['scfile'] = self.ft2
        evtbin['outfile'] = self.lcfile
        evtbin['algorithm'] = "LC"
        evtbin["tstart"] = self.t1
        evtbin["tstop"] = self.t2
        evtbin["dtime"] = dt
        if dt > 0 :
            evtbin["tbinalg"] = "LIN"
        else :
            evtbin["tbinalg"] = "FILE"
            evtbin["tbinfile"] = self.BinDef
        evtbin['clobber'] = self.clobber
        #evtbin.run()
        self.run_retry_compress(evtbin)

    def GtCcube(self):
        """Run gtbin with the CCUBE option"""
        if (self.clobber=="no" and os.path.isfile(self.ccube)):
            #print("File exists and clobber is False")
            return(0)
        Nbdecade = log10(self.Emax_ext)-log10(self.Emin_ext)#Compute the number of decade
        evtbin['evfile'] = self.mktimefile
        evtbin['scfile'] = self.ft2
        evtbin['outfile'] = self.ccube
        evtbin['algorithm'] = "CCUBE"
        evtbin['nxpix'] = self.npixCntMp
        evtbin['nypix'] = self.npixCntMp
        evtbin['binsz'] = self.binsz
        evtbin['coordsys'] = self.Configuration['space']['coordsys']
        evtbin['xref'] = self.ra
        evtbin['yref'] = self.dec
        evtbin["emin"] = self.Emin_ext
        evtbin["emax"] = self.Emax_ext
        evtbin["tstart"] = self.t1
        evtbin["tstop"] = self.t2
        evtbin['ebinalg'] = "LOG"
        evtbin['axisrot'] = 0
        evtbin['proj'] = self.Configuration['space']['proj'] #"AIT"
        #The number of bin is the number of decade * the number of bin
        #per decade (given by the users). The +0.5 rounds it properly
        evtbin["enumbins"] = max(2,int(Nbdecade*self.Configuration['energy']['enumbins_per_decade']+0.5))
        evtbin['clobber'] = self.clobber
        #evtbin.run()
        self.run_retry_compress(evtbin)

    def GtBinnedMap(self):
        """Run the gtexpcube2 tool for binned analysis"""
        if (self.clobber=="no" and os.path.isfile(self.BinnedMapfile)):
            #print("File exists and clobber is False")
            return(0)
        Nbdecade = log10(self.Emax_ext)-log10(self.Emin_ext)#Compute the number of decade
        expcube2 = GtApp('gtexpcube2', 'Likelihood')
        expcube2['infile'] = self.Cubename
        expcube2['outfile'] = self.BinnedMapfile
        expcube2['cmap'] = self.ccube
        #if  self.irfs != 'CALDB':
        expcube2['evtype']= self.Configuration['event']['evtype']
        expcube2['irfs'] = self.irfs
        expcube2['emin'] = self.Emin_ext
        expcube2['emax'] = self.Emax_ext
        expcube2['xref'] = "INDEF"
        expcube2['yref'] = "INDEF"
        expcube2['nxpix'] = "INDEF"
        expcube2['nypix'] = "INDEF"
        expcube2['binsz'] = "INDEF"
        app = expcube2
        if 'edisp_bins' in list(app.pars.keys()):
            if self.use_edisp:
                app['edisp_bins'] = -min(3,int(Nbdecade*0.2+0.5))
            else:
                app['edisp_bins'] = 0
        elif 'edisp' in list(app.pars.keys()):
            app['edisp'] = True
                
        expcube2['enumbins'] = max(2,int(Nbdecade*self.Configuration['energy']['enumbins_per_decade']+0.5))
        expcube2['coordsys'] = self.Configuration['space']['coordsys']
        expcube2['proj'] = self.Configuration['space']['proj'] #"AIT"
        expcube2['clobber'] = self.clobber
        #expcube2.run()
        self.run_retry_compress(expcube2)

    def FirstCut(self):
        """Run gtselect tool"""
        if (self.clobber=="no" and os.path.isfile(self.eventcoarse)):
            #print("File exists and clobber is False")
            return(0)
        filter['infile'] = self.ft1
        filter['outfile'] = self.eventcoarse
        if (self.roicuts == True):
            filter['ra'] = self.ra
            filter['dec'] = self.dec
            filter['rad'] = self.roi
        else:
            filter['ra']  = 0
            filter['dec'] = 0
            filter['rad'] = 180
        if (self.timecuts == True):
            filter['tmin'] = self.t1
            filter['tmax'] = self.t2
        else:
            filter['tmin'] = "INDEF"
            filter['tmax'] = "INDEF"
        filter['emin'] = 0        #self.Emin Do not cut based on energy for now
        filter['emax'] = 10000000 #self.Emax
        filter['zmax'] = self.Configuration['analysis']['zmax']
        filter['evclass'] = self.Configuration['event']['evclass']
        filter['evtype'] = "INDEF"
        filter['clobber'] = self.clobber
        #filter.run()
        self.run_retry_compress(filter)

    def SelectEvents(self):
        """Run gtselect tool"""
        if (self.clobber=="no" and os.path.isfile(self.eventfile)):
            #print("File exists and clobber is False")
            return(0)
        filter['infile'] = self.eventcoarse
        filter['outfile'] = self.eventfile
        filter['ra'] =   0            #self.ra
        filter['dec'] =  0            #self.dec
        filter['rad'] =  180          #self.roi
        filter['tmin'] = "INDEF"      #self.t1
        filter['tmax'] = "INDEF"      #self.t2
        filter['emin'] = self.Emin_ext
        filter['emax'] = self.Emax_ext
        filter['zmax'] = self.Configuration['analysis']['zmax']
        filter['evclass'] = self.Configuration['event']['evclass']
        filter['evtype'] = self.Configuration['event']['evtype']
        filter['clobber'] = self.clobber
        #filter.run()
        self.run_retry_compress(filter)

    def gen_filter_fits_file(self):
        # Convert any set of time cuts into a fits file with the list there (as in 4FGL)
        
        data=None
        out = self.gtifitsfile
        for ext in ['.fits','.fts','fit']:
            if ext in self.Configuration['time']['file']:
                data = pyfits.open(self.Configuration['time']['file'])
                c1 = list(data['GTI'].data['START'])
                c2 = list(data['GTI'].data['STOP'])
                break
        
        if data==None:
            data = aascii.read(self.Configuration['time']['file'],names=['START','STOP'])
            c1 = np.asarray(list(data['START']))
            c2 = np.asarray(list(data['STOP']))

        if self.Configuration['time']['type']=='MJD':
            c1 = utils.MJD_to_met(c1)
            c2 = utils.MJD_to_met(c2)
        elif self.Configuration['time']['type']=='JD':
            c1 = utils.MJD_to_met(c1)
            c2 = utils.MJD_to_met(c2)

        header = pyfits.Header()
        header['COMMENT'] = "Fermi-LAT/Enrico GTI file"
        primary = pyfits.PrimaryHDU(header=header)
        c1 = pyfits.Column(name='START', array=np.array(c1), format='D')
        c2 = pyfits.Column(name='STOP',  array=np.array(c2), format='D')
        gtitab = pyfits.BinTableHDU.from_columns([c1,c2],name="GTI")
        hdu = pyfits.HDUList([primary,gtitab])
        hdu.writeto(out,overwrite=True)

    def MkTime(self):
        import os.path
        """compute GTI"""
        ## Maketime does not listen to clobber variables
        if (self.clobber=="no" and os.path.isfile(self.mktimefile)):
            return(0)

        selstr = self.Configuration['analysis']['filter']
        if self.Configuration['time']['file'] != '':
            self.gen_filter_fits_file()
            selstr = "{0} && gtifilter(\'{1}[GTI]\',START) && gtifilter(\'{1}[GTI]\',STOP)".format(selstr,self.gtifitsfile)
        
        outfile = self.mktimefile#+".tmp"
        self._RunMktime(selstr,outfile,self.Configuration['analysis']['roicut'])
        #os.system("mv "+outfile+" "+self.mktimefile)

    def _RunMktime(self,selstr,outfile,roicut):
        """run gtmktime tool"""
        if (self.clobber=="no" and os.path.isfile(self.mktimefile)):
            #print("File exists and clobber is False")
            return(0)
        maketime['scfile']  = self.ft2
        maketime['filter']  = selstr #self.Configuration['analysis']['filter']
        maketime['roicut']  = roicut
        maketime['tstart']  = self.t1
        maketime['tstop']   = self.t2
        maketime['evfile']  = self.eventfile
        maketime['outfile'] = outfile
        maketime['clobber'] = self.clobber
        #maketime.run()
        self.run_retry_compress(maketime)

    def DiffResps(self):
        """run gtdiffresp"""
        if (self.clobber=="no" and os.path.isfile(self.diffrspflag)):
            #print("File exists and clobber is False")
            return(0)
        diffResps['evfile']=self.mktimefile
        diffResps['scfile']=self.ft2
        diffResps['srcmdl']=self.xmlfile
        diffResps['evtype']= 'INDEF'
        diffResps['evclass']= 'INDEF'
        diffResps['irfs']=self.irfs
        diffResps['convert']="no"

        diffResps['clobber'] = self.clobber
        #diffResps.run()
        self.run_retry_compress(diffResps)
        with open(self.diffrspflag,"w") as diffrspflag:
            diffrspflag.write("")

        print("\ndone")

    def ExpCube(self):
        "Run gtltcube tool to produce livetime cube"
        if (self.clobber=="no" and os.path.isfile(self.Cubename)):
            #print("File exists and clobber is False")
            return(0)
        expCube['evfile']=self.mktimefile
        expCube['scfile']=self.ft2.lstrip('@') # @ allows for weekly SC files
        expCube['outfile']=self.Cubename
        expCube['dcostheta']=0.025
        expCube['binsz']=1
        expCube['zmax']=self.Configuration['analysis']['zmax']
        expCube['phibins']=self.Configuration['space']['phibins']
        expCube['clobber'] = self.clobber
        #expCube.run()
        self.run_retry_compress(expCube)

    def ExpMap(self):
        "Run gtexpmap for unbinned analysis"
        if (self.clobber=="no" and os.path.isfile(self.Mapname)):
            #print("File exists and clobber is False")
            return(0)
        Nbdecade = log10(self.Emax_ext)-log10(self.Emin_ext)#Compute the number of decade
        expMap['evfile'] = self.mktimefile
        expMap['scfile'] = self.ft2
        expMap['expcube'] = self.Cubename
        expMap['outfile'] = self.Mapname
        if  self.irfs != 'CALDB':
            expMap['evtype'] = self.Configuration['event']['evtype']
        else :
            expMap['evtype']= 'INDEF'
        expMap['irfs'] = self.irfs
        expMap['evtype'] = self.Configuration['event']['evtype']
        expMap['srcrad'] = self.roi+15
        #The number of bin is the number of decade * the number of bin per decade (given by the users)
        expMap['nenergies'] =  max(2,int(Nbdecade*self.Configuration['energy']['enumbins_per_decade']+0.5))
        expMap['clobber'] = self.clobber
        #expMap.run()
        self.run_retry_compress(expMap)

    def Obssim(self):
        obsSim = GtApp('gtobssim', 'observationSim')
        """Run gtobssim tool"""
        if (self.clobber=="no" and os.path.isfile(self.srcMap)):
            #print("File exists and clobber is False")
            return(0) 
        obsSim["infile"] = self.SimXmlfile     
        obsSim["srclist"] = self.srcList     
        obsSim['scfile'] = self.ft2
        obsSim['evroot'] = self.Configuration['file']['tag']
        obsSim['simtime'] = self.t2-self.t1
        obsSim['tstart'] = self.t1
        obsSim['use_ac'] = "no"
        obsSim["emin"]    = self.Emin_ext
        obsSim["emax"]    = self.Emax_ext
        obsSim["edisp"]    = "yes"
        # if  self.irfs != 'CALDB':
        #     obsSim['evtype']= self.Configuration['event']['evtype']
        # else :
        #     obsSim['evtype']= 'INDEF'
        obsSim['irfs']= self.irfs
        obsSim['seed']= int(random()*100000)
        self.run_retry_compress(obsSim)

    def SrcMap(self):
        """Run gtsrcmaps tool for binned analysis"""
        if (self.clobber=="no" and os.path.isfile(self.srcMap)):
            #print("File exists and clobber is False")
            return(0)
        Nbdecade = log10(self.Emax_ext)-log10(self.Emin_ext)#Compute the number of decade
        srcMaps['scfile'] = self.ft2
        srcMaps['expcube'] = self.Cubename
        srcMaps['cmap'] = self.ccube
        srcMaps['bexpmap'] = self.BinnedMapfile
        srcMaps['srcmdl']=self.xmlfile
        if  self.irfs != 'CALDB':
            srcMaps['evtype']= self.Configuration['event']['evtype']
        else :
            srcMaps['evtype']= 'INDEF'
        srcMaps['irfs']= self.irfs
        srcMaps['outfile'] = self.srcMap
        
        # energy dispersion correction
        app = srcMaps
        if 'edisp_bins' in list(app.pars.keys()):
            if self.use_edisp:
                #app['edisp_ebins'] = -min(3,int(Nbdecade*0.2+0.5)) # adaptative.
                app['edisp_bins'] = -2 # lets keep it simple
            else:
                app['edisp_bins'] = 0
        elif 'edisp' in list(app.pars.keys()):
            app['edisp'] = True
        
        srcMaps['emapbnds']='no'
        if (self.Configuration['analysis']['keep_all_srcmaps'] == 'yes'):
            # should speed up future re-fitting, at the cost of disk space
            srcMaps['copyall']='yes' 
        else:
            # default behavior, compute them on the fly.
            srcMaps['copyall']='no' 
        ### TODO: test this flag to see, speed up? (default: disabled)
        srcMaps['clobber'] = self.clobber
        #srcMaps.run()
        self.run_retry_compress(srcMaps)

    def ModelMap(self,xml):
        """Run gtmodel tool for binned analysis and make a subtraction of the produced map
         with the count map to produce a residual map"""
        if (self.clobber=="no" and os.path.isfile(self.ModelMapFile)):
            #print("File exists and clobber is False")
            return(0)
        model_map['expcube'] = self.Cubename
        model_map['srcmaps'] = self.srcMap
        model_map['bexpmap'] = self.BinnedMapfile
        model_map['srcmdl'] = xml
        model_map["irfs"] = self.irfs
        model_map['outfile'] = self.ModelMapFile
        model_map['clobber'] = self.clobber
        self.run_retry_compress(model_map)
        #Compute the residual map
        utils.SubtractFits(self.cmapfile,
                           self.ModelMapFile,
                           self.Configuration,
                           tag=self.inttag,
                           rel_diff_file=self.rel_diff_file,
                           abs_diff_file=self.abs_diff_file)

    def FindSource(self):
        """Run the gtfindsrc tool"""
        outfile = utils._dump_findsrcout(self.Configuration)
        if (self.clobber=="no" and os.path.isfile(outfile)):
            #print("File exists and clobber is False")
            return(0)
        findsrc = GtApp('gtfindsrc', 'Likelihood')
        findsrc['evfile'] = self.mktimefile
        findsrc['scfile'] = self.ft2
        if  self.irfs != 'CALDB':
            findsrc['evtype']= self.Configuration['event']['evtype']
        else :
            findsrc['evtype']= 'INDEF'
        findsrc['irfs'] = self.irfs
        findsrc['expcube'] = self.Cubename
        findsrc['expmap'] = self.Mapname
        findsrc['srcmdl'] = utils._dump_xml( self.Configuration)
        findsrc['coordsys'] = self.Configuration['space']['coordsys']
        findsrc['target'] = self.srcname
        findsrc['optimizer'] = self.Configuration["fitting"]["optimizer"]
        findsrc['ftol'] = self.Configuration["fitting"]["ftol"]
        findsrc['clobber'] = self.clobber
        findsrc['reopt'] = self.Configuration["findsrc"]["Refit"]
        findsrc['outfile'] = outfile
        #findsrc.run()
        self.run_retry_compress(findsrc)

    def SrcProb(self):
        """Run the gtsrcprob tool"""
        if (self.clobber=="no" and os.path.isfile(self.Probfile)):
            #print("File exists and clobber is False")
            return(0)
        srcprob = GtApp('gtsrcprob', 'Likelihood')
        srcprob['evfile'] = self.mktimefile
        srcprob['scfile'] = self.ft2
        if  self.irfs != 'CALDB':
            srcprob['evtype']= self.Configuration['event']['evtype']
        else :
            srcprob['evtype']= 'INDEF'
        srcprob['irfs'] = self.irfs
        srcprob['srcmdl'] = self.xmlfile
        srcprob['outfile'] = self.Probfile
        srcprob['srclist'] = self.Configuration['srcprob']['srclist']
        srcprob['clobber'] = self.clobber
        #srcprob.run()
        self.run_retry_compress(srcprob)

    def GtPSF(self):
        if (self.clobber=="no" and os.path.isfile(self.psf)):
            #print("File exists and clobber is False")
            return(0)

        #Compute the number of decade
        Nbdecade = log10(self.Emax_ext)-log10(self.Emin_ext)
        irfs,_ = utils.GetIRFS(self.Configuration['event']['evclass'],
                               self.Configuration['event']['evtype'])
        psf = GtApp('gtpsf', 'Likelihood')
        psf["expcube"] = self.Cubename
        psf["outfile"] = self.psf
        psf["irfs"]    = irfs
        psf["evtype"]  = self.Configuration['event']['evtype']
        psf["ra"]      = self.ra
        psf["dec"]     = self.dec
        psf["emin"]    = self.Emin_ext
        psf["emax"]    = self.Emax_ext
        psf["nenergies"] = max(2,\
          int(Nbdecade*self.Configuration['energy']['enumbins_per_decade']+0.5))
        psf["thetamax"] = 5.
        #psf.run()
        self.run_retry_compress(psf)

    def GtDRM(self):
        if (self.clobber=="no" and os.path.isfile(self.drmfile)):
            #print("File exists and clobber is False")
            return(0)
        irfs,_ = utils.GetIRFS(self.Configuration['event']['evclass'],
                               self.Configuration['event']['evtype'])
        drm = GtApp('gtdrm', 'DRM')
        drm["cmap"]    = self.srcMap
        drm["outfile"] = self.drmfile
        drm["irfs"]    = irfs
        drm["expcube"] = self.Cubename
        drm["bexpmap"] = self.BinnedMapfile
        drm["chatter"] = 0
        self.run_retry_compress(drm)
    
    def GtEffBkg(self,efact=2):
        '''
        DataCube of Effective background for a point source
        '''
        if (self.clobber=="no" and os.path.isfile(self.effbkgfile)):
            #print("File exists and clobber is False")
            return(0)
        irfs,_ = utils.GetIRFS(self.Configuration['event']['evclass'],
                               self.Configuration['event']['evtype'])
        effbkg = GtApp('gteffbkg', 'EffBkg')
        effbkg["cmap"]    = self.ccube
        effbkg["outfile"] = self.effbkgfile
        effbkg["irfs"]    = irfs
        effbkg["expcube"] = self.Cubename
        effbkg["bexpmap"] = self.BinnedMapfile
        effbkg["efact"]   = efact
        self.run_retry_compress(effbkg)

    def GtAlphaBkg(self,epsilon=0.03):
        '''
        Hypercube - relative contributions to likelihood weights from different analysis components
        needs to be run once ALL components are done
        '''
        if (self.clobber=="no" and os.path.isfile(self.alphabkgfile)):
            #print("File exists and clobber is False")
            return(0)
        irfs,_ = utils.GetIRFS(self.Configuration['event']['evclass'],
                               self.Configuration['event']['evtype'])
        
        effbkg_files = glob.glob(\
                            self.folder+'/'+\
                            self.srcname+"*_"+\
                            self.modelname+\
                            "_effbkgfile.fits"+\
                            self.gzflag)
        effbkg_textfile = self.folder+'/'+self.srcname+"_"+self.modelname+"_effbkgfile.list"
        with open(effbkg_textfile, "w") as f:
            for effbkg_f in effbkg_files:
                f.write('{}\n'.format(effbkg_f))

        alphabkg = GtApp('gtalphabkg', 'AlphaBkg')
        alphabkg["inputs"]  = effbkg_textfile
        alphabkg["outfile"] = self.alphabkgfile
        alphabkg["epsilon"] = epsilon
        self.run_retry_compress(alphabkg)
    
    def GtWtsMap(self,epsilon=0.03):
        '''
        Cube of likelihood weight factors
        needs to be done once ALL components are done, component by component.
        '''
        if (self.clobber=="no" and os.path.isfile(self.alphabkgfile)):
            #print("File exists and clobber is False")
            return(0)
        irfs,_ = utils.GetIRFS(self.Configuration['event']['evclass'],
                               self.Configuration['event']['evtype'])
        wtsmap = GtApp('gtwtsmap', 'WtsMap')
        wtsmap["effbkgfile"] = self.effbkgfile
        wtsmap["alphafile"]  = self.alphabkgfile
        wtsmap["epsilon"]    = epsilon
        wtsmap["outfile"]     = self.wtsmapfile
        self.run_retry_compress(wtsmap)
        
