"""
gtfunction.py written by David Sanchez : david.sanchez@lapp.in2p3.fr
Collection of function to run the ST tools.
The description of each ST tool is given in the dedicated NASA website
The Observation class contains all the variables needed  to run the ST like file path, energy, time etc, ...
begun October 2010
"""
import os
from math import sqrt, log10
from gt_apps import evtbin, maketime, diffResps, expCube, expMap, srcMaps, model_map, filter
from GtApp import GtApp
from enrico import utils


class Observation:
    # init function of the Observation class. 
    # folder : folder where the produced fits files will be stored.
    # configuration is the confi of enrico (contains the variable)
    # convtyp is the convertion type of the Fermi events (0 : front, 1 : back, -1 : all) 
    def __init__(self,folder,Configuration,convtyp=-1,tag=""):
        
        #Read the configuration object and init all the variable
        self.Configuration = Configuration
        inttag = "_"+Configuration['file']['tag']
        if not(tag==""):
            inttag+="_"+tag

        self.srcname   = Configuration['target']['name']
        self.ft1       = Configuration['file']['event']
        self.ft2       = Configuration['file']['spacecraft']
        self.xmlfile   = Configuration['file']['xml']

        #Fits files
        self.eventfile = folder+'/'+self.srcname+inttag+"_Evt.fits"
        self.Cubename  = folder+'/'+self.srcname+inttag+"_ltCube.fits"
        self.Mapname   = folder+'/'+self.srcname+inttag+"_ExpMap.fits"
        self.BinnedMapfile = folder+'/'+self.srcname+inttag+"_BinnedMap.fits"
        self.cmapfile  = folder+'/'+self.srcname+inttag+"_CountMap.fits"
        self.lcfile    = folder+'/'+self.srcname+inttag+"_applc.fits"
        self.ccube     = folder+'/'+self.srcname+inttag+"_CCUBE.fits"
        self.scrMap    = folder+'/'+self.srcname+inttag+"_srcMap.fits"
        self.ModelMap  = folder+'/'+self.srcname+inttag+"_ModelMap.fits"
        self.BinDef  = folder+'/'+self.srcname+inttag+"_BinDef.fits"

        #Variables
        self.t1        = float(Configuration['time']['tmin'])
        self.t2        = float(Configuration['time']['tmax'])
        self.Emin      = float(Configuration['energy']['emin'])
        self.Emax      = float(Configuration['energy']['emax'])
        self.ra        = float(Configuration['space']['xref'])
        self.dec       = float(Configuration['space']['yref'])
        self.roi       = float(Configuration['space']['rad'])
        self.irfs      = Configuration['analysis']['irfs']

        #convtyp : the name of the irfs is updated
        if convtyp==0 :
            self.irfs += "::FRONT"
        if convtyp==1 :
            self.irfs += "::BACK"

        #Maps binning
        self.binsz     = self.Configuration['space']['binsz']
        self.npix      = int(2*self.roi/sqrt(2.)/self.binsz)
        self.convtyp   = convtyp

        #tool options
        self.clobber = self.Configuration['clobber']

    def printSum(self):
        """Print a summary of the value stored in the class"""
        print "Source\t=\t",self.srcname
        print "RA\t=\t",self.ra," degrees"
        print "Dec\t=\t",self.dec," degrees"
        print "Start\t=\t",self.t1,"  MET (s)"
        print "Stop\t=\t",self.t2,"  MET (s)"
        print "ROI\t=\t",self.roi," degrees"
        print "E min\t=\t",self.Emin," MeV"
        print "E max\t=\t",self.Emax," MeV"
        print "IRFs\t=\t",self.irfs

    def Gtbin(self):
        """Run gtbin with the CMAP option. A count map is produced"""
        evtbin['evfile'] = self.eventfile
        evtbin['scfile'] = self.ft2
        evtbin['outfile'] = self.cmapfile
        evtbin['algorithm'] = "CMAP"
        evtbin['nxpix'] = self.npix 
        evtbin['nypix'] = self.npix 
        evtbin['binsz'] = self.binsz
        evtbin['coordsys'] = self.Configuration['space']['coordsys']
        evtbin["emin"] = self.Emin
        evtbin["emax"] = self.Emax
        evtbin['xref'] = self.ra 
        evtbin['yref'] = self.dec
        evtbin['axisrot'] = 0
        evtbin['proj'] = self.Configuration['space']['proj'] 
        evtbin['clobber'] = self.clobber
        evtbin.run()

    def GtBinDef(self,filename):
        bindef = GtApp('gtbindef', 'Likelihood')
        bindef['bintype'] = 'T'
        bindef['binfile'] = filename
        bindef['outfile'] = self.BinDef
        bindef.run()

    def GtExposure(self):
        exposure = GtApp('gtexposure', 'Likelihood')
        exposure['infile'] = self.lcfile
        exposure['scfile'] = self.ft2
        exposure['irfs'] = self.irfs
        exposure['srcmdl'] = "none"
        exposure['specin'] = -self.Configuration['AppLC']['index']
        exposure['clobber'] = self.clobber
        exposure.run()

    def GtLCbin(self,dt=60):
        """Run gtbin with the LC option. the default dt is 60 sec and data can be rebinned after.
          Can also take a file as input to define the time bins"""
        evtbin['evfile'] = self.eventfile
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
        evtbin.run()

    def GtCcube(self):
        """Run gtbin with the CCUBE option"""
        Nbdecade = log10(self.Emax)-log10(self.Emin)#Compute the number of decade
        evtbin['evfile'] = self.eventfile
        evtbin['scfile'] = self.ft2
        evtbin['outfile'] = self.ccube
        evtbin['algorithm'] = "CCUBE"
        evtbin['nxpix'] = self.npix
        evtbin['nypix'] = self.npix
        evtbin['binsz'] = self.binsz
        evtbin['coordsys'] = self.Configuration['space']['coordsys']
        evtbin['xref'] = self.ra 
        evtbin['yref'] = self.dec
        evtbin["emin"] = self.Emin
        evtbin["emax"] = self.Emax
        evtbin['ebinalg'] = "LOG"
        evtbin['axisrot'] = 0
        evtbin['proj'] = self.Configuration['space']['proj'] #"AIT"
        #The number of bin is the number of decade * the number of bin per decade (given by the users)
        evtbin["enumbins"] = int(Nbdecade*self.Configuration['energy']['enumbins_per_decade'])
        evtbin['clobber'] = self.clobber
        evtbin.run()

    def GtBinnedMap(self):
        """Run the gtexpcube2 tool for binned analysis"""
        expcube2 = GtApp('gtexpcube2', 'Likelihood')
        expcube2['infile'] = self.Cubename
        expcube2['outfile'] = self.BinnedMapfile
        expcube2['cmap'] = self.ccube
        expcube2['irfs'] = self.irfs
        expcube2['emin'] = self.Emin
        expcube2['emax'] = self.Emax
        expcube2['coordsys'] = self.Configuration['space']['coordsys']
        expcube2['proj'] = self.Configuration['space']['proj'] #"AIT"
        expcube2['clobber'] = self.clobber
        expcube2.run()

    def FirstCut(self):
        """Run gtselect tool"""
        filter['infile'] = self.ft1
        filter['outfile'] = self.eventfile
        filter['ra'] = self.ra 
        filter['dec'] = self.dec 
        filter['rad'] = self.roi
        filter['emin'] = self.Emin
        filter['emax'] = self.Emax
        filter['tmin'] = self.t1
        filter['tmax'] = self.t2
        filter['zmax'] = self.Configuration['analysis']['zmax'] 
        filter['evclsmin'] = self.Configuration['analysis']['evclass']
        filter['evclsmax'] = 4
        filter['convtype'] = self.convtyp
        filter['clobber'] = self.clobber
        filter.run()

    def time_selection(self):
        """
        Do a GTI selection based on a file of time spans

        CFITSIO won't allow filenames (including filter expression) longer than
        ~1100 chars, so for selections that require very long filters (i.e.,
        more than ~30 time spans covered) we split the gtmktime calls into
        chunks of ~20 time spans.
        """
        eventlist = []
        last = False
        numbin = None
        while not last:
            selstr,numbin,last = utils.time_selection_string(self.Configuration,numbin)
            outfile = self.eventfile.replace('.fits','_{}'.format(numbin))
            self._RunMktime(selstr,outfile,'no')
            eventlist.append(outfile+'\n')
            maketime['outfile'] = outfile
            maketime.run()

        evlist_filename = self.eventfile.replace('.fits','.list')
        with open(evlist_filename,'w') as evlistfile:
            evlistfile.writelines(eventlist)

        # Redo first-cut to consolidate into single fits file (gtmktime does not accept lists!)
        ft1 = self.ft1 # Store FT1 to restore it later
        self.ft1 = evlist_filename
        self.FirstCut()
        self.ft1 = ft1

        # Clean cruft: all temp event files and event file list
        os.unlink(evlist_filename)
        for file in eventlist:
            os.unlink(file.strip()) # strip of endline char

    def MkTime(self):
        """compute GTI"""
        if self.Configuration['time']['file'] != '':
            self.time_selection()
        selstr = self.Configuration['analysis']['filter']
        outfile = self.eventfile+".tmp"
        self._RunMktime(selstr,outfile,'yes')
        os.system("mv "+self.eventfile+".tmp "+self.eventfile)

    def _RunMktime(self,selstr,outfile,roicut):
        """run gtmktime tool"""
        maketime['scfile'] = self.ft2
        maketime['filter'] = selstr#self.Configuration['analysis']['filter']
        maketime['roicut'] = roicut
        maketime['tstart'] = self.t1
        maketime['tstop'] = self.t2
        maketime['evfile']= self.eventfile
        maketime['outfile']=outfile
        maketime['clobber'] = self.clobber
        maketime.run()

    def DiffResps(self):
        """run gtdiffresp"""
        diffResps['evfile']=self.eventfile
        diffResps['scfile']=self.ft2
        diffResps['srcmdl']=self.xmlfile
        diffResps['irfs']=self.irfs
        diffResps['convert']="no"
        diffResps['evclsmin']=self.Configuration['analysis']['evclass']
        diffResps['clobber'] = self.clobber
        diffResps.run()
        print "\ndone"

    def ExpCube(self): 
        "Run gtltcube tool to produce livetime cube"
        expCube['evfile']=self.eventfile
        expCube['scfile']=self.ft2
        expCube['outfile']=self.Cubename
        expCube['dcostheta']=0.025
        expCube['binsz']=1
        expCube['phibins']=self.Configuration['space']['phibins']
        expCube['clobber'] = self.clobber
        expCube.run()

    def ExpMap(self):
        "Run gtexpmap for unbinned analysis"
        Nbdecade = log10(self.Emax)-log10(self.Emin)#Compute the number of decade
        expMap['evfile'] = self.eventfile
        expMap['scfile'] = self.ft2
        expMap['expcube'] = self.Cubename
        expMap['outfile'] = self.Mapname
        expMap['irfs'] = self.irfs
        expMap['srcrad'] = self.roi+10
        #The number of bin is the number of decade * the number of bin per decade (given by the users)
        expMap['nenergies'] = int(Nbdecade*self.Configuration['energy']['enumbins_per_decade'])
        expMap['clobber'] = self.clobber
        expMap.run() 

    def SrcMap(self):
        """Run gtsrcmap tool for binned analysis"""
        srcMaps['scfile'] = self.ft2
        srcMaps['expcube'] = self.Cubename
        srcMaps['cmap'] = self.ccube
        srcMaps['bexpmap'] = self.BinnedMapfile
        srcMaps['srcmdl']=self.xmlfile
        srcMaps['irfs']= self.irfs
        srcMaps['outfile'] = self.scrMap
        srcMaps['emapbnds']='no'
        srcMaps['clobber'] = self.clobber
        srcMaps.run()

    def ModelMaps(self,xml):
        """Run gtmodelmap tool for binned analysis and make a subtraction of the produced map
         with the count map to produce a residual map"""
        model_map['expcube'] = self.Cubename
        model_map['srcmaps'] = self.scrMap
        model_map['bexpmap'] = self.BinnedMapfile
        model_map['srcmdl'] = xml
        model_map["irfs"]=self.irfs
        model_map['outfile'] = self.ModelMap
        model_map['clobber'] = self.clobber
        model_map.run()
        #Compute the residual map
        utils.SubtractFits(self.cmapfile,self.ModelMap,self.Configuration)


