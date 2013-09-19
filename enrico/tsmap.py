#!/usr/bin/env python
import os,sys,logging
import pyfits
import SummedLikelihood
import utils
from submit import call
import environ
from RunGTlike import Analysis,GenAnalysisObjects
from gtfunction import Observation

class TSMap:
    # This class groups all the needed functions and 
    # variables to compute a TS map
    # Variables : RAref, DECref and binsz are related to the cmap used
    # to define the grip for the TS map and are respectively 
    # the center of the map and the size of 1 bin in degrees
    # tsfolder is the location where the produced files file be stored
    # npix is the number of pixel of the TS map
    # infile is the used configuration file

    # A TS map can be computed by submiting jobs to a cluster
    # There is two way to compute a TS map : Each job
    # computes the TS in one pixel or in one row
    def __init__(self,config,infile):
        self.config = config
        self.config['Spectrum']['FitsGeneration'] = 'no'
        self.tsfolder = self.config['out']+"/TSMap"
        self.TSfits = self.config['target']['name']+'_'+self.config['file']['tag']+"_TSMap.fits"
        self.infile = infile
        self.npix = self.config['TSMap']['npix']
        # Read the cmap produced before to get the grid for the TS map
        FitRunner = Observation(self.config['out'], self.config)
        try :
             cmap = pyfits.open(FitRunner.cmapfile)
        except :
             logging.error('Count map not found.')
             sys.exit(1)
        
        npix_im = min(cmap[0].header['NAXIS1'],cmap[0].header['NAXIS2'])
        self.npix = min(self.npix,npix_im)
        self.RAref = cmap[0].header['CRVAL1']
        self.DECref = cmap[0].header['CRVAL2']
        self.binsz = cmap[0].header['CDELT1']

    def _launch(self,ra,dec,i,j):
        """ Launch a job (either pixel evaluation or row evaluation). 
        Can be the submittion of a job to a cluster """
        enricodir = environ.DIRS.get('ENRICO_DIR')
        fermidir = environ.DIRS.get('FERMI_DIR')
        cmd = enricodir+"/enrico/tsmap.py "+os.getcwd()+"/"+self.infile +" "+ str(ra) +" "+ str(dec) +" "+ str(i) +" "+ str(j) #cmd line to send

        if self.config['Submit'] == 'yes':
            prefix = self.tsfolder + "/TSMap_" + str(i) +"_"+ str(j)
            scriptname = prefix + "_Script.sh"
            JobLog = prefix + "_Job.log"
            JobName = (self.config['target']['name'] + "_TSMap_" + str(i) +"_"+ str(j))
            call(cmd, enricodir, fermidir, scriptname, JobLog, JobName) #submition
        else : 
            os.system(cmd) #run directly 

    def _PixelFile(self,i,j):
        """ return the name of a file where the result of 1 pixel 
        evaluation will be stored"""
        return self.tsfolder+'/Pixel_'+str(i)+'_'+str(j)

    def FitOnePixel(self,ra,dec,i,j) :
        """Run a evaluation of the pixel (i,j) corresponding to position (ra,dec)"""
        outXml = utils._dump_xml(self.config)
        folder = self.config['out']
        _,Fit = GenAnalysisObjects(self.config,xmlfile=outXml) #get the Fit object

        src = GetSrc(Fit,ra,dec) # get the Source object at postion ra dec

        if self.config['TSMap']['RemoveTarget'] : # remove the target is asked
            Fit.deleteSource(self.config['target']['name'])

        if self.config['TSMap']['Re-Fit']: # reoptimze before is asked
            Fit.fit(0,optimizer=self.config['fitting']['optimizer'])

        for par in xrange(Fit.logLike.getNumParams()): # freeze all the source parameters
            Fit[par].setFree(0)

        Fit.addSource(src)# add a spurious source
        # dump the *new* xml file
        Fit.writeXml(self.tsfolder+"/model_"+str(ra)+"_"+str(dec)+".xml") 

        # a new Fit object with the new xml file is needed.
        # just changing the position of the spurious source does not work 
        _,TSFit = GenAnalysisObjects(self.config,xmlfile=self.tsfolder+"/model_"+str(ra)+"_"+str(dec)+".xml") #get the Fit object

        TSFit.fit(0,optimizer=self.config['fitting']['optimizer'])

        # save the result
        fsave = open(self._PixelFile(i,j),'w')
        fsave.write(str(ra)+"\t"+str(dec)+"\t"+
                    str(TSFit.Ts("Spurious"))+"\t"+
                    str(TSFit.logLike.value()))
        fsave.close()

    def FitOneRow(self,ra,i) :
        """ function which run the evaluation of 1 row of the TS map
        using a loop and calling the fit for 1 pixel"""
        for j in xrange(self.npix):
            dec = self.DECref + self.binsz*(j-self.npix/2.)
            print 'FitOneRow ',dec
            self.FitOnePixel(ra,dec,i,j)

    def runTSMap(self,row=-1,column=-1) :
        """ Run a TS map using the configuration file given"""
        folder = self.config['out']
        os.system('mkdir -p ' + self.tsfolder)

        # This part is used to rerun either a row or a pixel.
        if row>0:#rerun only 1 row
            ra = self.RAref + self.binsz*(row-self.npix/2.)
            if column>0: #rerun only 1 pixel
                dec = self.DECref + self.binsz*(column-self.npix/2.)
                print 'Run Pixel evaluation at ',ra,' ',dec
                self._launch(ra,dec,row,column) 
            else :
                print 'Run Row evaluation at ',ra
                self._launch(ra,0,row,0)
            return 

        # Normal operation : all row and piwel are computed
        for i in xrange(self.npix): #loop over the X axis
            ra = self.RAref + self.binsz*(i-self.npix/2.)
            if self.config['TSMap']['method'] == 'row' : # a row is evaluated in one job
#                if row<0 or i==row:
                 print 'Run Row evaluation at ',ra
                 self._launch(ra,0,i,0)
            else : # each pixel is evaluated by one job
                for j in xrange(self.npix): #loop over the Y axis
#                    if (row<0 and column<0) or (i==row and column<0) or (i==row and j==column):
                     dec = self.DECref + self.binsz*(j-self.npix/2.)
                     print 'Run Pixel evaluation at ',ra,' ',dec
                     self._launch(ra,dec,i,j) 

    def PlotTSmap(self) :
        """ Gather the results of the evaluation of 
        each pixel and fill a fits file"""
        folder = self.config['out']

        # Read the cmap produced before to get the grid for the TS map
        FitRunner = Observation(folder, self.config)
        try :
             header = pyfits.getheader(FitRunner.cmapfile)
        except :
             logging.error('Count map not found.')
             sys.exit(1)
        data = pyfits.getdata(FitRunner.cmapfile)*0.
        npix_im = min(header['NAXIS1'],header['NAXIS2'])
        npix = min(self.config['TSMap']['npix'],npix_im)
        Xref = header['CRPIX1']
        Yref = header['CRPIX2']
        binsz = header['CDELT1']

        import string # read the results
        for i in xrange(npix):
            for j in xrange(npix):
                try : 
                    lines = open(self._PixelFile(i,j),"r").readlines()
                    Value = float(string.split(lines[0])[2])
                except :
                    print "Cannot find, open or read ",self._PixelFile(i,j)
                    Value = 0.
                data[Xref+ (i-npix/2.)][Yref+ (j-npix/2.)] = Value

        # save in a fits files
        pyfits.writeto(folder+"/"+self.TSfits,data,header)
        print "TS Map saved in "+folder+"/"+self.TSfits


def GetSrc(Fit,ra,dec):
    """ return a Source object by cloning a pointlike source from the Fit object
    the source is rename, move to (ra,dec) and the spetral model is change to PowerLaw"""
    ind = 0
    for name in  Fit.logLike.srcNames() :
        src = Fit.logLike.getSource(name).clone()
        if str(type(src))=="<type 'NoneType'>" :
            continue
        if src.getType() == 'Point' :
            src.setName("Spurious")
            src.setSpectrum("PowerLaw")
            src.getSrcFuncs()['Position'].getParam('RA').setValue(ra)
            src.getSrcFuncs()['Position'].getParam('DEC').setValue(dec)

            src.getSrcFuncs()['Spectrum'].getParam('Prefactor').setBounds(1e-5,1e5)
            src.getSrcFuncs()['Spectrum'].getParam('Prefactor').setScale(1e-9)
            src.getSrcFuncs()['Spectrum'].getParam('Index').setBounds(-5,0)

            src.getSrcFuncs()['Spectrum'].getParam('Scale').setValue(300)
            src.getSrcFuncs()['Spectrum'].getParam('Scale').setBounds(1e-5,1e5)

            return src



if __name__ == '__main__':
    try:
        infile = sys.argv[1]
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)
	

    from enrico.config import get_config
    config = get_config(infile)
    TSmap = TSMap(config,infile)

    if len(sys.argv)== 6 :
        if TSmap.config['TSMap']['method'] == 'row' :
            TSmap.FitOneRow(float(sys.argv[2]),int(sys.argv[4]))
        else :
            TSmap.FitOnePixel(float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]))
    else :
        print "Wrong number of arguments"
        sys.exit(1)
