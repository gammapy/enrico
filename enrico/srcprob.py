#!/usr/bin/env python
import string
import numpy as np
from enrico.gtfunction import Observation
from enrico.utils import calcAngSepDeg
from enrico.config import get_config
import pyfits
from enrico import Loggin

def Print(indices,config,ra,dec,proba,energy,time):
      mes = Loggin.Message()
      mes.info("Energy\tAngular Sep\tProba\tTime")
      for i in xrange(min(10,indices.size)):
          angSep = calcAngSepDeg(config['target']["ra"],config['target']["dec"],ra[indices[indices.size-1-i]],dec[indices[indices.size-1-i]])
          print "%2.3e\t%2.3f\t\t%2.5f\t%2.1f"%(energy[indices[indices.size-1-i]],angSep,proba[indices[indices.size-1-i]],time[indices[indices.size-1-i]])


def Runsrcprob(config):
    config['space']['rad'] = config['srcprob']['rad']
    Obs = Observation(config['out'], config, config['analysis']['convtype'], tag="srcprob")
    if config['srcprob']['FitsGeneration'] =='yes':
        Obs.FirstCut()
        if config['analysis']['ComputeDiffrsp'] == 'yes':
            Obs.DiffResps()
    Obs.SrcProb()
    probfile=pyfits.open(Obs.Probfile)
    srclist = open(config['srcprob']['srclist'],"r").readlines()
    for src in srclist:
      proba = probfile[1].data.field(string.split(src)[0])
      energy = probfile[1].data.field("ENERGY")
      ra = probfile[1].data.field("RA")
      dec = probfile[1].data.field("DEC")
      time = probfile[1].data.field("TIME")
      indices = energy.argsort()
      mes = Loggin.Message()
      mes.info( "Results sorted by decreasing energy")
      Print(indices,config,ra,dec,proba,energy,time)
      print 
      mes.info( "Results sorted by decreasing probability")
      indices = proba.argsort()
      Print(indices,config,ra,dec,proba,energy,time)


