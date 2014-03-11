#!/usr/bin/env python
import string
import numpy as np
from enrico.gtfunction import Observation
from enrico.utils import calcAngSepDeg
from enrico.config import get_config
import pyfits

def Runsrcprob(config):
    config['space']['rad'] = config['srcprob']['rad']
    Obs = Observation(config['out'], config, config['analysis']['convtype'], tag="srcprob")
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
      indices = energy.argsort()
      print "Energy\tAngular Sep\tProba"
      for i in xrange(min(10,indices.size)):
          angSep = calcAngSepDeg(config['target']["ra"],config['target']["dec"],ra[indices[i]],dec[indices[i]])
          print "%2.1f\t%2.3f\t\t%2.5f"%(energy[indices[i]],angSep,proba[indices[i]])



