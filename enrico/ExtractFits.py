#!/usr/bin/env python
from RunGTlike import Analysis
import sys
from enrico import utils
from enrico import Loggin
from enrico.config import get_config

if __name__ == '__main__':

    mes = Loggin.Message()
    try:
        infile = sys.argv[1]
    except:
        print('Usage: '+sys.argv[0]+' <config file name>')
        mes.error('Config file not found.')
    """Run an  Fermi analysis to generate FITS files by reading a config file"""
    config = get_config(infile)
    folder = config['out']
    utils.mkdir_p(folder)
    Analyse = Analysis(folder, config, \
                configgeneric=config,\
                tag="", verbose = 1)
