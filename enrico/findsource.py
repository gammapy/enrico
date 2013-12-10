#!/usr/bin/env python
import string,os
from enrico import utils
from enrico.fitmaker import FitMaker
from enrico.config import get_config
from enrico.gtfunction import Observation

def update_reg(config):
    lines=open(utils._dump_findsrcout(config),"r").readlines()[-4]
    word = string.split(lines)
    text = "\ncircle("+word[0]+","+word[1]+","+word[3]+")  # color=red text={Best Fit Position}"
    lines = open(utils._dump_reg(config),"r").readlines()
    lines.append(text)
    file_out = open(utils._dump_reg(config),"w")
    file_out.write(lines)
    file_out.close()

def FindSrc(infile):
    config = get_config(infile)
    folder = config['out']

    Obs = Observation(folder, config)
    utils._log('SUMMARY: ')
    Obs.printSum()

    FitRunner = FitMaker(Obs, config)

    if config["findsrc"]["FitsGeneration"]== "yes":
        config['analysis']['likelihood'] = 'unbinned'
        FitRunner.GenerateFits()

    FitRunner._log('gtfindsrc', 'Optimize source position')
    os.system("rm "+utils._dump_findsrcout(config))
    Obs.FindSource()
    update_reg(config)


if __name__ == '__main__':
    import sys
    try:
        infile = sys.argv[1]
    except:
        print('FATAL: Config file not found.')
        sys.exit(1)

    FindSrc(infile)

