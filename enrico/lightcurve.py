import os
from os.path import join
import numpy as np
import ROOT
import utils
import root_style
import plotting
import RunGTlike
import environ
from config import get_config
from submit import call


def PrepareLC(infile):
    """@todo: document me"""
    config = get_config(infile)

    Tag = config['file']['tag']
    Nbin = config['LightCurve']['NLCbin']
    tmin = config['time']['tmin']
    tmax = config['time']['tmax']

    config['UpperLimit']['TSlimit'] = config['LightCurve']['TSLightCurve']
    config['out'] += '/LightCurve'

    config['Spectrum']['ResultPlots'] = 'no'
    config['Ebin']['NumEnergyBins'] = 0
    config['UpperLimit']['envelope'] = 'no'

    AllConfigFile = []
    for i in xrange(Nbin):
        dt = (tmax - tmin) / Nbin
        config['time']['tmin'] = tmin + i * dt
        config['time']['tmax'] = tmin + (i + 1) * dt
        config['file']['tag'] = Tag + '_LC_' + str(i)
        filename = (config['out'] + "/Config_" +
                    str(config['time']['tmin']) + "_" +
                    str(config['time']['tmax']))

        if config['LightCurve']['Prepare'] == 'yes':
            config.write(open(filename, 'w'))

        AllConfigFile.append(filename)
    return AllConfigFile


def WriteToAscii(Time, TimeErr, Flux, FluxErr, TS, Npred, filename):
    """@todo: document me"""
    flc = open(filename, 'w')
    flc.write('Time (MET) Delta_Time Flux(ph cm-2 s-1) '
              'Delta_Flux TS Npred\n')
    for i in xrange(len(Time)):
        flc.write(str(Time[i]) + "\t" + str(TimeErr[i]) + "\t" +
                  str(Flux[i]) + "\t" + str(FluxErr[i]) + "\t" +
                  str(TS[i]) + "\t" + str(Npred[i]) + "\n")
    flc.close()


def MakeLC(infile):
    """@todo: document me"""
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    root_style.RootStyle()

    Doplot = True

    enricodir = environ.DIRS.get('ENRICO_DIR')

    config = get_config(infile)

    folder = config['out']
    os.system('mkdir -p ' + folder + '/LightCurve')

    AllConfigFile = PrepareLC(infile)

    Nbin = config['LightCurve']['NLCbin']

    Time = []
    TimeErr = []
    Flux = []
    FluxErr = []
    Npred = []
    TS = []

    # @todo: Do you mean '=' instead of '==' here?
    # config['Spectrum']['FitsGeneration'] == config['LightCurve']['FitsGeneration']
    if (config['LightCurve']['FitsGeneration'] == 'yes' or
        config['LightCurve']['Re-Fit'] == 'yes'):

        for i in xrange(Nbin):
            if config['LightCurve']['Submit'] == 'yes':
                Doplot = False
                cmd = "enrico_fit " + AllConfigFile[i]
                scriptname = folder + "/LightCurve/LC_Script_" + str(i) + ".sh"
                JobLog = folder + "/LightCurve/LC_Job_" + str(i) + ".log"
                JobName = config['target']['name'] + "_" + str(i) + ".log"

                call(cmd, enricodir, scriptname, JobLog, JobName)
            else:
                RunGTlike.run(AllConfigFile[i])

    if config['LightCurve']['Plot'] == 'yes' and Doplot:
        print "Reading output files"
        for i in xrange(Nbin):
            CurConfig = get_config(AllConfigFile[i])
            try:
                result = utils.ReadResult(CurConfig)
            except:
                print('WARNING : fail reading the configuration file : ',
                      AllConfigFile[i])
                print "Job Number : ", i
                print "Please have a look at this job log file"
                continue

            Time.append((result.get("tmax") + result.get("tmin")) / 2.)
            TimeErr.append((result.get("tmax") - result.get("tmin")) / 2.)
            if 'Ulvalue' in result:
                Flux.append(result.get("Ulvalue"))
                FluxErr.append(0)
            else:
                Flux.append(result.get("Flux"))
                FluxErr.append(result.get("dFlux"))

            Npred.append(result.get("Npred"))
            TS.append(result.get("TS"))

        TS = np.array(TS)
        Npred = np.array(Npred)
        Time = np.array(Time)
        TimeErr = np.array(TimeErr)
        Flux = np.array(Flux)
        FluxErr = np.array(FluxErr)

        filebase = join(folder, 'LightCurve', config['target']['name'])

        if config['LightCurve']['DiagnosticPlots'] == 'yes':
            gTHNpred, TgrNpred = plotting.PlotNpred(Npred, Flux, FluxErr)
            CanvNpred = ROOT.TCanvas()
            gTHNpred.Draw()
            TgrNpred.Draw('zP')
            CanvNpred.Print(filebase + '_Npred.eps')
            CanvNpred.Print(filebase + '_Npred.C')

            gTHTS, TgrTS = plotting.PlotTS(Time, TimeErr, TS)
            CanvTS = ROOT.TCanvas()
            gTHTS.Draw()
            TgrTS.Draw('zP')
            CanvTS.Print(filebase + '_TS.eps')
            CanvTS.Print(filebase + '_TS.C')

        gTHLC, TgrLC, ArrowLC = plotting.PlotLC(Time, TimeErr, Flux, FluxErr)
        CanvLC = ROOT.TCanvas()
        gTHLC.Draw()
        TgrLC.Draw('zP')

        for i in xrange(len(ArrowLC)):
            ArrowLC[i].Draw()

        CanvLC.Print(filebase + '_LC.eps')
        CanvLC.Print(filebase + '_LC.C')

        filename = filebase + "_results.dat"
        print "Write to Ascii file : ", filename
        WriteToAscii(Time, TimeErr, Flux, FluxErr, TS, Npred, filename)
