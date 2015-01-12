"""Central place for config file handling"""
import sys
from os.path import join
from enrico.extern.configobj import ConfigObj, flatten_errors
from enrico.extern.validate import Validator
from enrico.environ import CONFIG_DIR, DOWNLOAD_DIR
from enrico import Loggin


def get_config(infile, configspec=join(CONFIG_DIR, 'default.conf')):
    """Parse config file, and in addition:
    - include default options
    - exit with an error if a required option is missing"""
    config = ConfigObj(infile, configspec=configspec,
                       file_error=True)
    validator = Validator()
    # @todo: I'm not sure we always want to copy all default options here
    results = config.validate(validator, copy=True)
    mes = Loggin.Message()
    if results != True:
        for (section_list, key, _) in flatten_errors(config, results):
            if key is not None:
                mes.warning('The "%s" key in the section "%s" failed validation' %
                      (key, ', '.join(section_list)))
            else:
                mes.warning('The following section was missing:%s ' %
                      ', '.join(section_list))
        mes.warning('   Please check your config file for missing '
              'and wrong options!')
        mes.error('Config file is not valid.')
    return config


# @todo: This doesn't work because missing values are invalid!!!
# Maybe fill those values by hand?
def get_default_config(configspec=join(CONFIG_DIR, 'default.conf')):
    return ConfigObj(None, configspec=configspec)


def query_config():
    import os
    """Make a new config object, asking the user for required options"""
    config = ConfigObj(indent_type='\t')
    mes = Loggin.Message()
    mes.info('Please provide the following required options [default] :')
    config['out'] = os.getcwd()
    out = raw_input('Output directory ['+config['out']+'] : ')
    if not(out=='') :
      config['out'] = out

#    Informations about the source
    config['target'] = {}
    config['target']['name'] = raw_input('Target Name : ')
    config['target']['ra'] = raw_input('Right Ascension: ')
    config['target']['dec'] = raw_input('Declination: ')

    config['target']['redshift'] = '0'
    redshift = raw_input('redshift, no effect if null [0] : ')
    if not(redshift=='') :
        config['target']['redshift'] = redshift
        config['target']['ebl_model'] = raw_input('ebl model to used\n'
                                                   '0=Kneiske, 1=Primack05, 2=Kneiske_HighUV, 3=Stecker05, ' 
                                                   '4=Franceschini, 5=Finke, 6=Gilmore : ')
        
    message = ('Options are : PowerLaw, PowerLaw2, LogParabola, '
               'PLExpCutoff, Generic\nGeneric is design to allow the user to fit with non-supported models\n'
                'EBL absorption can be added for PowerLaw2, LogParabola, PLExpCutoff\n'
                'Spectral Model [PowerLaw] : ')
    config['target']['spectrum'] = 'PowerLaw'
    model = raw_input(message)
    if not(model=='') :
      config['target']['spectrum'] = model

#    informations about the ROI
    config['space'] = {}
    config['space']['xref'] = config['target']['ra']
    config['space']['yref'] = config['target']['dec']
    config['space']['rad'] = '15'
    roi = raw_input('ROI Size [15] : ')
    if not(roi=='') :
      config['space']['rad'] = roi

#    informations about the input files
    config['file'] = {}
    config['file']['spacecraft'] = DOWNLOAD_DIR+'/lat_spacecraft_merged.fits'
    ft2 = raw_input('FT2 file ['+config['file']['spacecraft']+'] : ')
    if not(ft2=='') :
      config['file']['spacecraft'] = ft2
    config['file']['event'] = DOWNLOAD_DIR+'/events.lis'
    ft1list = raw_input('FT1 list of files ['+config['file']['event']+'] : ')
    if not(ft1list=='') :
      config['file']['event'] = ft1list
    config['file']['xml'] = config['out']+'/'+config['target']['name']+'_'+config['target']['spectrum']+'_model.xml'
    tag = raw_input('tag [LAT_Analysis] : ')
    if not(tag=='') :
      config['file']['tag'] = tag
    else :
      config['file']['tag'] = 'LAT_Analysis'

#    informations about the time
    config['time'] = {}
    tmin = raw_input('Start time [239557418] : ')
    if not(tmin=='') :
      config['time']['tmin'] = tmin
    else :
      config['time']['tmin'] = '239557418'
    tmax = raw_input('End time [334165418] : ')
    if not(tmax=='') :
      config['time']['tmax'] = tmax
    else :
      config['time']['tmax'] = '334165418'

#    informations about the energy
    config['energy'] = {}
    emin = raw_input('Emin [100] : ')
    if not(emin=='') :
      config['energy']['emin'] = emin
    else :
      config['energy']['emin'] = '100'
    emax = raw_input('Emax [300000] : ')
    if not(emax=='') :
      config['energy']['emax'] = emax
    else :
      config['energy']['emax'] = '300000'

    return get_config(config)
