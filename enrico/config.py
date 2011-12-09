"""Central place for config file handling"""
import sys
from os.path import join
from configobj import ConfigObj, flatten_errors
from validate import Validator
from environ import CONFIG_DIR

def get_config(infile, configspec=join(CONFIG_DIR, 'default.conf')):
    """Parse config file, and in addition:
    - include default options
    - exit with an error if a required option is missing"""
    config = ConfigObj(infile, configspec=configspec,
                       file_error=True)
    validator = Validator()
    # @todo: I'm not sure we always want to copy all default options here
    results = config.validate(validator, copy=True)
    if results != True:
        for (section_list, key, _) in flatten_errors(config, results):
            if key is not None:
                #label = 
                print('The "%s" key in the section "%s" failed validation' % 
                      (key, ', '.join(section_list)))
            else:
                print('The following section was missing:%s ' % 
                      ', '.join(section_list))
        print('   Please check your config file for missing and wrong options!')
        print('FATAL: Config file is not valid.')
        sys.exit(1)
    return config

# @todo: This doesn't work because missing values are invalid!!!
# Maybe fill those values by hand?
def get_default_config(configspec=join(CONFIG_DIR, 'default.conf')):
    return ConfigObj(None, configspec=configspec)

def query_config():
    """Make a new config object, asking the user for required options"""
    config = ConfigObj(indent_type='\t')
    print('Please provide the following required options:')
    config['out'] = raw_input('Output directory: ')
    config['target'] = {}
    config['target']['name'] = raw_input('Target Name : ')
    config['target']['ra'] = raw_input('Right Ascension: ')
    config['target']['dec'] = raw_input('Declination: ')
    config['target']['spectrum'] = raw_input('"Options are : PowerLaw, PowerLaw2, LogParabola, PLExpCutoff\nSpectral Model : ')
    config['space'] = {}
    config['space']['xref'] = config['target']['ra'] 
    config['space']['yref'] = config['target']['dec'] 
    config['space']['rad'] = raw_input('ROI Size: ')
    return get_config(config)
    
