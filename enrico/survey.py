"""Utilities for making survey maps"""
import os.path
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

from environ import *

class MapMaker(object):
    """Make survey maps using binned Fermi tools"""
    # Available preprocessing steps
    STEPS = ['count_map']

    def __init__(self, input_dir, output_dir,
                 chatter=4, clobber='no', debug='no'):
        self.input_dir = input_dir
        try:
            os.makedirs(output_dir)
            log.info('MKDIR: %s' % output_dir)
        except OSError:
            log.info('%s exists alreay' % output_dir)
        os.chdir(output_dir)
        self.chatter = chatter
        self.clobber = clobber
        self.debug = debug
        self.nxpix = 3600
        self.nypix = 1800
        self.binsz = 0.1
        self.coordsys = 'GAL'
        self.xref = 0
        self.yref = 0
        self.proj = 'AIT'

    def process(self, steps=None):
        """Execute steps"""
        log.info('input_dir: %s' % self.input_dir)
        if steps == None:
            steps = self.STEPS
        log.info('steps: %s' % steps)
        if 'count_map' in steps:
            self._count_map()

    def _count_map(self):
        """Run gtbin"""
        from gt_apps import evtbin as tool
        self._set_common_tool_options(tool)
        tool['evfile'] = os.path.join(self.input_dir, 'gtmktime.fits')
        tool['scfile'] = 'NONE'
        tool['outfile'] = 'map.fits'
        tool['algorithm'] = 'CMAP'
        tool['nxpix'] = self.nxpix
        tool['nypix'] = self.nypix
        tool['binsz'] = self.binsz
        tool['coordsys'] = self.coordsys
        tool['xref'] = self.xref
        tool['yref'] = self.yref
        tool['axisrot'] = 0
        tool['proj'] = self.proj
        tool.run()

    def _set_common_tool_options(self, tool):
        tool['chatter'] = self.chatter
        tool['clobber'] = self.clobber
        tool['debug'] = self.debug
        tool['gui'] = 'no'
        tool['mode'] = 'ql'
