"""
simulationmaker.py written by David Sanchez : david.sanchez@lapp.in2p3.fr
begun Jan 2022
"""
import numpy as np
import string
from enrico import utils
from enrico import Loggin
from enrico import environ
from enrico.gtfunction import Observation

class SimulationMaker(Loggin.Message):

    """Collection of functions to run the simulation"""
    def __init__(self, folder, config, tag=""):
        super(SimulationMaker,self).__init__()
        Loggin.Message.__init__(self)
        self.obs = Observation(folder, config, tag=tag)


    def run(self):
        self.info("Run SimulationMaker")
        self.obs.Obssim()