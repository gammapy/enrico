import pyLikelihood
import pyIrfLoader

pyIrfLoader.Loader_go()
irf = pyIrfLoader.IrfsFactory.instance().create('P7SOURCE_V6::FRONT')
psf = irf.psf()
