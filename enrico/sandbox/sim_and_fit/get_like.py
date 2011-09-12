"""Quick access to the like object for interactive testing:
$ python -i get_like.py
>>> <do with obs and like what you want>"""
from fermi import SCFILE, DATA, IRFS, OPTIMIZER
from UnbinnedAnalysis import UnbinnedObs, UnbinnedAnalysis

print 'Creating UnbinnedObs'
obs = UnbinnedObs(eventFile='gtmktime.fits',
                  scFile=SCFILE,
                  expMap='gtexpmap.fits',
                  expCube=DATA + '/11month_ltcube.fits',
                  irfs=IRFS)

print 'Creating UnbinnedAnalysis'
like = UnbinnedAnalysis(obs,
                        srcModel='crab_1FGL.xml',
                        optimizer=OPTIMIZER)

print '-'*60
print obs
print '-'*60
print like
print '-'*60