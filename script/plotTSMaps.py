
try:
  import matplotlib
except:
  sys.exit("Module astropy missing")

matplotlib.use('Agg')
import matplotlib.pyplot as mpl

try:
  from astropy.io import fits
except:
  sys.exit("Module astropy missing")

from astropy.wcs import WCS
import os,sys

try:
  from aplpy import FITSFigure
except:
  print 
  sys.exit("Module aplpy missing")

#enrico imports
try:
  from enrico.config import get_config
except:
  sys.exit("Enrico Module missing?")

try :
  vmax = float(sys.argv[2])
  print "use user's vmax for the color scale: "+sys.argv[2]
except:
  print "use default vmax for the color scale: 40"
  print "you can change it by typing python "+sys.argv[0]+" config [vmax] [stretch]\n stretch can be linear (default), log, sqrt, arcsinh or power"
  vmax=40

try :
  st = [sys.argv[3]]
except:
  st = ['linear']

stretch =[ 'linear', 'log', 'sqrt', 'arcsinh', 'power']
if not(bool(sum(map(lambda x: x in st, stretch)))):
  print "use user's stretch for the color scale: "+st[0]
else:
  print "use default stretch value: "+st[0]


def set_hgps_style(f):
    """Set HGPS style for a f = aplpy.FITSFigure"""
    #f.set_tick_labels_font(size='small')
    #f.set_axis_labels_font(size='small')
    f.ticks.set_xspacing(2)
    f.ticks.set_yspacing(2)
    f.ticks.set_linewidth(1.5)
    f.tick_labels.set_xformat('dd')
    f.tick_labels.set_yformat('dd')
    f.tick_labels.set_style('colons')
    #f.tick_labels.set_font(size='small')
    f.axis_labels.set_xtext('Right Ascension (deg)')
    f.axis_labels.set_ytext('Declination (deg)')

try:
  infile = sys.argv[1]
except:
  "Please provide a configuration file"
  sys.exit(1)

print 
config = get_config(infile)
tsfile = config['out']+'/'+config['target']['name']+'_'+config['file']['tag']+"_TSMap.fits"

tsimage = dict(label='TSMap', filename=tsfile)

# Determine image center and width / height
dpi = 2000
header = fits.getheader(tsimage['filename'])
wcs = WCS(header)
header['NAXIS1'] / dpi
header['NAXIS2'] / dpi
lon, lat = header['NAXIS1'] / 2., header['NAXIS2'] / 2.
x_center, y_center = wcs.wcs_pix2world(lon, lat, 0)
radius = header['CDELT2'] * header['NAXIS2'] / 2.

# Computing the sub-figure sizes is surprisingly hard
figsize=(5, 5)
figure = mpl.figure(figsize=figsize)

f = FITSFigure(tsimage['filename'], figure=figure)
f.recenter(x_center, y_center, 0.95 * radius)
set_hgps_style(f)
f.show_colorscale(vmin=1e-5, vmax=vmax, stretch=st[0], exponent=1, cmap='jet') #vmid=-3, stretch='log', )

f.show_colorbar()

filename = config["out"]+'/TSMaps.eps'
print('Writing {}'.format(filename))
figure.savefig(filename)
filename = config["out"]+'/TSMaps.png'
print('Writing {}'.format(filename))
figure.savefig(filename)
