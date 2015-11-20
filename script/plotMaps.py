import os,sys
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


try:
  from aplpy import FITSFigure
except:
  print 
  sys.exit("Module aplpy missing")

#enrico imports
try:
  from enrico.config import get_config
  from enrico.gtfunction import Observation
except:
  sys.exit("Enrico Module missing? Fermi ST installed?")

try :
  vmax = float(sys.argv[2])
  print "use user vmax for the color scale: "+sys.argv[2]
except:
  print "use default vmax for the color scale: 8"
  print "you can change it by typing python "+sys.argv[0]+" config [vmax]"
  vmax=8

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
obs = Observation(config["out"], config)

counts = dict(label='Counts', filename=obs.cmapfile)
model = dict(label='Model', filename=obs.ModelMap)
residuals = dict(label='Residuals', filename=config['out'] + "/" + config['target']['name']+'_Residual_Model_cmap.fits')
images = [counts, model, residuals][::-1]


# Determine image center and width / height
dpi = 2000
header = fits.getheader(images[0]['filename'])
wcs = WCS(header)
header['NAXIS1'] / dpi
header['NAXIS2'] / dpi
lon, lat = header['NAXIS1'] / 2., header['NAXIS2'] / 2.
x_center, y_center = wcs.wcs_pix2world(lon, lat, 0)
radius = header['CDELT2'] * header['NAXIS2'] / 2.

# Computing the sub-figure sizes is surprisingly hard
figsize=(5, 15)
figure = mpl.figure(figsize=figsize)
axis_ratio = figsize[0] / float(figsize[1])
edge_margin_x = 0.12
edge_margin_y = edge_margin_x * axis_ratio
edge_margin_x_up = 0.01
edge_margin_y_up = edge_margin_x_up * axis_ratio
inner_margin_x = 0.1
inner_margin_y = inner_margin_x * axis_ratio
size_x = (1 - edge_margin_x - edge_margin_x_up)
size_y = (1 - edge_margin_y - edge_margin_y_up - 2 * inner_margin_y) / 3

for ii, image in enumerate(images):
    subplot = [edge_margin_x, edge_margin_y + ii * (size_y + inner_margin_y), size_x, size_y]
    f = FITSFigure(image['filename'], figure=figure, subplot=subplot)
    f.recenter(x_center, y_center, 0.95 * radius)
    set_hgps_style(f)
    f.show_colorscale(vmin=-1, vmax=vmax, stretch='power', exponent=1, cmap='jet') #vmid=-3, stretch='log', )
    # TODO: overplot sources  
#    f.show_regions("sources.reg")

filename = config["out"]+'/Maps.eps'
print('Writing {}'.format(filename))
figure.savefig(filename)
filename = config["out"]+'/Maps.png'
print('Writing {}'.format(filename))
figure.savefig(filename)
