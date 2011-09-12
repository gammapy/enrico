# Simulate some data according to sim_model.xml
# Fit it according to fit_model.xml
# Then you should compare sim_model.xml to fit_result.xml

# Files
SIM_MODEL = sim_model.xml
SIM_LIST = sim_model.lis
SIM_EVROOT = mc
SIM_EVFILE = $(SIM_EVROOT)_events_0000.fits
SCFILE = $(SIM_EVROOT)_scData_0000.fits
EVFILE = events.fits
LTCUBE = gtltcube.fits
EXPMAP = gtexpmap.fits
FIT_MODEL = fit_model.xml
RESULT = result

# Parameters
SIMTIME = 7.88e4 # 7.88e7
SEED = 0

RA = 0
DEC = 0
RAD = 5 # 15 ROI fov (deg)
SRCRAD = 7 # 20 expmap fov (deg)
SRCPIX = 120 # 120 expmap fov (pix)
BINSZ = 0.1
NXPIX = 150
NYPIX = 150

EMIN = 1e4
EMAX = 1e6
NENERGIES = 5
ZMAX = 105
EVCLSMIN = 0
EVCLSMAX = 4
IRFS = P6_V3_DATACLEAN
FILTER = 'DATA_QUAL==0' # bogus to avoid error when running gtmktime

TIME = time
CHATTER = 4
CLOBBER = yes
DEBUG = yes


all: $(SIM_EVFILE) $(RESULT).xml 

$(SIM_EVFILE): $(SIM_MODEL) $(SIM_LIST)
	$(TIME) gtobssim infile=$(SIM_MODEL) srclist=$(SIM_LIST) scfile=none \
	evroot=$(SIM_EVROOT) simtime=$(SIMTIME) ltfrac=1 tstart=239557417 \
	nevents=no maxtime=3e8 startdate=2001-01-01 \
	use_ac=yes ra=$(RA) dec=$(DEC) radius=$(RAD) \
	emin=$(EMIN) emax=$(EMAX) irfs=$(IRFS) seed=$(SEED) \
	chatter=$(CHATTER) clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtobssim.log 2>&1

# Simulated events have class 0. We don't have to apply a cut.
gtselect.fits:
	$(TIME) gtselect infile=$(SIM_EVFILE) outfile=gtselect.fits \
	ra=$(RA) dec=$(DEC) rad=$(RAD) tmin=INDEF tmax=INDEF \
	emin=$(EMIN) emax=$(EMAX) \
	zmax=$(ZMAX) evclsmin=$(EVCLSMIN) evclsmax=$(EVCLSMAX) \
	chatter=$(CHATTER) clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtselect.log 2>&1

$(EVFILE): gtselect.fits
	$(TIME) gtmktime scfile=$(SCFILE) filter=$(FILTER) roicut=yes \
	evfile=gtselect.fits outfile=$(EVFILE) \
	chatter=$(CHATTER) clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtmktime.log 2>&1

gtbin.fits: $(EVFILE)
	time gtbin evfile=$(EVFILE) scfile=none outfile=gtbin.fits algorithm=CMAP \
	nxpix=$(NXPIX) nypix=$(NYPIX) binsz=$(BINSZ) c
	oordsys=CEL xref=$(RA) yref=$(DEC) axisrot=0 proj=CAR \
	chatter=$(CHATTER) clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtbin.log 2>&1

$(LTCUBE): $(EVFILE) $(SCFILE)
	time gtltcube evfile=$(EVFILE) scfile=$(SCFILE) outfile=$(LTCUBE) \
	dcostheta=0.025 binsz=1 zmax=$(ZMAX) \
	chatter=$(CHATTER) clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtltcube.log 2>&1

$(EXPMAP): $(EVFILE) $(SCFILE) $(LTCUBE)
	time gtexpmap evfile=$(EVFILE) scfile=$(SCFILE) expcube=$(LTCUBE) outfile=$(EXPMAP) \
	irfs = $(IRFS) srcrad=$(SRCRAD) nlong=$(SRCPIX) nlat=$(SRCPIX) nenergies=$(NENERGIES) \
	chatter=1 clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtexpmap.log 2>&1

$(RESULT).xml: $(LTCUBE) $(FIT_MODEL) $(EVFILE)
	time gtlike irfs=$(IRFS) expcube=$(LTCUBE) srcmdl=$(FIT_MODEL) sfile=$(FIT_RESULT) \
	results=$(RESULT).dat specfile=$(RESULT).fits statistic=UNBINNED optimizer=MINUIT \
	evfile=$(EVFILE) scfile=$(SCFILE) expmap=$(EXPMAP) \
	cmap=none bexpmap=none \
	chatter=$(CHATTER) clobber=$(CLOBBER) debug=$(DEBUG) \
	| tee gtlike.log 2>&1

gttsmap.fits: gtexpmap.fits 
	time gttsmap evfile=gtmktime.fits scfile=$(SCFILE) \
	expmap=gtexpmap.fits expcube=$(LTCUBE) \
	srcmdl=background.xml outfile=gttsmap.fits irfs=$(IRFS) optimizer=MINUIT \
	ftol=1e-3 nxpix=20 nypix=20 binsz=0.1 coordsys=CEL \
	xref=$(RA) yref=$(DEC) proj=CAR \
	clobber=$(CLOBBER) chatter=4 # > gttsmap.log 2>&1

clean:
	rm *.fits *.log results.dat $(RESULT) $(SIM_EVROOT)_srcIds.txt

