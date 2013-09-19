.. _tools:

Description of each tool
========================

ScienceTools are described `here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/references.html>`__

enrico_setupcheck : check your installation

enrico_download: download data and auxiliary files (backgrounds, catalog, etc...)

enrico_config: produce a configuration file

enrico_xml: produce an xml file that is use to model the ROI using the 2FGL catalog.

enrico_sed: Run gtlike afer having produced all the need fits files is asked.

enrico_testmodel : compute the log(likelihood) of the models `POWERLAW`, `LogParabola` and `PLExpCutoff`. An ascii file is then produced in the Spectrum folder with the value of the log(likelihood) for each model. You can then use the Wilk's theorem to decide which model best describe the data.

enrico_lc: produce a light-curve by running gtlike in different time bins

enrico_applc: produce a light-curve using the aperture photometry technique (see `here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html>`__)

enrico_tsmap : produce a TS map with or without the source of interest.

enrico_plot_sed: plot the SED resulting from enrico_sed

enrico_plot_lc: plot the LC resulting from enrico_lc

enrico_plot_tsmap: plot the TS Map resulting from enrico_tsmap

enrico_help: print help
