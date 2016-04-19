.. _tools:

Tools
=====

ScienceTools are described `here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/references.html>`__

* ``enrico_help`` : print help
* ``enrico_setupcheck`` : check your installation
* ``enrico_download`` : download data and auxiliary files (backgrounds, catalog, etc...)
* ``enrico_config`` : produce a configuration file
* ``enrico_gui`` : run the GUI of enrico

* ``enrico_xml`` : produce an xml file that is use to model the ROI using the 3FGL (default) catalog. it can use the 2FGL and the 1FHL.
* ``enrico_sed`` : Run gtlike afer having produced all the need fits files is asked.
* ``enrico_testmodel`` : compute the log(likelihood) of the models `POWERLAW`, `LogParabola` and `PLExpCutoff`.
   An ascii file is then produced in the Spectrum folder with the value of the log(likelihood) for each model.
   You can then use the Wilk's theorem to decide which model best describe the data.
* ``enrico_lc`` : produce a light-curve by running gtlike in different time bins
* ``enrico_foldedlc`` : produce a folded light-curve by running gtlike in different time bins. see the specific section in the config file
* ``enrico_applc``: produce a light-curve using the aperture photometry technique (see `here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html>`__)
* ``enrico_tsmap`` : produce a TS map with or without the source of interest.
* ``enrico_plot_sed`` : plot the SED resulting from enrico_sed
* ``enrico_plot_lc`` : plot the LC resulting from enrico_lc
* ``enrico_plot_tsmap`` : plot the TS Map resulting from enrico_tsmap
* ``enrico_scan`` : make a profile likelihood for each free parameter of the target
* ``enrico_findsrc`` : run gtfindsource
* ``enrico_contour`` : make a confidence contour of 2 parameters
* ``enrico_lrt`` : test custom spectral shapes by calculating their likelihoods
