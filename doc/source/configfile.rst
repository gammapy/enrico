.. _configfile:

Configuration Files
===================


Generalities
------------

The config file is there to defined the parameters of your analysis. A config file is divided into sections. Each section is a logical group of option and there is some constants across the sections like names, functionalities.


This page will go through all the section, one by one. As for the :ref:`tutorial`, we will asume that the working directory is  `~/myanalysis`.

The `out` option gives the main folder where enrico will put the produced files and the results. The `verbose` option (default is 'yes') allow enrico to print information on the main output like TS, fluxes, predicted numbers of gammas, etc... Such values need computation time which can be saved by turning the option to 'no'. The `Submit` option, if turn to 'yes', will send the jobs to a cluster. `clobber` is a wrapper for the ST tool option and allow to overwrite existing output files is 'yes' (default)

.. code-block:: ini

   out = ~/myanalysis
   verbose = yes
   Submit = no
   clobber = no

Target
------

The target options gives the name, position and spectral model you want to use for the source you are interested in.

Note :

 * Coordinate are in degrees.
 * Available models are 'PowerLaw', 'PowerLaw2', 'LogParabola', 'PLExpCutoff'

The parameter `spectrum` is used for the generation of the sky model. Of course, you can change the model describing the spectrum of any sources (by hand...). All the models supported by the ST are described `here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Likelihood/Model_Selection.html>`_. For now, the supported models by `enrico` are PowerLaw, PowerLaw2, LogParabola, PLExpCutoff, Generic. `Generic` is for the non-supported models. By supported, we mean that enrico can produce a xml file with such model and has some optional features. It is completely possible to run an analyze with enrico and a non-supported model.

.. code-block:: ini

   [target]
      name = PG155+113
      ra = 238.92935
      dec = 11.190102
      spectrum = PowerLaw2

Space
-----

The space option defined the properties of the region of interest (ROI). The names are the same as for the ScienceTools.
The center of the ROI is xref, yref and the size is rad (in degrees). By default xref and yred are the target coordinates but this might be changed.

.. code-block:: ini

   [space]
      xref = 238.92935
      yref = 11.190102
      rad = 10.0
      binsz = 0.1
      coordsys = CEL
      proj = AIT
      phibins = 0

.. note:: 
   if you used binned analysis or for the generation of a TS map, the ROI
   will be divided in nxpix and nypix pixel with a bin size of binsz.

File
----

Here you defined where the FT2 and FT1 files are (this can be an ascii list). The xml option gives the location of the sky model. The `tag` is added to all files produced by enrico. When generating a config file, the value of spacecraft and event are set in such way that is point towards the the file downloaded with `enrico_download`


.. code-block:: ini

   [file]
      spacecraft = ~/myanalysis/FT2.fits
      event = ~/myanalysis/data.list
      xml = ~/myanalysis/XML_model.xml
      tag = MyTag


Time
----

Start and stop time of your analysis in MET. The file option allow the analysis (Lc or SED) to be performed in disjoint time bins. This can be useful for e.g. MWL campaigns or non-constant time bins LC. The file must be an ascii file with 2 columns (start and stop) and each line is a time bin


.. code-block:: ini


   [time]
      tmin = 239557417.0
      tmax = 256970880.0
      file = ""
      type = 'MET'

Energy
------

Minimal and maximal energy of your analysis in MeV. `enumbins_per_decade` is the number of bins per decade for the BINNED analysis chain.


.. code-block:: ini

   [energy]
      emin = 200.0
      emax = 300000.0
      enumbins_per_decade = 10



Environ
-------

Here are defined some directories. They are also defined as environment variables which can be over-writted using the configuration file.

.. code-block:: ini

   [environ]
      # Analysis environment configuration
      # Can also be done via shell environment variables
      FERMI_DATA_DIR = ""
      FERMI_CATALOG_DIR = ""
      FERMI_CATALOG = ""
      FERMI_DIFFUSE_DIR = ""
      FERMI_PREPROCESSED_DIR = ""


Analysis
--------

This part is used to defined how enrico should select the event. You can defined the event class (evclass : 1, 2 , etc..), the zenith angle cut (zmax) and the filter for gtmktime (filter). Also the IRFS used to describe the instrument are defined here (irfs). 

Convtype is use to select either the front (0), back (1) or both (-1) events. If convtype =0 or 1, an ::FRONT of ::BACK is happened at the end of the irfs string automatically allowing to use the good IRFS.

.. code-block:: ini

   [analysis]
      # General analysis options
      likelihood = binned
      evclass = 2
      zmax = 100.0
      roicut = no
      filter = DATA_QUAL==1&&LAT_CONFIG==1&&ABS(ROCK_ANGLE)<52
      irfs = P7REP_SOURCE_V15
      # if convtype =0 or 1, an ::FRONT of ::BACK is happend at the end of the irfs string automatically
      convtype = -1


fitting
-------

Option for the minimizer. You can use MINUIT, NEWMINUIT, DRMGB, etc. ftol is the tolerance that the minimizer should reach.

.. code-block:: ini

   [fitting]
      optimizer = MINUIT
      ftol = 1e-06


model
-----

This section is about the sky model generation. If you have set correctly you environment variables, then enrico is able to find the galactic and extragalactic model. If you want to use other model, you can specify here, their names and locations.

The 3FGL is used by default to find the source in the ROI. All the source with a significance greater than `min_significance` will be added. All sources within `max_radius` (in degrees) have their parameters free to vary in the fitting procedure. The other sources have their parameters frozen to the 3FGL value. 

You can use also the 2FGL or the 1FHL by specifying their name and location.

.. code-block:: ini

   [model]
      # The following options determine the xml model
      diffuse_gal_dir = ""
      diffuse_iso_dir = ""
      diffuse_gal = gal_2yearp7v6_v0.fits
      diffuse_iso = iso_p7v6source.txt
      
      # user points sources for diffuse catalog sources
      point_only = True
      # freeze spectral parameters for weak and far away sources:
      min_significance = 4.0
      max_radius = 3.0



Spectrum
--------

Options for `enrico_sed` which run all the ST tool to make an pointlike analysis.

 * FitsGeneration, if yes, enrico will make all the steps before running gtlike and generated all the fits files needed. If the files have already been generated, change FitsGeneration to no and enrico will only run gtlike

 * ResultPlots : Compute the SED (butterfly) and the model map (in the case of an binned analysis)

 * FrozenSpectralIndex : froze the spectral index of the source (works for POWERLAW and POWERLAW2 models)

 * SummedLike : you can use the summed likelihood method, then front and back event are treated separately and the likelihood which is minimized is the the sum of the front likelihood and back likelihood. This feature is provided by the ScienceTools.

 * Submit : submit the job to a cluster or run it in the current shell.

.. code-block:: ini

   [Spectrum]
      #Generates fits files or not?
      FitsGeneration = no
      #Generates plots (SED, model map)
      ResultPlots = yes
      #Freeze the spectral index of the source
      FrozenSpectralIndex = 0.0
      #Use the summed likelihood method
      SummedLike = no


UpperLimit
----------

This section allows to set up the upper limit computation. During the
computation, the spectral index of the source (it is assumed that a POWERLAW or
POWERLAW2 model is used) is frozen to `SpectralIndex`. Two methods can be used,
Profile of Integral, see the Fermi web site for more informations.

An upper limit, at the confidence level `cl`, is computed if the TS is below TSlimit. This hold only for `enrico_sed`


.. code-block:: ini

   [UpperLimit]
      #Assumed Spectral index
      SpectralIndex = 1.5
      # UL method could be Profile or Integral (provided by the fermi collaboration)
      Method = Profile
      envelope = no
      #Compute an UL if the TS of the sources is <TSlimit
      TSlimit = 25.0
      # Confidence level for the Ul computation
      cl = 0.95

LightCurve
----------

Option for enrico_lc which run an entire analysis in time bins and produce all the fits files needed to use gtlike.

 * FitsGeneration, if yes, enrico will make all the steps before running gtlike and generated all the fits files needed. If the files have already been generated, change FitsGeneration to no and enrico will only run gtlike

 * NLCbin : number of time bins

 * MakeConfFile : enrico_lc will produce config file readable by enrico for each time bin. You can ask the tool to not do so, if you want to use/modify the config files.

 * Submit : submit the job to a cluster or run it in the current shell.

 * TSLightCurve : an upper limit is computed is the TS in a time bin is below this value.

 * DiagnosticPlots : ask enrico_plot_lc to generate diagnostic plot (TS vs time, Npred vs flux ...)

.. code-block:: ini

   [LightCurve]
      #Generates fits files or not?
      FitsGeneration = yes
      #Number of points for the LC
      NLCbin = 20
      MakeConfFile = no
      #Compute an UL if the TS of the sources is <TSLightCurve
      TSLightCurve = 9.0
      #Generates control plots
      DiagnosticPlots = yes


Folded LightCurve
-----------------

This section is devoted to the folded LC. This is designed for binary system analysis.

  * NLCbin : number of time bins

  * epoch: Epoch of phase=0 in MJD, equal to tmin is 0

  * Period: Orbital period in days

.. code-block:: ini

   [FoldedLC]
      #Number of bins for the orbitally folded LC
      NLCbin = 10
      #Epoch of phase=0 in MJD, equal to tmin is 0
      epoch = 0
      #Orbital period in days
      Period = 10


Ebin
----

 * FitsGeneration, if yes, enrico will make all the steps before running gtlike and generated all the fits files needed. If the files have already been generated, change FitsGeneration to no and enrico will only run gtlike

 * NumEnergyBins :  number of bins in energy

 * TSEnergyBins : an upper limit is computed is the TS in an energy bin is below this value.

 * Submit : submit the job to a cluster or run it in the current shell.

.. code-block:: ini

   [Ebin]
      #Generates fits files or not?
      FitsGeneration = yes
      NumEnergyBins = 7
      #Compute an UL if the TS of the sources is <TSEnergyBins
      TSEnergyBins = 9

Option for enrico_tsmap

TSMap
-----

This section is used to configured `enrico_tsmap` and `enrico_plot_tsmap` 

 * Re-Fit : use rerun gtlike in order to have the best fit parameters in your model.

 * npix : number of pixels of you map. Remember that the TS map grid is based on the other maps (like count map) produced before and centred to the coordinates xref,yref.

 * RemoveTarget : remove your source of interest form the map by freezing its parameters.

 * Submit : submit the job to a cluster or run it in the current shell.

In order to speed up the process, parallel computation can be used. Either each pixel can be a job by itself (option [TSMap]/method = pixel) or a job can regroup an entire row of pixel (option [TSMap]/method = row)

.. code-block:: ini

   [TSMap]
      #Re-fit before computing the TS map
      Re-Fit = no
      #Numbers of pixel in x and y
      npix = 10
      #Remove or not the target from the model
      RemoveTarget = yes
      #Generate the TS map pixel by pixel or by grouping the pixels by row.
      #(reduce the numbers of jobs but each job are longer)
      method = row


If a pixel (or a row) has failed you can rerun it. For the pixel 49,4 :

.. code-block:: ini

   enrico_tsmap myanalysis.conf 49 4


For the entire row 49 :

.. code-block:: ini

   enrico_tsmap myanalysis.conf 49



Finding the position of a source
--------------------------------

This section is used to configured `enrico_findsrc`. It run the tool gtfindsource and update the file Roi_model.reg with the fitted position in red.

 * FitsGeneration, if yes, enrico will make all the steps before running gtfindsource and generated all the fits files needed. If the files have already been generated, change FitsGeneration to no and enrico will only run gtfindsource

 * Refit :  re-run the optimizer before (use the option reopt)

.. code-block:: ini

   [findsrc]
      #Generates fits files or not?
      FitsGeneration = option('yes', 'no', default='yes')
      #Reoptimize before
      Refit = option('yes', 'no', default='yes')

