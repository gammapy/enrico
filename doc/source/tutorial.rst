.. _tutorial:

Tutorial
========

This tutorial will walk you through the steps to repeat the analysis
of the AGN TODO done in the official `Fermi collaboration python tutorial
<http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/python_tutorial.html>`__
and results published here.

We strongly recommend that you run this analysis once, so that you
can be sure to produce correct results before analyzing your source.

First init the Fermi and Enrico tools

.. code-block:: bash

   source ~deil/code/enrico/init.sh

Then run this check to make sure your setup is ok:

.. code-block:: bash

   enrico_setup --check
   
Enrico will tell you if you forgot something, 
in which case you should add it now.

Make a new directory now where you want all of your intermediate
and final analysis results to be placed and go there. 
For illustration purposes we will assume in this tutorial this
directory is `~/myanalysis`:

.. code-block:: bash

   mkdir ~/myanalysis
   cd ~/myanalysis

Make a config file
------------------

You can use the `enrico_config` tool to quickly make a config file
called `enrico.conf`. It will ask you for the required options
and copy the rest from a default config file `enrico/cata/config/default.conf`:

.. code-block:: bash

   $ enrico_config enrico.conf
   Please provide the following required options:
   Output directory: .
   Right Ascension: 43
   Declination: 44

Now you can edit this config file by hand to make further adjustments.

.. note:: 
   If you know exactly how the analysis steps work you can also make
   adjustments later on. But we have not put in a gadzillion of
   checks for each step to make sure that parameters are consistent
   with previous steps, so it is best to only adjust parameters
   at the beginning.

Make a model xml file
---------------------

Next you run enrico_xml to make a model xml file.

.. code-block:: bash

   $ enrico_xml enrico.conf 
   This is make2FGLxml version 04.
   Writing ./model.xml
   Added 29 point sources, note that any extended sources in ROI were modeled as point sources becaue psForce option was set to True

Note that you give options for this step simply by mentioning
the config file. For the `enrico_xml` tool, the following options
are relevant (@todo: the way this works is work in progress!):

.. code-block:: ini

   [target]
      ra = 43.0
      dec = 44.0
      # Target and modelling options
      # If you put a 2FGL source name here, you don't have to specify the position
      name = MySource
      spectrum = PL
      # If the source of interest corresponds to one or more 2FGL sources,
      # list them here (comma separated) to not include them in the xml model
      2FGL_Sources = ,

   [model]
      # The following options determine the xml model
      diffuse_gal = gal_2yearp7v6_v0
      diffuse_iso = iso_p7v6source
      
      # user points sources for diffuse catalog sources
      point_only = True
      # freeze spectral parameters for weak and far away sources:
      min_significance = 4.0
      max_radius = -1.0

Get data
--------

There are two possibilities:

* Download data by hand for this target.
* If weekly files are set up, just extract data from there automatically
  as a first step in `enrico_like`.

Run global fit
--------------

Make flux points
----------------

Make a light curve
------------------

Check results
-------------

