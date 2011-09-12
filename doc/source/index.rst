Welcome
=======

..
   I could not make this image apprear besides the welcome text,
   so I included it as html_logo in conf.py for now.
   .. |enrico| image:: _static/Enrico.jpg
      :width: 200
      :align: bottom

Hi, I am Enrico, and I will help you run your 
`Fermi <http://fermi.gsfc.nasa.gov/>`_ data analysis!

`Fermi <http://fermi.gsfc.nasa.gov/>`_

After you have done the :ref:`setup`, running an analysis is as simple as
typing these commands:

.. code-block:: bash

   enrico_config <answer a few questions like the position of your target>
   enrico_xml <answer a few questions about how you'd like to model sources in your ROI>

Now you have to get the data from the 
`FSSC LAT Data Server <http://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi>`_
or if you set up weekly files you can apply selections once so that subsequent analysis will be faster

.. code-block:: bash

   enrico_get_data (runs gtselect on the weekly photon files)

Now run the likelihood analysis for the full energy band:

.. code-block:: bash

   enrico_like (runs the analysis for the full energy band)

If you like, you can additionally compute flux points, a light curve,
a TS map:

.. code-block:: bash
   
   enrico_fluxpoints (run gtlike in energy bins) 
   enrico_lightcurve (run gtlike in time bins)
   enrico_tsmap (make a TS or residual TS map)

Note that if you have a computing cluster at your disposal, these
compute-time-intensive tasks can be run in parallel.

Finally at the end you should run the summary command to make
inspecting your results easy:

.. code-block:: bash

   enrico_summary (make plots and a html summary of your analysis)

Features
--------

* Get your results easy and fast by using the enrico command line tools.
* Results are reproducible, because config files and logs are used. 
* The enrico command line tools are just frontends for functions and
  classes in the enrico python package, so if you know the Fermi tools
  and some python it is easy for you to modify things to your needs.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   setup
   tutorial
   sed_fit
   developer

.. code-block:: python


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

