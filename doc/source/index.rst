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

Enrico is based on a configuration file which contains all the setup for your analysis. For each enrico tool, you just have to type

.. code-block:: bash

   enrico_'tool' Configuration.conf

A simple analysis
--------

After you have done the :ref:`setup`, running an analysis is as simple as
typing these commands:

.. code-block:: bash

   enrico_config Myconfig.conf <answer a few questions like the position of your target>

It will create an file named Myconfig.conf. Be sure that it contains all your setup. After you should generate an xml file for the sky model : 

.. code-block:: bash

   enrico_xml Myconfig.conf


Now run the likelihood analysis for the full spectrum and data points:

.. code-block:: bash

   enrico_sed Myconfig.conf

If you like,  a light curve or a TS map:

.. code-block:: bash

   enrico_lc Myconfig.conf
   enrico_tsmap Myconfig.conf

Note that if you have a computing cluster at your disposal, these
compute-time-intensive tasks can be run in parallel.

Finally at the end you should plot your results using: 

.. code-block:: bash

   enrico_plot_sed Myconfig.conf (plot the SED)
   enrico_plot_lc Myconfig.conf (plot the lightcurve)
   enrico_plot_tsmap Myconfig.conf (generate a fits file for the TS map)


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

