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

* Code: https://github.com/gammapy/enrico
* Issues: https://github.com/gammapy/enrico/issues
* Documentation: http://enrico.readthedocs.org/
* Mailing List: http://groups.google.com/group/gammapy_enrico
* Fermi Science Tools: http://fermi.gsfc.nasa.gov/ssc/


Features
--------

* Get your results easy and fast by using the enrico command line tools.
* Results are reproducible, because config files and logs are used. 
* The enrico command line tools are just frontends for functions and
  classes in the enrico python package, so if you know the Fermi tools
  and some python it is easy for you to modify things to your needs.
* A :ref:`gui` has been developed to simply the use of the package.

Enrico is based on a configuration file which contains all the setup for your
analysis. For each enrico tool, you just have to type

.. code-block:: bash

   enrico_'tool' Configuration.conf

Quick setup
-----------

You need to have the Fermi ScienceTools installed on you machine. 

The provided scripts enrico-init.sh will set up enrico for you. Before you need
to define few global variables.

.. code-block:: bash

   export FERMI_DIR=< location of your Fermi software installation >
   export FERMI_DATA_DIR=< location of your Fermi weekly and preprocessed data > 
   # (not mandatory)
   source $FERMI_DIR/fermi-init.sh

   export ENRICO_DIR=< location of your Enrico software checkout >
   source $ENRICO_DIR/enrico-init.sh

Sourcing init_enrico.sh will setup your environment. You can check your
installation by typing 

.. code-block:: bash

  enrico_setupcheck

For more informations, go to :ref:`setup`

A simple analysis
-----------------

After you have done the :ref:`setup`, running an analysis is as simple as
typing these commands:

.. code-block:: bash

   enrico_config Myconfig.conf <answer a few questions like the position of your target>

It will create an file named Myconfig.conf. Be sure that it contains all your
setup. After you should generate an xml file for the sky model : 

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


For more information, see :ref:`tutorial`


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   setup
   tutorial
   configfile
   gui
   tools
   script
   developer
   acknowledge

Other tools and tutorials
-------------------------

You can find further examples how to use Enrico in the `fermi-hero <https://fermi-hero.readthedocs.org/>`__
tutorial given at the `IMPRS 2013 summer school on high-energy astrophysics <http://www.mpia.de/imprs-hd/SummerSchools/2013/>`__.

Other tools to simplify the use of the Fermi ScienceTools have been collected on this page:
`Fermi user-contributed tools <http://fermi.gsfc.nasa.gov/ssc/data/analysis/user/>`__

Other packages useful for gamma-ray astronomy data analysis in general are
`Astropy <http://www.astropy.org/>`__,
`Gammapy <https://gammapy.readthedocs.org/>`__,
`Gammafit <http://gammafit.readthedocs.org/>`__
and `Gammalib <http://gammalib.sourceforge.net/>`__.  

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

