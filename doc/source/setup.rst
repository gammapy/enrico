.. _setup:

Setup
=====

The most difficult tasks is how to set up your analysis
environment correctly. Once that's done running the
analysis will be simple, because Enrico will help you.

You need to:

* Install the Fermi ScienceTools
* Install Enrico
* Install additional (optional) python packages
* Download Fermi data (photon and spacecraft files)
* Download the diffuse model and 3FGL catalog files


Here we give some instructions on how to do these steps,
but to be honest, if you don't know how to install python
packages (e.g. know what PYTHONPATH, setup.py and pip are),
you might fail and will have to ask someone at your institute for help.

Install the Fermi ScienceTools
------------------------------

The first step is to get the Fermi Science Tools

Download and install the Fermi Science Tools as described 
`here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/>`__. Best option
is to download the binnary and unzip them


Set up your shell (we are assuming `bash` throughout) for Fermi analysis. The
first thing to do is to setup a environment variable 'FERMI_DIR'. This variable
must point to the folder where the 'init_fermi.sh' script is. For exemple for
the v10r0p5 version of the ST and the SL6 binary you have :


.. code-block:: bash

   export FERMI_DIR = /Where/I/put/the/ST/ScienceTools-v10r0p5-fssc-20150518-x86_64-unknown-linux-gnu-libc2.12/x86_64-unknown-linux-gnu-libc2.12
   source $FERMI_DIR/init_fermi.sh

An 'ls' command in the FERMI_DIR should give you something like :

.. code-block:: bash

    BUILD_DIR  build  etc             fermi-init.sh  help   include  macros  refdata  syspfiles  tutorials
    bin        cint   fermi-init.csh  fonts

Check that you have access to the Fermi command line tools and python packages:

.. code-block:: bash

   gtirfs

If there is no error message, the ST are installed.

.. code-block:: bash
		
   python
   >>> import UnbinnedAnalysis

If there is no error message, the python support of ST is installed.

Install Enrico
--------------

Get the enrico package 

.. code-block:: bash

   git clone https://github.com/gammapy/enrico.git


Again you have to setup an environment variable name ENRICO_DIR which point to
the enrico folder. After the clone type


.. code-block:: bash
   
   cd enrico
   export ENRICO_DIR=$PWD

or alternatavely 

.. code-block:: bash

   export ENRICO_DIR=/Where/I/put/Enrico/enrico/


An 'ls' command in the ENRICO_DIR should give you something like :

.. code-block:: bash

    CHANGES.txt  LICENSE.txt  README.rst  bin  doc  enrico  enrico-init.csh  enrico-init.sh  script

The last step is to source the init file:

.. code-block:: bash

   source $ENRICO_DIR/enrico-init.sh


This command will setup you PATH and PYTHONPATH variable to have access to the
enrico tools. Run the following command to check the status of your analysis
environment:

.. code-block:: bash

   enrico_setupcheck


Build the documentation if you like:

.. code-block:: bash

   cd doc
   make html
   firefox build/html/index.html

Install additional (optional)  python packages
----------------------------------------------

.. note::
   You don't have to install all of the following packages,
   but if you do you'll have a much nicer and more powerful
   python environment.
   
   configobj is used throughout and you really need it,
   other packages are optional or come with the ST

   You'll get an `ImportError` with the name of the missing package
   once you try to use part of the code that relies on that package.

First of all you should install `setuptools <http://pythonhosted.org//setuptools/>`__ 
and `pip <https://pip.pypa.io/en/latest/>`__ as described
`here <https://pip.pypa.io/en/latest/installing.html#install-pip>`__, because
pip makes it easy to install additional packages. To install both just run:

.. code-block:: bash

   curl -O https://bootstrap.pypa.io/get-pip.py
   python get-pip.py
   which pip # should be located in the Fermi software
   pip # should print a help message

Next install `ipython <http://ipython.org/>`__, which is a much nicer interactive 
python shell than the default python shell and 
`configobj <http://www.voidspace.org.uk/python/configobj.html>`__,
which is a more powerful config file reader and is user
by Enrico instead of the `ConfigParser <http://docs.python.org/library/configparser.html>`_ 
from the python standard library. `nose <http://readthedocs.org/docs/nose/en/latest/>__
is a python test runner, used e.g. by `numpy.test()`. `Sphinx <http://sphinx.pocoo.org/>`__
is the python documentation generator and we also use it for this project:

.. code-block:: bash

   pip install ipython
   pip install configobj
   pip install nose
   pip install sphinx
   
Now update to a recent `Numpy and Scipy <http://www.scipy.org/>`__. The Fermi tools
ship with a very old Numpy (version 1.4.1) and no Scipy (even though
scipy is used e.g. in `IntegralUpperLimits.py`.

.. code-block:: bash

   pip install numpy
   pip install scipy

.. note::
   Numpy and Scipy have many C and Fortran extensions and compiling
   those can fail. In that case you have to download the packages
   and build them yourself, adjusting some build options to your system.

   .. code-block:: bash
   
      git clone https://github.com/numpy/numpy/
      cd numpy
      python setup.py build <options for your system here>

Finally install some nice and useful python packages:

* `Kapteyn <http://www.astro.rug.nl/software/kapteyn-beta/>`__
  is great for working with coordinates and plotting images,
* `ATpy <http://atpy.github.com/>`__
  has a nicer API for working with tables than pyfits
* `uncertainties <http://packages.python.org/uncertainties/>`__
  makes error propagation dead simple.


.. code-block:: bash

   pip install http://www.astro.rug.nl/software/kapteyn-beta/kapteyn-2.1.1b9.tar.gz
   pip install atpy   
   pip install uncertainties



Download Fermi data (photon and spacecraft files)
-------------------------------------------------

There are two options. If you are only analysing one ore two
targets, you can download the data for these targets specifically
from the `FSSC dat server <http://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi>`__.

If you are doing many analyses or survey work, you should download
the complete data set, i.e. one global spacecraft file and
weekly photon files from the `FSSC FTP server <ftp://legacy.gsfc.nasa.gov/fermi/data/>`__.

Actually Enrico will help you working with the weekly files.
Just set the following environment variable to 
wherever you'd like the spacecraft file and weekly photon files to be:

.. code-block:: bash

   FERMI_DATA = <somewhere with ~20 GB storage space>
   

Then running the following command will download the data in an incremental manner

.. code-block:: bash

   enrico_download --download_data

This will run wget to update only the weekly files that are necessary and download a 
spacecraft file for the whole mission (~ 500 MB). There is no documented method
to combine weekly spacecraft files.

Obviously you should share one software and data installation per institute and
not hit the FSSC servers without need.

Download the diffuse model and 3FGL catalog files
-------------------------------------------------

The diffuse model and 3FGL catalog files can be downloaded from the `FSSC <http://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html>`__

Enrico uses the following environment variables to find
the catalog and diffuse model files

.. code-block:: bash

   FERMI_CATALOG_DIR
   FERMI_DIFFUSE_DIR
   FERMI_DOWNLOAD_DIR
   FERMI_PREPROCESSED_DIR

They are set automatically but you can change the default value and run the following command to download any missing files from the FFSC

.. code-block:: bash

   enrico_download  --download_aux

This will also download the Template files for the analysis of extended sources.


Issues
------

* Building from source doesn't work on the MPIK cluster or on my Mac.

* Importing pyIrfLoader might fail if pyLikelihood hasn't been
  imported first. So if you ever see that error, look at the
  traceback where it happens and replace

.. code-block:: python

   >>> import pyIrfLoader

with 
   
.. code-block:: python

   >>> import pyLikelihood      
   >>> import pyIrfLoader


