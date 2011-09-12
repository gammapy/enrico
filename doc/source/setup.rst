.. _setup:

Setup
=====

The most difficult tasks is how to set up your analysis
environment correctly. Once that's done running the
analysis will be simple, because Enrico will help you.

You need to:

* Install the Fermi ScienceTools
* Install additional python packages
* Install Enrico
* Download the diffuse model and 2FGL catalog files
* Download Fermi data (photon and spacecraft files)

Here we give some instructions on how to do these steps,
but to be honest, if you don't know how to install python
packages (e.g. know what PYTHONPATH, setup.py and pip are),
you might fail and will have to ask someone at your institute for help.

Install the Fermi ScienceTools
------------------------------

Download and install the Fermi Science Tools as described 
`here <http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/>`__.

Then do this to set up your shell (we are assuming `bash` throughout) for Fermi analysis:

.. code-block:: bash

   export FERMI_DIR = <...>
   source $FERMI_DIR/bin/init_fermi.sh

Check that you have access to the Fermi command line tools and python packages:

.. code-block:: bash

   gtirfs
   python
   >>> import UnbinnedAnalysis

Install additional python packages
----------------------------------

.. note::
   You don't have to install all of the following packages,
   but if you do you'll have a much nicer and more powerful
   python environment.
   
   configobj is used throughout and you really need it,
   other packages are only used in some parts (e.g. uncertainties
   to compute the spectrum butterfly or minuit for multiwavelength
   spectral fitting) and you can use the rest of Enrico without it.

   You'll get an `ImportError` with the name of the missing package
   once you try to use part of the code that relies on that package.

First of all you should install `distribute <http://packages.python.org/distribute/>`__ 
and `pip <http://www.pip-installer.org>`__ as discribed
`here <http://www.pip-installer.org/en/latest/installing.html>`__, because
pip makes it easy to install additional packages:

.. code-block:: bash

   curl -O http://python-distribute.org/distribute_setup.py
   python distribute_setup.py
   curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
   python get-pip.py
   which pip (should be located in the Fermi software)
   pip (should print a help message)

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

Install Enrico
--------------

Get the enrico package and build the documentation if you like:

.. code-block:: bash

   git clone ~deil/git/enrico
   cd doc
   make html
   firefox build/html/index.html

Make sure your PATH contains the enrico command line tools (in the scripts directory)
and PYTHONPATH contain the enrico python package:

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:<...>
   export PATH=$PATH:<...>

If you like, you can also install Enrico into your `site-packages`:

.. code-block:: bash

   cd enrico
   python setup.py install --user 
   (or however you install python packages on your system)

Run the following command to check the status of your analysis environment:

.. code-block:: bash

   enrico_setup --check

Download the diffuse model and 2FGL catalog files
-------------------------------------------------

Enrico uses the following environment variables to find
the catalog and diffuse model files

.. code-block:: bash

   FERMI_CATALOG_DIR
   FERMI_DIFFUSE_DIR

After setting these you can run this command to download any missing files from the 
You can use the `enrico_check` tool to check that you have all
you need and to download missing files from the 
`FSSC <http://fermi.gsfc.nasa.gov/ssc/>`__:

.. code-block:: bash

   enrico_setup --download

Download Fermi data (photon and spacecraft files)
-------------------------------------------------

There are two options. If you are only analyzing one ore two
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
   
and then run

.. code-block:: bash

   enrico_setup --update_weekly

This will run wget to update only the weekly files that are necessary and download a 
spacecraft file for the whole mission (~ 500 MB). There is no documented method
to combine weekly spacecraft files.

Obviously you should share one software and data installation per institute and
not hit the FSSC servers without need.

Issues
------

* Building from source doesn't work on the MPIK cluster or on my Mac.

* On my Mac importing pyIrfLoader fails if pyLikelihood hasn't been
  imported first. So if you ever see that error, look at the
  traceback where it happens and replace

.. code-block:: python

   >>> import pyIrfLoader

   with 
   
.. code-block:: python

   >>> import pyLikelihood      
   >>> import pyIrfLoader

* On the MPIK cluster we have SUSE Linux, which is not officially supported by Fermi.
  Nevertheless is is possible to use the binary, but there are a few shared libs for python missing,
  which I soft-linked to the system (64 bit) versions:

.. code-block: bash

   deil@lfs1:/lfs/l1/hess/users/deil/bin/HESSSurveyBin/fermi/latest/x86_64-unknown-linux-gnu-libc2.5/lib/python2.6/lib-dynload> pwd
   /lfs/l1/hess/users/deil/bin/HESSSurveyBin/fermi/latest/x86_64-unknown-linux-gnu-libc2.5/lib/python2.6/lib-dynload
   lrwxrwxrwx 1 deil hfm   40 Sep 12 02:23 _md5.so -> /usr/lib64/python2.6/lib-dynload/_md5.so
   lrwxrwxrwx 1 deil hfm   40 Sep 12 02:24 _sha.so -> /usr/lib64/python2.6/lib-dynload/_sha.so
   lrwxrwxrwx 1 deil hfm   43 Sep 12 02:24 _sha256.so -> /usr/lib64/python2.6/lib-dynload/_sha256.so
   lrwxrwxrwx 1 deil hfm   43 Sep 12 02:24 _sha512.so -> /usr/lib64/python2.6/lib-dynload/_sha512.so

To build numpy more the Fermi readline away, then

.. code-block: bash
   python setup.py build --fcompiler=gnu95

* I did not manage to install ipython, because it uses sqlite3, which is also missing from the
  Fermi binary installation, and the system _sqlite.so is also incompatible because it uses another unicode encoding.

.. code-block:: bash
 
  ImportError: /lfs/l1/hess/users/deil/bin/HESSSurveyBin/fermi/latest/x86_64-unknown-linux-gnu-libc2.5/lib/python2.6/lib-dynload/_sqlite3.so: undefined symbol: PyUnicodeUCS4_DecodeUTF8
