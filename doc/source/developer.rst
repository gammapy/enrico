.. _developer:

Developer Information
=====================

Any improvements to Enrico are welcome!

Please report bugs and contribute features using github.

To be implemented and tested
----------------------------

Here is a list of things that need to be done.

* Access 2FGL and plot lightcurve / spectra with 
  flux points / upper limits / butterflies
* TS map
* Nice analysis result summary web page
* MPIK find data & diffuse emission & catalog files automatically
* Simultaneous Fermi - HESS spectral fit
* Write documentation -> tutorial
* Model the 12 ExtendedSources (extension 4 in the catalog) correctly.
* Include methods to fit Gaussian and disk sources (i.e. make and use FITS template)

Open Issues
-----------

* I'm not sure what the best way is to include data files in
  the python package.
  Maybe use this solution instead, which has the advantage
  that the package doesn't have to be installed?
  http://www.no-ack.org/2010/09/including-data-files-into-python.html

How to contribute?
------------------

Please file bug reports and feature requests on github.

Or even better (for us :-), fork us on github, make the improvement yourself and
make a pull request for your change to be included in the main repo.

If you don't know how to use git and github, check out
http://astropy.readthedocs.org/en/latest/development/
