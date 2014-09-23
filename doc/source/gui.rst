.. _gui:

Graphical User Interface (GUI)
==============================

This page describe the windows of the `enrico` GUI. To run the GUI, type 

.. code-block:: bash

   enrico_'gui' Configuration.conf

The GUI aims to allow easy configuration file management and to run tools. Each page, arranged in tabs, roughly correspond to a section of the configuration file (see :ref:`configfile`). At the bottom of the GUI are buttons to save the configuration file and close the GUI. The convention is such that a crossed box stand for a `yes` in the configuration file.

Main page
---------

The first page is the Main page of the GUI. Here the main options can be defined and tools can be run using the buttons (run a tool save the configuration file on disk).

.. figure::  _static/GuiMain.png
   :align:   center


Files page
----------

The second page manage the files definition (event, FT2 and xml) as well as the tag of your analysis. 


.. figure::  _static/GuiFile.png
   :align:   center

Target/Space page
-----------------

The target (name, position and spectral model) is defined is this page (first frame). The second frame defined the ROI (center  (Xref,  Yref), size). The `sync` button update the value Xref and Yref with the target position. It is then possible to have a ROI not centered on the target. The projection type and coordinate system is defined here also.

.. figure::  _static/GuiTargetSpace.png
   :align:   center

Analysis page
-------------

The analysis tab deal with the analysis chain (binned or unbinned), some cuts (zenith angle (zmax), filter for the GTI definition). The IRFs are also defined here.

The fitting algorithm (MINUIT, NEWMINUIT, etc) and the tolerance are setup here.

.. figure::  _static/GuiAnalysis.png
   :align:   center


Energy/Time page
----------------

This page define the energy and time ranges.

.. figure::  _static/GuiET.png
   :align:   center


Spectrum/Ebin page
------------------

This page is used to manage the spectrum and energy bins generation as in the configuration file. The buttons `Re-run Ebin` can be used to only rerun the bin (by running as many jobs as the number of bin) 


.. figure::  _static/GuiEbin.png
   :align:   center

   Spectrum page of the GUI

Upper Limits page
-----------------

This page allows the definition of the UL computation parameters : assumed index, minimal TS  and confidence level

.. figure::  _static/GuiUL.png
   :align:   center


Light Curves page
-----------------

The light Curves are setup here.

.. figure::  _static/GuiLC.png
   :align:   center



Aperture/FoldedLC page
----------------------

The first frame of the page is for the aperture LC and the second for the Folded LC.

.. figure::  _static/GuiAppLC.png
   :align:   center


TS Map page
-----------

TS Map parameters are managed here

.. figure::  _static/GuiTSMap.png
   :align:   center


Findsrc/srcprob page
--------------------


The findsrc and the srcprob tools parameters are managed here

.. figure::  _static/GuiFindSRC.png
   :align:   center

Plots page
----------

This page allow to draw the produced plots. Using the buttons, you have access to 

* the SED and corresponding debug plots
* the LC and corresponding debug plots

If the plot has not been produced, Enrico Fermi picture is display.

.. figure::  _static/Guiplot.png
   :align:   center


