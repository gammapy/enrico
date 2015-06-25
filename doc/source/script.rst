.. _script:

Scripts
========

Few nice and helpful scripts are available. They required the additional python packages and are located under the script folder.


Plot Count Maps
------------------
The script is plotMaps.py and allows count, model and residuals map to be plotted. The user must provide a configuration file.

.. code-block:: ini

  python plotMaps.py config [vmax]

Options:

 * vmax is used to defined the max of the color scale



.. figure::  _static/Maps.png
   :align:   center

   Maps of  PG 1553+113, from top to bottom: counts map, model map, residuals map.


Plot TS Maps
------------------

The script is plotTSMaps.py. The user must provide a configuration file.

.. code-block:: ini

  python plotTSMaps.py  config [vmax] [stretch]

Options:

 * vmax is used to defined the max of the color scale

 * stretch can be linear (default), log, sqrt, arcsinh or power


.. figure::  _static/TSMaps.png
   :align:   center

   TS Map of PG 1553+113.

Convert Time
------------------

The script ConvertTime.py convert a given time to another units

Type can be : MJD, JD, MET or ISO

.. code-block:: ini

  python ConvertTime.py  type value

Here are some examples:

.. code-block:: ini

    python  ConvertTime.py  MJD 56101
    python  ConvertTime.py  ISO 2012-02-05
    python  ConvertTime.py  ISO 2012-02-05 16:15:18
