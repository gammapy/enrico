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

Plot TS Maps
------------------

The script is plotTSMaps.py. The user must provide a configuration file.

.. code-block:: ini

  python plotTSMaps.py  config [vmax] [stretch]

Options:

 * vmax is used to defined the max of the color scale

 * stretch can be linear (default), log, sqrt, arcsinh or power
