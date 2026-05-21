prep_site_forcing
=================

Prepare input data to run ecosystem models.

|License|


About
-----

The scripts of ``prep_site_forcing`` take input data, gap-fill it, and
produce files in the format and units suitable to run specific
ecosystem models. They can take data from ascii input files, download
data from the ICOS carbon portal, or can use ERA5 and ERA5-Land
data. Data can be imputed (gap-filled) by linear interpolation or
using bias-corrected ERA5(-Land) data.

The script reads information from a configuration file and produces a
csv file with imputed data. The script then runs a model specific
version of `ascii2netcdf`, e.g. `ascii2musica`, that produces forcing
data for the ecosystem model in the very model-specific
form. Currently there are the two model-specific scripts
`ascii2isba.py` and `ascii2musica.py` extending the base class
`ascii2Netcdf` from `class_ascii2netcdf.py`.

The driver script `prep_site_forcing.py` is controlled by a
configuration file, which is in the simple Python configparser format.


Installation
------------

There is no installation. Simply clone the repository:

.. code-block:: bash

   git clone https://github.com/mcuntz/prep_site_forcing.git

and run the script:

.. code-block:: bash

   python prep_site_forcing.py -p bias_correction.pdf -t pdf FR-Hes.cfg

Requirements
   * cftime_
   * matplotlib_
   * netcdf4_
   * numpy_
   * pandas_
   * pyjams_
   * scipy_
   * xarray_
as well as
   * cdsapi_  if ERA5 must be downloaded
   * icoscp_core_  if ICOS data is used







License
-------

``prep_site_forcing`` is distributed under the MIT License. See the LICENSE_ file
for details.

Copyright (c) 2026- Matthias Cuntz


.. |License| image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/mcuntz/prep_site_forcing/blob/main/LICENSE

.. _cftime: https://github.com/Unidata/cftime
.. _matplotlib: https://matplotlib.org/
.. _netCDF4: https://github.com/Unidata/netcdf4-python
.. _numpy: https://numpy.org/
.. _pandas: https://pandas.pydata.org/
.. _pyjams: https://github.com/mcuntz/pyjams
.. _scipy: https://scipy.org/
.. _xarray: https://xarray.dev/
.. _cdsapi: https://cds.climate.copernicus.eu/how-to-api
.. _icoscp_core: https://icos-carbon-portal.github.io/pylib/icoscp/getting_started/
