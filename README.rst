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

The script reads information and options from a configuration file,
then reads the input data and produces a csv file with imputed
data. The script then calls a model specific version of
`ascii2netcdf`, e.g. `ascii2musica`, that produces forcing data for
the ecosystem model in the very model-specific form. Currently there
are the two model-specific scripts `ascii2isba.py` and
`ascii2musica.py` extending the base class `ascii2Netcdf` from
`class_ascii2netcdf.py`.


Installation
------------

There is no installation. Simply clone the repository:

.. code-block:: bash

   git clone https://github.com/mcuntz/prep_site_forcing.git

and run the script:

.. code-block:: bash

   python prep_site_forcing.py FR-Hes.cfg

You can see the check the results of the linear bias-correction of the
ERA5(-Land) data by plotting them into a file:

.. code-block:: bash

   python prep_site_forcing.py -p bias_correction.pdf -t pdf FR-Hes.cfg

See `python prep_site_forcing.py -h` for help.

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


Config file
-----------

The driver script `prep_site_forcing.py` is controlled by a
configuration file, which is in the simple Python configparser
format. It has sections with case-sensitive names in brackets,
e.g. `[Site]`, and case-insensitive options within each
section. Mandatory sections are `[Model]`, `[Site]`, and
`[VarNames]`. Further, the option `input` in `[Options]` must be set.

See `FR-Hes.cfg` for an example that can be used as a template for
other sites. The example files is highly commented and should be
rather self-explaining. Or see the explanations below.


Section [Model]
^^^^^^^^^^^^^^^

The `model` name must be given. This is mostly used to call the
correct `ascii2netcdf` script.

For example:

.. code-block:: python

   [Model]
   model = MuSICA


Section [Options]
^^^^^^^^^^^^^^^^^

There are several general options.

`input` tells, which data should be taken. It can be `file`, `ICOS`,
or `ERA5`. If it is `file`, the section `[Input]` will be used and
data will be read from a local file. If `input` is set to `ICOS`, the
section `[ICOS]` will be used and data is downloaded from the ICOS
Carbon Portal or read from a local file in ICOS convention. If `input`
is `ERA5`, the section `ERA5` will be used and data downloaded from
the Copernicus archive or read from local files will be used directly
as forcing data. There is no default value.

A Mean Absolute Deviation (MAD) filter can be applied to the input
data to catch outliers. Set `mad_z` greater than zero to apply an MAD
filter with `z = mad_z`. Default is not to apply an MAD filter (`mad_z
= 0`).

The imputation method can be chosen. Currently there is linear
interpolation (`imputation_method = 0`) to fill gaps or using
ERA5(-Land) data that gets bias-corrected with the existing local data
(`imputation_method = 1`; Vuichard and Papale, ESSD 2013,
https://doi.org/10.5194/essd-7-157-2015). Default is
`imputation_method = 1`.

`make_netcdf` controls if a function `ascii2netcdf` from a file
`ascii2netcdf.py` will be called. 'netcdf' is thereby replaced by the
lowercase model name. If `make_netcdf = False`, the script stops after
production of the csv file (this implies that the next option is True,
i.e. `keep_csv = True`). Default is `make_netcdf = True`.

If `keep_csv` is True, then the csv file will not be deleted after
calling `ascii2netcdf`. Default is `keep_csv = False`.

For example:

.. code-block:: python

   [Options]
   input = file
   mad_z = 7
   imputation_method = 1
   make_netcdf = True
   keep_csv = True


Section [Site]
^^^^^^^^^^^^^^

Data in the `[Site]` section is mostly used to write it into the
netcdf files.  The site `name` is basically for information but it is
also used as site id for ICOS stations. `latitude` and `longitude` are
also used to download or extract ERA5(-Land) data. `latitude` is
from -90 to 90 and `longitude` is from -180 to 180. Units of
`altitude` and `reference_height` are meter (m). `reference_height` is
the height where temperature and humidity were measured. There are no
defaults except `altitude = 0` and `reference_height = 2`.

For example:

.. code-block:: python

   [Site]
   name = FR-Hes
   latitude = 48.6742166667
   longitude = 7.06461666667
   altitude = 300.
   reference_height = 22.5


Section [ISBA]
^^^^^^^^^^^^^^

There a three options specific to the land surface model
ISBA. `reference_height_wind` is the height where wind speed was
measured (m). `slope` and `aspect` are the terrain slope and aspect in
degrees. Defaults are that `reference_height_wind` is the same as
`reference_height`, as well as `slope = 0` and `aspect = 0`.

For example:

.. code-block:: python

   [ISBA]
   reference_height_wind = 22.5
   slope =
   aspect =


Section [MuSICA]
^^^^^^^^^^^^^^^^

There a two options specific to the ecosystem model MuSICA. `time2gmt`
is difference between local time and GMT (-12 to 12). The information
will be passed to `ascii2musica`, which adds to the attribute
`time_origin` of the time variable of the netCDF file. If `rsl_yoyo`
is True, then the boundary layer height of ERA5 will be added to the
output file(s). Note that the variable `boundary_layer_height` or
`h_sbl` is not in ERA5-Land but only in ERA5. So forcing files for
runs with the so-called 'yoyo' in MuSICA must use ERA5 data instead of
ERA5-Land data. Defaults are `time2gmt = 0` and `rsl_yoyo = False`.

For example:

.. code-block:: python

   [MuSICA]
   time2gmt = 1
   rsl_yoyo = False


Section [Input]
^^^^^^^^^^^^^^^^

If `input = file` in the section `[Options]`, data will be read from
`inputfile` using information the options `delimiter`, `skiprows`,
`date_columns`, and `date_formats`. `delimiter` is the field separator in the input file. `skiprows` can be set


For example:

.. code-block:: python

   [Input]
   # Input ascii file
   inputfile = /Users/cuntz/data/inrae/hesse/BD_Hesse/DB2/Hesse_DB2_1997.csv
   delimiter = ;
   skiprows = 1
   # number or comma-separated list
   date_columns = 0
   # use space between year (%Y), month (%m), etc. if in separate date columns
   date_format = %d.%m.%Y %H:%M
   # Fraction of time step at which input time is given (default: 0.5):
   #   0   = beginning of time step,
   #   0.5 = middle of time step,
   #   1   = end of time step.
   # Note that MuSICA is using the middle of the time step,
   # and ISBA is using the end of the time step in their forcing files,
   # so output times will be shifted to the middle or end, respectively.
   ftimestep = 1.0
   # Missing values in input data (default: None).
   undef = -9999.




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
