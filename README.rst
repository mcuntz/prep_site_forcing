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
``ascii2netcdf``, e.g. ``ascii2musica``, that produces forcing data for
the ecosystem model in the very model-specific form. Currently there
are the two model-specific scripts ``ascii2isba.py`` and
``ascii2musica.py`` extending the base class ``ascii2Netcdf`` from
``class_ascii2netcdf.py``.


Installation
------------

There is no installation. Simply clone the repository:

.. code-block:: bash

   git clone https://github.com/mcuntz/prep_site_forcing.git

and run the script:

.. code-block:: bash

   python prep_site_forcing.py FR-Hes.cfg

You can check the results of the linear bias-correction of the
ERA5(-Land) data by plotting them into a file:

.. code-block:: bash

   python prep_site_forcing.py -p bias_correction.pdf -t pdf FR-Hes.cfg

See ``python prep_site_forcing.py -h`` for help.

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

   * cdsapi_  if ERA5 must be downloaded; needs a Copernicus account (see Section [ERA5])
   * `icoscp_core`_  if ICOS data is used (see Section [ICOS])


Config file
-----------

The driver script ``prep_site_forcing.py`` is controlled by a
configuration file, which is in the simple Python configparser
format. It has sections with case-sensitive names in brackets,
e.g. ``[Site]``, and case-insensitive options within each
section. Mandatory sections are ``[Model]``, ``[Site]``, and
``[VarNames]``. Further, the option ``input`` in ``[Options]`` must be set.

See ``FR-Hes.cfg`` for an example that can be used as a template for
other sites. The example files is highly commented and should be
rather self-explaining. Or see the explanations below.


Section [Model]
^^^^^^^^^^^^^^^

The ``model`` name must be given. This is mostly used to call the
correct ``ascii2netcdf`` script.

For example:

.. code-block:: python

   [Model]
   model = MuSICA


Section [Options]
^^^^^^^^^^^^^^^^^

There are several general options.

``input`` tells, which data should be taken. It can be ``file``, ``ICOS``,
or ``ERA5``. If it is ``file``, the section ``[Input]`` will be used and
data will be read from a local file. If ``input`` is set to ``ICOS``, the
section ``[ICOS]`` will be used and data is downloaded from the ICOS
Carbon Portal or read from a local file in ICOS convention. If ``input``
is ``ERA5``, the section ``ERA5`` will be used and data downloaded from
the Copernicus archive or read from local files will be used directly
as forcing data. There is no default value.

A Mean Absolute Deviation (MAD) filter can be applied to the input
data to catch outliers. Set ``mad_z`` greater than zero to apply an MAD
filter with ``z = mad_z``. Default is not to apply an MAD filter (``mad_z
= 0``).

The imputation method can be chosen. Currently there is linear
interpolation (``imputation_method = 0``) to fill gaps or using
ERA5(-Land) data that gets bias-corrected with the existing local data
(``imputation_method = 1``; Vuichard and Papale, ESSD 2013,
https://doi.org/10.5194/essd-7-157-2015). Default is
``imputation_method = 1``.

``make_netcdf`` controls if a function ``ascii2netcdf`` from a file
``ascii2netcdf.py`` will be called. 'netcdf' is thereby replaced by the
lowercase model name. If ``make_netcdf = False``, the script stops after
production of the csv file (this implies that the next option is True,
i.e. ``keep_csv = True``). Default is ``make_netcdf = True``.

If ``keep_csv`` is True, then the csv file will not be deleted after
calling ``ascii2netcdf``. Default is ``keep_csv = False``.

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

Data in the ``[Site]`` section is mostly used to write it into the
netcdf files.  The site ``name`` is basically for information but it is
also used as site id for ICOS stations. ``latitude`` and ``longitude`` are
also used to download or extract ERA5(-Land) data. ``latitude`` is
from -90 to 90 and ``longitude`` is from -180 to 180. Units of
``altitude`` and ``reference_height`` are meter (m). ``reference_height`` is
the height where temperature and humidity were measured. There are no
defaults except ``altitude = 0`` and ``reference_height = 2``.

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
ISBA. ``reference_height_wind`` is the height where wind speed was
measured (m). ``slope`` and ``aspect`` are the terrain slope and aspect in
degrees. Defaults are that ``reference_height_wind`` is the same as
``reference_height``, as well as ``slope = 0`` and ``aspect = 0``.

For example:

.. code-block:: python

   [ISBA]
   reference_height_wind = 22.5
   slope =
   aspect =


Section [MuSICA]
^^^^^^^^^^^^^^^^

There a two options specific to the ecosystem model MuSICA. ``time2gmt``
is difference between local time and GMT (-12 to 12). The information
will be passed to ``ascii2musica``, which adds to the attribute
``time_origin`` of the time variable of the netCDF file. If ``rsl_yoyo``
is True, then the boundary layer height of ERA5 will be added to the
output file(s). Note that the variable ``boundary_layer_height`` or
``h_sbl`` is not in ERA5-Land but only in ERA5. So forcing files for
runs with the so-called 'yoyo' in MuSICA must use ERA5 data instead of
ERA5-Land data. Defaults are ``time2gmt = 0`` and ``rsl_yoyo = False``.

For example:

.. code-block:: python

   [MuSICA]
   time2gmt = 1
   rsl_yoyo = False


Section [Input]
^^^^^^^^^^^^^^^^

If ``input = file`` in the section ``[Options]``, data will be read from
``inputfile`` using `pandas.read_csv`_. The options ``sep``, ``header``,
``index_col``, ``usecols``, ``skiprows``, ``na_values``, ``parse_dates``, and
``date_format`` are implemented, which allows to read a vast majority of
data files. Please see the documentation of `pandas.read_csv`_ for
details.

``ftimestep`` indicates which time point is represented by the time
information in the input file. ``0`` means that the times represent the
beginning of the time steps, ``0.5`` the middle, and ``1`` the end of time
steps. MuSICA, for example, is using the middle of the time step and
ISBA is using the end of the time step in their forcing files. Time
steps will hence be shifted appropriately in the ``ascii2netcdf``
routines (not yet in the csv file).

For example:

.. code-block:: python

   [Input]
   inputfile = /Users/cuntz/data/inrae/hesse/BD_Hesse/DB2/Hesse_DB2_1997.csv
   sep = ;
   header = 0
   index_col = 0
   usecols =
   skiprows = 1
   na_values = -9999.
   parse_dates = True
   date_format = %d.%m.%Y %H:%M
   ftimestep = 1.0


Section [Output]
^^^^^^^^^^^^^^^^

Section ``[Output]`` gives the netCDF output filename. Default is the
input file name with the suffix replaced by .nc if ``input = file``, or
``prep_site_forcing.nc`` otherwise. The name of the csv file is the name
of the output file with suffix replaced by .csv. ``fill_value`` is the
missing value used in the netCDF file, which are highly model
specific. For example, MuSICA is using the netCDF default ``_FillValue``
(``fill_value =``), while ISBA is using ``fill_value = -9999999``.

``startdate`` and ``enddate`` (both inclusive) can be given in ISO8601 format (e.g. 1994-12-31 23:30) to restrict the forcing file between the two dates. Defaults are the first and last dates of the input file. ``startdate``, ``enddate``, and ``timestep`` have to be given if ``input`` is not ``file``, i.e. either ICOS or ERA5(-Land). In this cases, ``timestep`` will be the time step of the output file (by linear interpolation from hourly ICOS or ERA5(-Land) data). Notation such as 1800s or 30min can be used (`pandas timeseries`_).

For example:

.. code-block:: python

   [Output]
   outputfile = ISBA_in_FR-Hes_1997-era5-3.nc
   fill_value = -9999999.
   startdate = 1997-05-01 00:00
   enddate = 1997-05-31 23:30
   timestep = 1800s


Section [ICOS]
^^^^^^^^^^^^^^

Input can also come directly from the ICOS Carbon Portal (``input = ICOS``) using the ``[Site].name`` as station id. One has to have the Python library `icoscp_core`_ installed and have it initialised with:

.. code-block:: python

   from icoscp_core.icos import auth
   auth.init_config_file()

Known ``icos_product`` are L2, NRT, and Fluxnet. For ``icos_meteo``, one
can use the original Meteosens stream for L2 and NRT, which gives the
individual meteorological variables such as SW_IN_1_1_1 and
SW_IN_1_1_2. Or one uses the Meteo stream, which are all redundant
variables aggregated, such as SW_IN. Or one uses the Fluxnet meteo
stream, which gives gap-filled aggregated variables such as SW_IN_F.
``icos_product = Fluxnet`` includes not only ICOS data but also pre-ICOS
data for the site in the same gap-filled Fluxnet meteo
format. ``icos_meteo`` is ignored in this case.

``icos_product`` can also name of a local file in a appropriate standard
format, i.e. it can be read with ``df = pd.read_csv(file, index_col=0,
parse_dates=True, date_format='ISO8601', na_values='-9999')``. In this case, the units must be part of the header row, e.g. SW_IN (W m-2). For example:

.. code-block:: bash

   TIMESTAMP_END,CO2 (µmol mol-1),...
   2025-12-31 23:00:00,432.47,...

Such a file can be produced, for example, using the function
``write_icos`` in ``icos.py``.

One can use the quality flag for aggregated variables. ``icos_qc = 2``
uses all data, and ``icos_qc = 0`` uses only the original, measured
data. Filling of data gaps will then be done with ERA5(-Land) data
(see below).

For example:

.. code-block:: python

   [ICOS]
   icos_product = L2
   icos_meteo = Meteosens
   icos_qc = 2


Section [ERA5]
^^^^^^^^^^^^^^

Missing data can be filled with ERA5 reanalysis data. There are the
products ERA5 and ERA5-land, which are on different resolutions (0.25
vs 0.1 degree). They are also stored in different formats, grib and
zarr. The latter is optimised for remote access. This is indicated by
'-ts' in the ``era5type`` name. It should generally be ``era5type =
era5-land-ts``. ERA5-Land has no boundary layer height. So if
``boundary_layer_height`` (or ``h_sbl``) is needed, such as in the case of
``rsl_yoyo = True``, then ``era5type = era5-ts`` is preferred. The script
will warn if ``era5type`` is ``era5-land-ts`` and set it to ``era5-ts``.

The script calls the function ``get_era5`` from the file
``get_era5.py``. One needs the Python library ``cdsapi`` and a Copernicus
account (see how-to-api_). The check the files in the directory
``era5path`` depending on the output format wanted, which is netcdf if
era5 or era5-land (``*.nc*``) and csv if era5-ts or era5-land-ts
(``*.csv``). It check for available variables and timesteps and
downloads the missing information.

If ``input = ERA5`` in ``[Options]``, then the ERA data will directly used
as forcing data.

For example:

.. code-block:: python

   [ERA5]
   era5path = era5
   era5type = era5-land-ts


Section [CO2]
^^^^^^^^^^^^^

Most ecosystem models need atmospheric CO2 concentrations. Missing
data will be filled using a simple input file ``co2file``. The CO2 data
provided (``cmip6_co2_hist-ssp370_1850-2100.csv``) are annual mean
values from 1850 to 2100 from the CMIP6 SSP3.7 scenario. They are
delimited by comma (``co2delimiter = ,``) with a decimal date in the
first column (``co2date_column = 0``) and CO2 mixing ratios in ppm in
the second column (``co2co2_column = 1``). More sophisticated filling
should replace the method ``fill_co2``.

For example:

.. code-block:: python

   [CO2]
   co2file = cmip6_co2_hist-ssp370_1850-2100.csv
   co2delimiter = ,
   co2date_column = 0
   co2co2_column = 1


Section [VarNames]
^^^^^^^^^^^^^^^^^^

Standard forcing variables for ecosystems models are shortwave
incoming radiation (swdown), longwave incoming radiation (lwdown),
atmospheric pressure (psurf), air temperature (tair) and humidity
(qair), wind speed (wind_speed), precipitation (precip), as well as
atmospheric CO2 concentrations (co2). In addition, wind direction
(wind_dir) can be used, and MuSICA also needs the atmospheric boundary
layer height for some applications. Precipitation can be separated in
liquid (rainf) and solid (snowf) precipitation.

So one has to give the names in the input file or data stream that
corresponds to the above variables. One has to give namely
``name_swdown``, ``name_lwdown``, ``name_psurf``, ``name_tair``, ``name_qair``,
``name_wind_speed``, ``name_co2air``, and optionally ``name_wind_dir``. One
has to give further either ``name_precip``, or both ``name_rainf`` and
``name_snowf``. Humidity (qair) can be relative humidity, specific
humidity, or vapour pressure deficit (VPD). The variable will be
identified bu its unit (see below).

The names can be regular expression such as ``TA_.*_1_1``. Columns will
be filtered, which uses ``re.search(name,
available_variables)``. Variables will be averaged over all columns
found. The filtering with ``re.search``_ implies that the name ``Var_1``,
for example, also finds columns named ``Var_2/Var_1``, ``Var_1_QC``, or
similar. In this case, one can start the variable name with ``^``,
i.e. ``^Var_1`` in this case, or end it with ``$``, i.e. ``Var_1$``, for
example.

The ICOS Fluxnet product has for each variable also quality control
columns, e.g. ``TA_F`` and ``TA_F_QC``. So one would end the variables
with ``$`` and the section ``[VarNames]`` would look like:

.. code-block:: python

   [VarNames]
   name_co2air =
   name_lwdown = LW_IN_F$
   name_psurf = PA_F$
   name_qair = VPD_F$
   name_swdown = SW_IN_F$
   name_tair = TA_F$
   name_wind_dir =
   name_wind_speed = WS_F$
   name_precip = P_F$

The ICOS Meteo products have aggregated variables, e.g. ``TA`` as an
average over all air temperature measurements on top of the tower. But
they also have a temperature profile with the variables ``TA_1``,
``TA_2``, ... So one would also end the variables with ``$`` and the
section ``[VarNames]`` would look like:

.. code-block:: python

   [VarNames]
   name_co2air = CO2$
   name_lwdown = LW_IN$
   name_psurf = PA$
   name_qair = RH$
   name_swdown = SW_IN$
   name_tair = TA$
   name_wind_dir = WD$
   name_wind_speed = WS$
   name_precip = P$

The ICOS Meteosens products have the individual sensors such as
``TA_1_1_1``. One can decide to average over all or some of the sensors
with the regular expressions, or rather use simply one sensor and use
the redundant measurements only when the primary sensor has missing
data (see section [AlternativeVarNames] below). For example, take
single columns for most variables, average all longwave radiation and
air temperature sensors, but average only sensors 1 and 2 for
precipitation:

.. code-block:: python

   [VarNames]
   name_co2air = CO2
   name_lwdown = LW_IN_1_._1
   name_psurf = PA_1_1_1
   name_qair = RH_1_1_1
   name_swdown = SW_IN_1_1_1
   name_tair = TA_.*_1_1
   name_wind_dir = WD_1_1_1
   name_wind_speed = WS_1_1_1
   name_precip = P_[1-2]_1_1

The files for the historical data at FR-Hes used individual names and
had also calculated columns such as the ratio between one PAR sensor
and one global radiation sensor ``PAR_JY/Rg_Kipp H1``, so that shortwave
incoming radiation would be ``^Rg_Kipp H1``:

.. code-block:: python

   [VarNames]
   name_co2air = CO2_EC H1
   name_lwdown = IRLinc H1
   name_psurf = Patm
   name_qair = H2O_Vais H1
   name_swdown = ^Rg_Kipp H1
   name_tair = T_Psy H1
   name_wind_dir = WD H1
   name_wind_speed = WS_EC H1
   name_precip = Prec H1

MuSICA might need more forcing variables such as isotopic forcing
data. Extra variables ``extra_vars`` can hence be extracted from the
input file and data stream and written into the output file with the
names ``extra_names``. Note that these variables will not be imputed. Adding soil moisture and soil heat flux from the historical FR-Hes data to the forcing file would be:

.. code-block:: python

   [VarNames]
   name_co2air = CO2_EC H1
   name_lwdown = IRLinc H1
   name_psurf = Patm
   name_qair = H2O_Vais H1
   name_swdown = ^Rg_Kipp H1
   name_tair = T_Psy H1
   name_wind_dir = WD H1
   name_wind_speed = WS_EC H1
   name_precip = Prec H1
   extra_vars = SWC-10_. H1, SWC-30_. H1, SWC-55_. H1, SWC-80_. H1, SWC-120_. H1, Gsol_. H1
   extra_names = SWC_10, SWC_30, SWC_55, SWC_80, SWC_120, G


Section [AlternativeVarNames]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of averaging variables such as primary and backup air
temperature ``TA_._1_1``, one could also use one sensor and only fill in
data from another sensor in case of missing data. This 'backup'
variables can be given in the section ``[AlternativeVarNames]``. The
alternative variables have to have the same units as the primary
variables. Using only the main sensors of the ICOS meteosens and only
filling in the backup sensors in case of missing data would give:

.. code-block:: python

   [VarNames]
   name_co2air = CO2
   name_lwdown = LW_IN_1_1_1
   name_psurf = PA_1_1_1
   name_qair = RH_1_1_1
   name_swdown = SW_IN_1_1_1
   name_tair = TA_.*_1_1
   name_wind_dir = WD_1_1_1
   name_wind_speed = WS_1_1_1
   name_precip = P_1_1_1

   [AlternativeVarNames]
   aname_co2air =
   aname_lwdown = LW_IN_1_2_1
   aname_psurf = PA_1_2_1
   aname_qair = RH_2_1_1
   aname_swdown = SW_IN_1_1_2
   aname_tair = TA_2_1_1
   aname_wind_dir =
   aname_wind_speed = WS_1_2_1
   aname_precip = P_2_1_1


Section [VarUnits]
^^^^^^^^^^^^^^^^^^

The script has to know the units of the variables, which are given in
the section ``[VarUnits]``. Known units are:
   * ['W/m2', 'W m-2'] for shortave and longwave radiation
   * ['C', 'degreeC', 'degree C', 'degC', 'deg C', '°C'] for air
     temperature (otherwise Kelvin assumed)
   * ['mol/mol', 'mol mol-1', 'mmol/mol', 'mmol mol-1'], ['kPa',
     'hPa', 'Pa'] and ['%', 'percent', '', '0-1'] for specific
     humidity, VPD, and relative humidity, respectively
   * ['hPa', 'kPa', 'mbar', 'bar'] for pressure
   * ['mm', 'mm/dt', 'mm dt-1','kg/m2/dt', 'kg m-2 dt-1'] for
     precipitation
   * ['m/s', 'm s-1', ''] for wind speed
   * ['degree', 'deg', '°', ''] for wind direction
   * ['ppm', 'ppmv', 'µmol/mol', 'µmol mol-1'] for CO2
   * ['m', ''] for boundary layer height

For example:

.. code-block:: python

   [VarUnits]
   unit_co2air = ppm
   unit_lwdown = W/m2
   unit_psurf = hPa
   unit_qair = mmol/mol
   unit_swdown = W/m2
   unit_tair = degC
   unit_wind_dir = degree
   unit_wind_speed = m/s
   unit_precip = mm
   unit_rainf = mm
   unit_snowf = mm

``[VarUnits]`` is ignored in case ``input`` is ``ICOS`` or ``ERA5``. If the
``icos_product`` is a local file, then the units must be in the header
row (see section ``[ICOS]``).


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
.. _icoscp_core: https://github.com/ICOS-Carbon-Portal/data/tree/master/src/main/python/icoscp_core
.. _pandas.read_csv: https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
.. _pandas timeseries: https://pandas.pydata.org/docs/user_guide/timeseries.html#timeseries-offset-aliases
.. _how-to-api: https://cds.climate.copernicus.eu/how-to-api
.. _re.search: https://docs.python.org/3/library/re.html#re.search
