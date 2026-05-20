#!/usr/bin/env python
'''
Download ERA5 or ERA5-Land data from Copernicus Climate Data Store.

It uses the CDS API. You have to have an account on ECMWF's Climate Data Store
(CDS) or Atmospheric Data Store (ADS):
https://cds.climate.copernicus.eu/how-to-api

Be aware that the request is processed once it is queued even if you abort
this script. Queued request can be deleted after login at:
   https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
You can (re-)download the data from that website later as well.


Examples
--------
# Hesse
python get_era5.py -r ERA5-ts -p era/ 48.6742166667,7.06461666667

# Tumbarumba
python get_era5.py -r era5-land -v 10m_u_component_of_wind,10m_v_component_of_wind,2m_temperature,2m_dewpoint_temperature,total_precipitation,surface_pressure,surface_solar_radiation_downwards,surface_thermal_radiation_downwards --area=-35.7833302,148.0166666 -d '2016/2018' -p era


History
-------
   * Written Matthias Cuntz, Jan 2019 - from get_era_interim.py
   * Default area (global) was not working: == instead of =,
     Matthias Cuntz, Dec 2019
   * Added optional reanalysis_model argument to download era5land,
     Stephan Thober, Mar 2020
   * Allow name era5-land, Matthias Cuntz, Jun 2020
   * Return correct file list and not the list of projected filenames for
     chosen area, Matthias Cuntz, Jun 2020
   * Input variable list with MuSICA variables as default,
     Matthias Cuntz, Jun 2020
   * Finalised era5-land capability, Matthias Cuntz, Jun 2020
   * Use numpydoc format, Matthias Cuntz, Jun 2020
   * Bug in single point -> lat1 > lat2, Matthias Cuntz, Feb 2021
   * Added grib format, Matthias Cuntz, Feb 2021
   * Download era5-land in grib format, Matthias Cuntz, Feb 2021
   * No default area: -a must be given, Matthias Cuntz, Sep 2021
   * Respect latency of provision of ERA5 and ERA5-Land for current year,
     Matthias Cuntz & Jerome Ogee, Sep 2022
   * Shift area by 0.01 so account for sites on grid borders of ERA,
     Matthias Cuntz, Oct 2022
   * Allow filenames starting with era5-land and era5land,
     Matthias Cuntz, Oct 2022
   * Half-yearly instead of yearly downloads,
     Klara Bouwen and Remi Lemaire-Patin, 2024
   * era5-land-ts, Matthias Cuntz, Apr 2026
   * area as argument not keyword argument of script, Matthias Cuntz, Apr 2026
   * years -> date, Matthias Cuntz, Apr 2026
   * check time in files, not only on filenames, Matthias Cuntz, Apr 2026
   * era5-ts, Matthias Cuntz, Apr 2026

--------------------------------------------------------

Script was originally adapted from CDS Web API:
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
Further help comes from:
    https://cds.climate.copernicus.eu/how-to-api
    https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+User+Guide
'''
import datetime as dt
import glob
import os
import numpy as np
import pandas as pd
import xarray as xr


__all__ = ['get_era5']


# --------------------------------------------------------------------
# Retrieval function
#

def _get_era5_single_level5(
        variables, date, time, area, target,
        grid=None, reanalysis_model='era5-land-ts',
        output_format='csv'):
    """

    Parameters
    ----------
    variables : list of str
        List of names of variables to retrieve.
    date : str
        Date range as date1/date2 in ISO8601 format.
    time : list of str
        Times to download per day.
    area : string
        Area as 'NorthLat/WestLon/SouthLat/EastLon' if era5 or era5-land,
        and 'lat,lon' if era5-ts or era5-land-ts
    target : string
        Output file name.
    grid : list of float, optional
        Output grid size with east-west (longitude) and
        north-south (latitude) resolution, e.g. [1.0, 1.0]
        Default: [0.25, 0.25] for era5, and [0.1, 0.1] for era5-land
    reanalysis_model : string, optional
        Reanalyis model to download. Can be era5, era5-land / era5land,
        era5-ts / era5-timeseries, or era5-land-ts / era5-land-timeseries.
        Default: 'era5-land-ts'.
    output_format : string, optional
        File format of output file.
        Default: 'grib' if era5 or era5-land, 'csv' if era5-land or
        era5-land-ts.
        Output filenames will be suffixed by '.nc', '.grb' or '.csv',
        respectively.

    Returns
    -------
    target
        output file name

    """
    import cdsapi

    # check reanalysis mode
    retrieve_name = {
        'era5': 'reanalysis-era5-single-levels',
        'era5-ts': 'reanalysis-era5-single-levels-timeseries',
        'era5-timeseries': 'reanalysis-era5-single-levels-timeseries',
        'era5land': 'reanalysis-era5-land',
        'era5-land': 'reanalysis-era5-land',
        'era5-land-ts': 'reanalysis-era5-land-timeseries',
        'era5-land-timeseries': 'reanalysis-era5-land-timeseries'
    }

    if reanalysis_model not in retrieve_name:
        raise ValueError(f'reanalysis_model not known: {reanalysis_model}')

    # Disable: InsecureRequestWarning: Unverified HTTPS request is being made.
    #          Adding certificate verification is strongly advised.
    try:
        import urllib3
        urllib3.disable_warnings()
    except:
        import requests
        from requests.packages.urllib3.exceptions import InsecureRequestWarning
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

    request = {
        'data_format': output_format,  # grib, netcdf, or csv
        'variable': variables,         # ['10m_v_component_of_wind', ...]
        'date': date                   # "2010-01-01", "2010-01-01/2015-12-31"
    }

    if ((reanalysis_model == 'era5') or
        (reanalysis_model == 'era5-land') or
        (reanalysis_model == 'era5land')):
        request.update({'product_type': ['reanalysis'],
                        'download_format': 'unarchived',
                        'time': time,  # "[00:00, 01:00, 02:00, ..., 23:00]
                        'area': area})  # North, West, South, East within
                                        # "90/-180/-90/180" or [90, -180, -90, 180]
    elif ((reanalysis_model == 'era5-ts') or
          (reanalysis_model == 'era5-timeseries') or
          (reanalysis_model == 'era5-land-ts') or
          (reanalysis_model == 'era5-land-timeseries')):
        if isinstance(area, str):
            area = area.split(',')
        request.update({'data_format': 'csv',
                        'location': {'latitude': area[0],
                                     'longitude': area[1]}})

    if grid is not None:
        # 'grid': [1.0, 1.0]
        # Latitude/longitude grid: east-west (longitude)
        # and north-south resolution (latitude).
        # Default: [0.25, 0.25] for era5, [0.1, 0.1] for era5-land
        request.update({'grid': grid})

    print(f'    Request: {request} for target {target}')
    client = cdsapi.Client()
    res = client.retrieve(retrieve_name[reanalysis_model], request)
    ofile = res.download()

    if ((reanalysis_model == 'era5-ts') or
        (reanalysis_model == 'era5-timeseries') or
        (reanalysis_model == 'era5-land-ts') or
        (reanalysis_model == 'era5-land-timeseries')):
        import zipfile
        ropt = {'sep': ',',
                'header': 'infer', 'skiprows': None,
                'index_col': 0, 'parse_dates': True}
        dfl = []
        with zipfile.ZipFile(ofile, 'r') as z:
            for ff in z.namelist():
                with z.open(ff, 'r') as f:
                    dfl.append(pd.read_csv(f, **ropt))
        df = pd.concat(dfl, axis=1)
        df = df.loc[:, ~df.columns.duplicated()]
        df.to_csv(target, sep=',', index=True, date_format='%Y-%m-%d %H:%M:%S',
                  index_label='valid_time')
        os.remove(ofile)
    else:
        if os.path.exists(target):
            os.remove(target)
        os.rename(ofile, target)

    return target


def _get_era5_single_level5_datestarget(
        variables, time, area, datestarget,
        grid=None, reanalysis_model='era5-land-ts',
        output_format='csv'):
    """
    Wrapper function for `_get_era5_single_level5` to parallelise with
    Python's `multiprocessing.Pool`
    
    """
    date, target = datestarget
    return _get_era5_single_level5(
        variables, date, time, area, target, grid,
        reanalysis_model=reanalysis_model,
        output_format=output_format)


# --------------------------------------------------------------------
# Main routine
#

def get_era5(varlist=[], area='90/-180/-90/180', date=None, path='.',
             override=False, reanalysis_model='era5-land-ts',
             output_format=''):
    """
    Download ERA5 or ERA5-Land data from Copernicus Climate Data Store.

    Be aware that the request is processed once it is queued even if you
    abort this function.
    Queued request can be deleted after login at:
       https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
    You can (re-)download the data from that website later as well.

    Parameters
    ----------
    varlist : list of str, optional
        List of names of variables to retrieve.
        Default are input variables for ecosystem models:

        `varlist=['10m_u_component_of_wind',
                  '10m_v_component_of_wind',
                  'total_precipitation',
                  'surface_pressure',
                  '2m_temperature',
                  '2m_dewpoint_temperature',
                  'surface_solar_radiation_downwards',
                  'surface_thermal_radiation_downwards']`

        'boundary_layer_height' is added if era5.

        Look up variable names on:
        https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
    area : string, optional
        Area. Can be one point given as lat,lon
        or a box as NorthLat/WestLon/SouthLat/EastLon.
        Default is the globe: '90/-180/-90/180'.

        Min and Max of Lat and Lon of box must be different by
        at least 0.25 degree for ERA5 and 0.1 degree for ERA5-Land.
    date : str, optional
        Date range as date1/date2 in ISO8601 format.
        if only years are given they default to 01 January to 31 December.
        Default: 1940/today for era5 and era5-ts, and 1950/today for
        era5-land and era5-land-ts.
    path : string, optional
        Output path; wil be created if not existent. Default: '.'
    override : bool, optional
        If True, there will be no check if data is already present in path.
        If False (default), checks files starting with era5_, era5-ts_,
        era5-land_ or era5-land-ts_ in path.
    reanalysis_model : string, optional
        Reanalyis model to download. Can be era5, era5-land / era5land,
        era5-ts / era5-timeseries, or era5-land-ts / era5-land-timeseries.
        Default: 'era5-land-ts'.
    output_format : string, optional
        File format of output file.
        Default: 'grib' if era5 or era5-land, 'csv' if era5-ts or era5-land-ts.
        Output filenames will be suffixed by '.nc', '.grb' or '.csv',
        respectively.

        Model output at ECMWF is stored in grib format. There are limitations
        on the conversion to netCDF using the current ECMWF infrastructure.
        One gets errors like 'One or more variable sizes violate format
        constraints.':
        https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+documentation

        Download in grib format. One can work with grib files using the
        climate data operators (cdo)

            `cdo -t ecmwf sinfov gribfile`
    
        or using xarray in Python with the library cfgrib

            `ds = xr.open_dataset('era5-land_48.71_7.02_48.63_7.1_2026-01-01to2026-04-11.grb', engine='cfgrib')`

    Returns
    -------
    list
        Returns filenames of the output files, either the newly written
        files or the files that contain the requested output:

    Warnings
    --------
    Existing files will only be checked by filename not by content.

    Examples
    --------
    >>> area  = '48/7/47/8'
    >>> date = '1995/2017'
    >>> ofile = get_era5(area=area, date=date, path='.', reanalysis_model='era5')
    >>> file1 = 'era5_48_7_47_8_1995-01-01to2017-12-31.grb')
    >>> if file1 != ofile[0]:
    ...     print(f'Filename {ofile} differs from {file1})

    """
    #
    # Check keywords
    #

    # reanalysis model
    rmodels = ['era5',
               'era5-ts', 'era5-timeseries',
               'era5-land', 'era5land',
               'era5-land-ts', 'era5-land-timeseries']
    rmodel = reanalysis_model.lower()
    assert rmodel in rmodels, (f'Reanalysis model {rmodel} must be in'
                               f' {rmodels}')

    # defaults for reanalysis models
    if rmodel == 'era5':
        resolution = 0.25
        mindate = '1940-01-01'
        # era5 needs days that are already present in the product.
        # It has a latency of about 5 days.
        latency = 10
        if not output_format:
            output_format = 'grib'
        omodel = 'era5'
    elif (rmodel == 'era5-ts') or (rmodel == 'era5-timeseries'):
        resolution = 0.25
        mindate = '1940-01-01'
        # era5 needs days that are already present in the product.
        # It has a latency of about 5 days.
        latency = 10
        if not output_format:
            output_format = 'csv'
        omodel = 'era5-ts'
    elif (rmodel == 'era5-land') or (rmodel == 'era5land'):
        resolution = 0.1
        mindate = '1950-01-02'
        # era5-land has a latency of about 2-3 months.
        # But it seems to accept days that do not exist in the product yet.
        latency = 0
        if not output_format:
            output_format = 'grib'
        omodel = 'era5-land'
    elif (rmodel == 'era5-land-ts') or (rmodel == 'era5-land-timeseries'):
        resolution = 0.1
        mindate = '1950-01-02'
        # era5-land has a latency of about 2-3 months.
        # But it seems to accept days that do not exist in the product yet.
        latency = 0
        if not output_format:
            output_format = 'csv'
        omodel = 'era5-land-ts'
    else:
        raise ValueError(f'ERA5 model unknown: {reanalysis_model}')
    output_format = output_format.lower()

    # output format
    if output_format == 'netcdf':
        suffix = 'nc'
        maxitem = 60000  # since April 2024
    elif output_format == 'grib':
        suffix = 'grb'
        maxitem = 100000
    elif output_format == 'csv':
        suffix = 'csv'
        maxitem = np.iinfo(int).max  # practically unlimited
    else:
        raise ValueError(f'Output format unknown: {output_format}')

    # variables
    if len(varlist) == 0:
        varlist = ['10m_u_component_of_wind', '10m_v_component_of_wind',
                   '2m_temperature', '2m_dewpoint_temperature',
                   'total_precipitation', 'surface_pressure',
                   'surface_solar_radiation_downwards',
                   'surface_thermal_radiation_downwards']
        if (rmodel == 'era5') or (rmodel == 'era5-ts'):
            varlist.append('boundary_layer_height')

    # area
    fstr = 'area format is lat,lon or NorthLat/WestLon/SouthLat/EastLon.'
    estr = 'area must be in 90/-180/-90/180'
    if '/' in area:
        sarea = [ round(float(ll) / resolution) * resolution
                  for ll in area.split('/') ]
        assert len(sarea) == 4, fstr
        assert sarea[0] >= sarea[2] + resolution, fstr
        assert sarea[1] + resolution <= sarea[3], fstr
        assert ((sarea[0] <= 90.) and (sarea[1] >= -180.) and
                (sarea[2] >= -90.) and (sarea[3] <= 180.)), estr
        area = f'{sarea[0]:.2f}/{sarea[1]:.2f}/{sarea[2]:.2f}/{sarea[3]:.2f}'
    else:
        assert ',' in area, fstr
        sarea = [ round(float(ll) / resolution) * resolution
                  for ll in area.split(',') ]
        assert len(sarea) == 2, fstr
        lat, lon = sarea
        if rmodel in ['era5-ts', 'era5-timeseries',
                      'era5-land-ts', 'era5-land-timeseries']:
            area = f'{lat:.2f},{lon:.2f}'
        else:
            # Single points do not work in era5 and era5-land but in
            # era-land-ts. So we need to encompass an actual grid
            # point. Take slightly more to the left and bottom to
            # take lower value when exactly on grid border.
            n = lat + resolution / 2 - 0.01
            w = lon - resolution / 2 + 0.01
            s = lat - resolution / 2 + 0.01
            e = lon + resolution / 2 - 0.01
            area = f'{n:.2f}/{w:.2f}/{s:.2f}/{e:.2f}'

    # date
    if date:
        if '/' in date:
            date1, date2 = date.split('/')
            try:
                yr = int(date1)
                startdate = f'{date1}-01-01'
            except ValueError:
                startdate = date1
            try:
                yr = int(date2)
                enddate = f'{date2}-12-31'
            except ValueError:
                enddate = date2
        else:
            try:
                yr = int(date)
                startdate = f'{date}-01-01'
                enddate = f'{date}-12-31'
            except ValueError:
                raise ValueError(
                    f'date range as date1/date2 in ISO8601 format.'
                    f' If no range (i.e. no /), date can only be a'
                    f' single year. Given: {date}')
    else:
        startdate = mindate
        today = dt.date.today()
        edate = today - np.timedelta64(latency, 'D')
        enddate = edate.strftime('%Y-%m-%d')
    mdate = pd.to_datetime(mindate, format='ISO8601')
    sdate = pd.to_datetime(startdate, format='ISO8601')
    edate = pd.to_datetime(enddate, format='ISO8601')
    assert sdate >= mdate, (
        f'Start date {startdate} must be greater {mindate}'
        f' for reanalysis model {reanalysis_model}.')
    today = dt.date.today()
    mxdate = today - np.timedelta64(latency, 'D')
    mxdate = pd.to_datetime(mxdate.strftime('%Y-%m-%d'), format='ISO8601')
    if edate > mxdate:
        edate = mxdate

    # output directory
    if not os.path.exists(path):
        os.makedirs(path)

    #
    # Check existing files
    #

    if output_format == 'netcdf':
        files = glob.glob(path + '/*.nc*')
    elif output_format == 'grib':
        files = glob.glob(path + '/*.gr*b*')
    elif output_format == 'csv':
        files = glob.glob(path + '/*.csv')
    ifiles = []
    for ff in files:
        base = os.path.basename(ff)
        if (base.startswith(f'{omodel}_') and (not base.endswith('.idx'))):
            ifiles.append(ff)
    files = ifiles

    # https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings
    vardict = {'10m_u_component_of_wind': 'u10',
               '10m_v_component_of_wind': 'v10',
               '2m_temperature' : 't2m',
               '2m_dewpoint_temperature': 'd2m',
               'total_precipitation': 'tp',
               'surface_pressure': 'sp',
               'surface_solar_radiation_downwards': 'ssrd',
               'surface_thermal_radiation_downwards': 'strd',
               'boundary_layer_height': 'blh'}
    if len(files) > 0:
        dvars = varlist.copy()
        dsdate = sdate
        dedate = edate
        ihavelat = True
        if (output_format == 'netcdf') or (output_format == 'grib'):
            ds = xr.open_mfdataset(files)
            # check for missing variables
            if 'u10m' in ds.variables:
                vardict.update({'10m_u_component_of_wind': 'u10m'})
            if 'v10m' in ds.variables:
                vardict.update({'10m_v_component_of_wind': 'v10m'})
            for dd in varlist:
                if vardict[dd] in ds.variables:
                    dvars.remove(dd)
            # check for missing time steps
            if 'time' in ds.dims:
                vtime = 'time'
            elif 'valid_time' in ds.dims:
                vtime = 'valid_time'
            else:
                ds.close()
                raise ValueError("'time' dimension not found")
            tmin = ds[vtime].values.min()
            tmax = ds[vtime].values.max()
            if dedate > tmax:
                if (dsdate <= tmax) and (dsdate >= tmin):
                    dsdate = tmax + np.timedelta64(1, 'h')
            elif dsdate < tmin:
                if (dedate <= tmax) and (dedate >= tmin):
                    dedate = tmin - np.timedelta64(1, 'h')
            else:
                dsdate = None
                dedate = None
            # check for missing lat,lon
            if 'lat' in ds.dims:
                vlat = 'lat'
            elif 'latitude' in ds.dims:
                vlat = 'latitude'
            else:
                ds.close()
                raise ValueError("'latitude' dimension not found")
            if 'lon' in ds.dims:
                vlon = 'lon'
            elif 'longitude' in ds.dims:
                vlon = 'longitude'
            else:
                ds.close()
                raise ValueError("'longitude' dimension not found")
            elat = np.round(ds[vlat].values, 2)
            elon = np.round(ds[vlon].values, 2)
            if ',' in area:
                lat, lon = [ float(ll) for ll in area.split(',') ]
                if lat not in elat:
                    ihavelat = False
                if lon not in elon:
                    ihavelat = False
            else:
                n, w, s, e = [ float(ll) for ll in area.split('/') ]
                res2 = resolution / 2
                if ((n >= elat.max() + res2) or
                    (w <= elon.min() - res2) or
                    (s <= elat.min() - res2) or
                    (e >= elon.max() + res2)):
                    ihavelat = False
            ds.close()
        else:
            ropt = {'sep': ',',
                    'header': 'infer',
                    'skiprows': None,
                    'index_col': 0,
                    'parse_dates': True}
            dfl = []
            for ff in files:
                dfl.append(pd.read_csv(ff, **ropt))
            # try to combine files with different variables and
            # files with different time steps
            df = dfl[0]
            for df1 in dfl[1:]:
                df, df1 = df.align(df1, join="outer")
                df1 = df1.fillna(df)
                df = df.fillna(df1)
            # check for missing variable
            if 'u10m' in df.columns:
                vardict.update({'10m_u_component_of_wind': 'u10m'})
            if 'v10m' in df.columns:
                vardict.update({'10m_v_component_of_wind': 'v10m'})
            for dd in varlist:
                if vardict[dd] in df.columns:
                    dvars.remove(dd)
            # check for missing time steps
            tmin = df.index.min()
            tmax = df.index.max()
            if dedate > tmax:
                if (dsdate <= tmax) and (dsdate >= tmin):
                    dsdate = tmax + np.timedelta64(1, 'h')
            elif dsdate < tmin:
                if (dedate <= tmax) and (dedate >= tmin):
                    dedate = tmin - np.timedelta64(1, 'h')
            else:
                dsdate = None
                dedate = None
            # check for missing lat,lon
            if 'lat' in df.columns:
                vlat = 'lat'
            elif 'latitude' in df.columns:
                vlat = 'latitude'
            else:
                raise ValueError("'latitude' column not found")
            if 'lon' in df.columns:
                vlon = 'lon'
            elif 'longitude' in df.columns:
                vlon = 'longitude'
            else:
                raise ValueError("'longitude' column not found")
            elat = np.round(df[vlat].values, 2)
            elon = np.round(df[vlon].values, 2)
            if ',' in area:
                lat, lon = [ float(ll) for ll in area.split(',') ]
                if lat not in elat:
                    ihavelat = False
                if lon not in elon:
                    ihavelat = False
            else:
                n, w, s, e = [ float(ll) for ll in area.split('/') ]
                res2 = resolution / 2
                if ((n >= elat.max() + res2) or
                    (w <= elon.min() - res2) or
                    (s <= elat.min() - res2) or
                    (e >= elon.max() + res2)):
                    ihavelat = False

        # Return if everything already present
        if ((len(dvars) == 0) and (dsdate is None) and ihavelat):
            return files
        else:
            if (dsdate is None):
                # download missing variables in time period of existing files
                sdate = tmin
                edate = tmax
            elif len(dvars) == 0:
                # download missing time steps for all variables
                sdate = dsdate
                edate = dedate
            else:
                # download all time steps for all variables
                pass

    #
    # Prepare download
    #

    # CDS allows to download 100000 items at a time.
    # Since April 2024, it allows only 60000 items if output format
    # is netcdf.
    nhours = int((edate - sdate + np.timedelta64(1, 'D'))
                 / np.timedelta64(1, 'h'))
    ndays = nhours // 24
    nvars = len(varlist)
    nitems = nhours * nvars

    if nitems > maxitem:
        idays = maxitem // (nvars * 24)
        iitems = idays * nvars * 24
    else:
        idays = ndays
        iitems = nitems

    nfiles = nitems // iitems
    if (nitems % iitems) > 0:
        nfiles += 1

    times = ['00:00', '01:00', '02:00',
             '03:00', '04:00', '05:00',
             '06:00', '07:00', '08:00',
             '09:00', '10:00', '11:00',
             '12:00', '13:00', '14:00',
             '15:00', '16:00', '17:00',
             '18:00', '19:00', '20:00',
             '21:00', '22:00', '23:00']

    isoform = '%Y-%m-%d'
    d1 = np.timedelta64(1, 'D')
    ddate = np.timedelta64(idays, 'D')

    dates = []
    targets = []
    iedate = pd.to_datetime(sdate) - d1
    for nn in range(nfiles):
        isdate = iedate + d1
        iedate = isdate + ddate
        if iedate > edate:
            iedate = pd.to_datetime(edate)
        idate = isdate.strftime(isoform) + '/' + iedate.strftime(isoform)
        dates.append(idate)
        
        # output filename
        if ',' in area:
            oarea = area.replace(',', '_')
        else:
            oarea = area.replace('/', '_')
        odate = idate.replace('/', 'to')
        target = f'{path}/{omodel}_{oarea}_{odate}.{suffix}'
        targets.append(target)

    #
    # Download
    #

    if len(dates) > 1:
        print('Pool')
        print('    dates', dates)
        print('    targets', targets)
        from functools import partial
        getit = partial(
            _get_era5_single_level5_datestarget, varlist, times, area,
            reanalysis_model=rmodel, output_format=output_format)
        from multiprocessing import Pool
        pool = Pool(processes=len(dates))
        pool.map(getit, zip(dates, targets), None)
    elif len(dates) == 1:
        print('Retrieve ', dates[0], targets[0])
        _ = _get_era5_single_level5(
            varlist, dates[0], times, area, targets[0],
            reanalysis_model=rmodel, output_format=output_format)
    else:
        raise ValueError('Dommage ...')

    return targets


# --------------------------------------------------------------------
# Script
#

if __name__ == "__main__":

    import argparse

    date = None
    output_format = ''
    override = False
    path = '.'
    reanalysis_model = 'era5-land-ts'
    varstr = ''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            'Download ERA5 or ERA5-Land data from Copernicus Climate'
            ' Data Store https://climate.copernicus.eu/climate-data-store.'))
    parser.add_argument('-d', '--date', action='store',
                        default=date, dest='date', metavar='date1/date2',
                        help=('date range as date1/date2 in ISO8601 format.'
                              ' Only years default to 01 January to'
                              ' 31 December. (default:'
                              ' 1940/today for era5 and era5-ts,'
                              ' 1950/today for era5-land and era5-land-ts.)'))
    parser.add_argument('-f', '--format', action='store',
                        default=output_format, dest='output_format',
                        metavar='format',
                        help=('Output format netcdf, grib, or csv.'
                              ' (default: grib if era5 or era5-land,'
                              ' csv if era5-ts or era5-land-ts).'))
    parser.add_argument('-o', '--override', action='store_true',
                        default=override, dest='override',
                        help=('Do not check that output file already exists'
                              ' that includes request. Override existing'
                              ' output file (default: False).'))
    parser.add_argument('-p', '--path', action='store',
                        default=path, dest='path', metavar='path',
                        help=('Output directory (default: current directory'
                              ' ".").'))
    parser.add_argument('-r', '--reanalyis-model', action='store',
                        default=reanalysis_model, dest='reanalyis_model',
                        metavar='reanalysis_model',
                        help=('Reanalyis model to download, either: era5,'
                              ' era5-ts / era5-timeseries,'
                              ' era5-land / era5land,'
                              ' era5-land-ts / era5-land-timeseries'
                              ' (default: era5-land-ts)'))
    parser.add_argument('-v', '--variables', action='store',
                        default=varstr, dest='varstr', metavar='variables',
                        help=('Comma-separated list of variables var1,var2,...'
                              ' default: shortwave and longwave radiation,'
                              ' wind speed, 2m temperature and dew point,'
                              ' air pressure, and total precipitation.'
                              ' Boundary layer height will be added if era5.'))
    parser.add_argument('area', nargs='?',
                        default=None, metavar='area',
                        help=('area format either as lat,lon or as'
                              ' NorthLat/WestLon/SouthLat/EastLon, e.g.'
                              ' global 90/-180/-90/180, mandatory.'))

    args = parser.parse_args()
    date = args.date
    output_format = args.output_format
    override = args.override
    path = args.path
    reanalysis_model = args.reanalyis_model
    varstr = args.varstr
    area = args.area

    del parser, args

    if area is None:
        raise ValueError(
            'area as either lat,lon or NorthLat/WestLon/SouthLat/EastLon'
            ' must be given, e.g. "48.67,7.06" or "90/-180/-90/180".')

    # variables str to list
    if ',' in varstr:
        varlist = varstr.split(',')
    elif varstr == '':
        varlist = []
    else:
        # single variable given
        varlist = [varstr]

    import time as ptime
    t1 = ptime.time()

    era5files = get_era5(
        varlist=varlist, area=area, date=date, path=path, override=override,
        reanalysis_model=reanalysis_model, output_format=output_format)
    print('ERA5 files: ', era5files)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.) if (t2 - t1) > 60. else
             '[s]: {:d}'.format(int(t2 - t1)))
    print('Time elapsed', strin)
