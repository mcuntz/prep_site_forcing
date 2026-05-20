#!/usr/bin/env python
"""
Write MuSICA netcdf forcing file from csv input file.

# Script can be called on the command line:
python ascii2musica.py -s 2005-01-01T00:30 -e 2006-01-01T00:30 -o MuSICA_site_in_FR-Hes_2005-01.nc -t 1.0 --site_latitude=48.6742166667 --site_longitude=7.06461666667 --time2gmt=0 --forcing_height=22.5 -f=-9999 MuSICA_site_in_FR-Hes_2005.csv

# The wrapper function can be called from within Python:
from ascii2musica import ascii2musica
ascii2musica(csvfile,
             namelist=namelist, site_latitude=site_latitude,
             site_longitude=site_longitude, site_altitude=site_altitude,
             forcing_height=forcing_height, time2gmt=time2gmt,
             startdate=startdate, enddate=enddate, interpolate=interpolate,
             outputfile=outputfile, fill_value=fill_value, ftimestep=ftimestep,
             verbose=False)

# The class can be used from within Python:
from ascii2musica import ascii2Musica
musica = ascii2Musica(namelist=namelist,
                      site_latitude=site_latitude,
                      site_longitude=site_longitude,
                      site_altitude=site_altitude,
                      forcing_height=forcing_height,
                      time2gmt=time2gmt,
                      verbose=False)
musica.read_data(csvfile,
                 startdate=startdate, enddate=enddate,
                 interpolate=interpolate)
musica.convert_units_musica()
musica.write_netcdf(outputfile, fill_value=fill_value,
                    ftimestep=ftimestep)


History
-------
   * Written Oct 2017 by Matthias Cuntz @ INRA
   * rm case sensitivity of header line, Jan 2018, Matthias Cuntz
   * added --sigvars, , Apr 2018, Matthias Cuntz
   * added -f, --fill_value, Apr 2018, Matthias Cuntz
   * added -t, --ftimestep, Jun 2018, Matthias Cuntz
   * changed calc time step cummulating hours in years,
     Jun 2018, Matthias Cuntz
   * make script as well as python routine, Nov 2018, Matthias Cuntz
   * allow water vapour mixing ratio (kg/kg) for humidity,
     Nov 2018, Matthias Cuntz
   * year,month,day,hour,minute,second dates did not work,
     Nov 2018, Matthias Cuntz
   * replace exec by use of dictionaries, Nov 2018, Matthias Cuntz
   * allow older netcdf4 implementation, , Nov 2018, Matthias Cuntz
   * bug in sigvars: search lowercase name, Nov 2018, Matthias Cuntz
   * bug enddate: used isyear instead of ieyear, Feb 2019, Julien Sainte-Marie
   * daily in/out possible, Nov 2019, Matthias Cuntz
   * special input from Dietrich et al. (Annals of Forest Science, 2019),
     Nov 2019, Matthias Cuntz
   * bug in year0, search for year column instead of assuming column 0,
     Jan 2020, Matthias Cuntz & James Ryder
   * Remove calc of Qair in case of daily input data: done in MuCSIA,
     Jan 2020, Matthias Cuntz & James Ryder
   * remove np.int, etc., Apr 2022, Matthias Cuntz
   * allow year,month,day,time, Aug 2022, Matthias Cuntz
   * interpolate, Aug 2022, Matthias Cuntz
   * assert type(dtime) == float if ftimestep != 0.5, Sep 2022, Matthias Cuntz
   * bug in creation of netcdf file: did not write z-dimension,
     Feb 2023, Matthias Cuntz
   * recalculate year0 if ftimestep != 0.5 to avoid times < 0,
     Apr 2023, Matthias Cuntz
   * bug in creation of netcdf variables Wind_N and Wind_E, were only 3D,
     Apr 2023, Matthias Cuntz
   * bug: use key 'tair' instead of 'wind' if wind_dir not given,
     May 2024, Matthias Cuntz
   * Use class_ascii2netcdf, Dec 2025, Matthias Cuntz
   * Removed --sigvars, Dec 2025, Matthias Cuntz
   * Added __docstring__, Apr 2026, Matthias Cuntz
   * Added --verbose, Apr 2026, Matthias Cuntz
   * Added ascii2musica function, Apr 2026, Matthias Cuntz

"""
import time as ptime
import cftime as cf
import netCDF4 as nc
import numpy as np
from class_ascii2netcdf import ascii2Netcdf


__all__ = ['ascii2Musica', 'ascii2musica', '__docstring__']


"""
Development

py ascii2musica.py -s 2005-01-01T00:30 -e 2006-01-01T00:30 -o MuSICA_site_in_FR-Hes_2005-01.nc -t 1.0 --site_latitude=48.6742166667 --site_longitude=7.06461666667 --time2gmt=0 --forcing_height=22.5 -f=-9999 MuSICA_site_in_FR-Hes_2005.csv

cdo --sortname copy MuSICA_site_in_FR-Hes_2005.nc test0.nc
cdo --sortname copy MuSICA_site_in_FR-Hes_2005-01.nc test1.nc
cdo diffn test0.nc test1.nc
cdo diffn,rellim=1e-09,abslim=1e-5 test0.nc test1.nc 

# mesocosm
meso='/Users/cuntz/prog/forge.inrae/musica/tests/mesocosm_Lolium_Regina_Juan'
py ascii2musica.py -s 2017-06-15T10:00:00 -n ${meso}/musica.nml -o ${meso}/inout_files/MuSICA_in_mesocosm-02.nc -f=-9999 ${meso}/inout_files/MuSICA_in_run2_ch2.csv
cdo --sortname copy ${meso}/inout_files/MuSICA_in_mesocosm.nc test0.nc
cdo --sortname copy ${meso}/inout_files/MuSICA_in_mesocosm.nc test2.nc
cdo diffn test0.nc test2.nc 

"""


__docstring__ = """
Write MuSICA netcdf forcing file from csv input file.

The script predicts the csv file structure from the first two lines, i.e.
variable names and units, field delimiter and date/time variables.

Apart from date/time information
(year,jday,time or year,month,day,hour,minute),
the script makes one netcdf variable for each column found in the csv file.

Additional information such as forcing height, site latitude, etc. can be given
as optional inputs or they can be read from the musica namelist (case
in-sensititve). If both are given, the information from the optional input
parameters will be taken.

Header
- - - 
Possible header structures have one or two lines:
  - 1 header line with only variables names that will be taken as netcdf
      variable names, e.g.
      year,jday,time,windspeed, ...
  - 1 header line with variables names followed by units with our without
      parentheses () or brackets [], e.g.
      year (years),jday (days),time (hours),windspeed (m/s), ...
  - 2 header lines, the first with the variables names and the second with
      the units with our without parentheses () or brackets [], e.g.
      year,jday,time,windspeed, ...
      years,days,hours,m/s, ...

Field Delimiter
- - - - - - - -
Commas, semicolons and whitespace are allowed as delimiters between columns.
If whitespace is used than only the first and last header structures are
allowed, no 'varname (unit)'.

Date/Time
- - - - -
Date/time information can be given in two formats:
  - year,jday,time.
    jday is the doy-of-year starting at 1 and time is the fractional time
    starting at 0.
    doy can be used in the header line instead of jday.
    For example for date/time 18 August 2009 at 13:15 h
    2009,230,13.25
  - year,month,day,hour,minute,second. Seconds can be omitted;
    assumed 0 in this case. That means for the 18 August 2009 at 13:15 h
    2009,8,18,13,15,00
    or
    2009,8,18,13,15
The script looks for month, and jday or doy (case in-sensitive), respectively.
If both are given, month takes precedence if the script finds also day, hour
and minute.
If year,month,day is given but no hour,minute,second, or if year,jday but not
time is given, then the script assumes daily input.

Mandatory Variables
- - - - - - - - - -
Apart from date/time information
(year,jday,time  or  year,month,day,hour,minute),
the script makes one netcdf variable for each column found in the csv file.
However, MuSICA needs the following forcing data (all at reference level) and
so these variables have to be present if sub-diurnal data:
       swdown | Incoming solar radiation   | W/m2
       lwdown | Incoming thermal radiation | W/m2
       precip | Precipitation              | kg/m2/s
         tair | Air temperature            | K
    windspeed | Wind speed                 | m/s
        psurf | Atmospheric pressure       | Pa
       co2air | CO2 mixing ratio           | ppmv
        rhair | Air relative humidity      | %
              or
         qair | water vapour mixing ratio  | kg/kg
If other units are given, the script tries to convert them.
If no units are given, the scripts assumes that the variables have already
the right units.
If precip is given, rain and snow will be split at 1 degree Celsius
air temperature.

Optional Variables
------------------
      winddir | Wind direction             | degree
can be used for northerly and easterly winds (together with windspeed).

Optional forcing data is needed by MuSICA dependent on model flags:
if neumann_water_bot_flag = .false.
      wtdepth | Water table depth          | m

if ozone_flag = .true.
   o3_air_ref | Ozone mixing ratio         | ppt

if ocs_flag = .true. (not used yet)
  ocs_air_ref | OCS mixing ratio           | ppt

if isotope_flag = .true.
  d_co2_air_c13_ref | 13C isotopic composition of atmospheric CO2                | permilVPDB
  d_co2_air_o18_ref | 18O isotopic composition of atmospheric CO2                | permilVPDB-CO2
  d_vapour_ox18_ref | 18O isotopic composition of atmospheric water vapour       | permilVSMOW
  d_vapour_deut_ref | Deuterium isotopic composition of atmospheric water vapour | permilVSMOW
  d_rain_ox18_ref   | 18O isotopic composition of precipitation                  | permilVSMOW
  d_rain_deut_ref   | Deuterium isotopic composition of precipitation            | permilVSMOW

Additional Information
----------------------
Additional information about the station must be known. These are:
  - forcing_height: the reference height where meteorological variables and
    fluxes were measured (in meters)
  - site_latitude: the site latitude (in degrees, -90 to +90)
  - site_longitude: the site longitude (in degrees, -180 to +180)
  - site_altitude: the site altitude (in m)
  - time2gmt: the time shift (in hours) between date/time in the csv file and
    GMT=UTC (from -12 to +12)
These can be given on the command line by their name preceded by two dashes,
e.g. --site_latitude 17.52 or --site_latitude=17.52 or they can be read from
the musica namelist (case in-sensititve) with the -n switch, i.e. -n
musica.nml. If both are given, the information from the command line will be
taken.

"""


class ascii2Musica(ascii2Netcdf):
    """
    Class for making MuSICA netCDF forcing file from ASCII data

    Parameters
    ----------
    namelist : str, optional
        musica.nml namelist file.
    site_latitude : float
        The site's latitude from -90 to 90.
    site_longitude : float
        The site's longitude from -180 to 180.
    site_altitude : float, optional
        The site's altitude (default: 0.).
    forcing_height : float, optional
        The height were temperature and humidity were measured in m
        (default: 2.).
    time2gmt : float, optional
        The time shift (in hours) between date/time in the csv file and
        GMT=UTC (from -12 to +12, default: 0).

    """

    def __init__(self, *args,
                 namelist='',
                 site_latitude=None, site_longitude=None, site_altitude=None,
                 forcing_height=None, time2gmt=None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        
        self.site_latitude = None
        self.site_longitude = None
        self.site_altitude = 0.
        self.forcing_height = 2.
        self.time2gmt = 0.

        if namelist:
            nlist = ['site_latitude', 'site_longitude', 'site_altitude',
                     'latitude', 'longitude', 'altitude',
                     'forcing_height', 'time2gmt']
            with open(namelist, 'r') as fn:
                for line in fn:
                    if '=' in line:
                        ll = line.strip().split('=')
                        ll0 = ll[0].strip().lower()
                        if ll0 in nlist:
                            if eval(ll0) is None:  # None is default
                                if ll0.endswith('latitude'):
                                    self.site_latitude = float(ll[1])
                                elif ll0.endswith('longitude'):
                                    self.site_longitude = float(ll[1])
                                elif ll0.endswith('altitude'):
                                    self.site_altitude = float(ll[1])
                                elif ll0 == 'forcing_height':
                                    self.forcing_height = float(ll[1])
                                elif ll0 == 'time2gmt':
                                    self.time2gmt = float(ll[1])

        if site_latitude is None:
            if self.site_latitude is None:
                raise ValueError('site_latitude must be given.')
        else:
            self.site_latitude = site_latitude

        if site_longitude is None:
            if self.site_longitude is None:
                raise ValueError('site_longitude must be given.')
        else:
            self.site_longitude = site_longitude

        if site_altitude is not None:
            self.site_altitude = site_altitude

        if forcing_height is not None:
            self.forcing_height = forcing_height

        if time2gmt is not None:
            self.time2gmt = time2gmt

        return

    #
    # Output variables names and units (replace standard)
    #

    def output_variables(self):
        """
        Output variable names, longnames, and units

        """
        # name of downwelling shortwave radiation 
        self.swdown = 'SWdown'
        self.lswdown = self.swdown.lower()
        # name of windspeed
        self.windspeed = 'Wind_N'
        self.lwindspeed = self.windspeed.lower()
        # output variables and units
        self.varnames = ['CO2air', self.swdown, 'LWdown', 'PSurf',
                         'Qair', 'Rainf', 'Snowf', 'Tair',
                         self.windspeed, 'Wind_E']  # , 'h_sbl']
        self.varlownames = [ vv.lower() for vv in self.varnames ]
        self.varunits = ['ppmv', 'W/m2', 'W/m2', 'Pa',
                         'kg/kg', 'kg/m2/s', 'kg/m2/s', 'K',
                         'm/s', 'm/s']  # ,'m']
        self.varlongnames = ['Near surface CO2 concentration',
                             'Surface incident shortwave radiation',
                             'Surface incident longwave radiation',
                             'Surface pressure',
                             'Near surface specific humidity',
                             'Rainfall rate',
                             'Snowfall Rate',
                             'Near surface air temperature',
                             'Near surface northerly wind speed',
                             'Near surface easterly wind speed']  # ,
                             # 'Surface boundary layer height']

        # dictionary with lowercase variables as keys
        self.dvarunits = dict(zip(self.varlownames, self.varunits))
        self.dvarnames = dict(zip(self.varlownames, self.varnames))
        self.dvarlongnames = dict(zip(self.varlownames, self.varlongnames))

    #
    # Convert variable units
    #

    def convert_units_musica(self):
        """
        Convert variables to output units

        """
        # standard conversions
        self.convert_units()

        # rename SWdown variables
        if 'swdown' in self.livars:
            self.livars[self.livars.index('swdown')] = self.lswdown

        # wind_n / wind_e
        if 'wind_n' in self.livars:
            if 'wind_e' not in self.livars:
                self.df['wind_e'] = np.zeros_like(self.df['wind_n'])  # easterly
        else:
            for vv in ['wind', 'windspeed', 'wind_speed']:
                if vv in self.livars:
                    vwind = vv
                    break
            vdir = ''
            for vv in ['wind_dir', 'winddir']:
                if vv in self.livars:
                    vdir = vv
                    break
            if vdir:
                # northerly / easterly
                self.df['wind_n'] = (self.df[vwind] *
                                     np.sin(np.deg2rad(self.df[vdir])))
                self.df['wind_e'] = (self.df[vwind] *
                                     np.cos(np.deg2rad(self.df[vdir])))
                #MC comment keep Wind_dir at FR-Hes
                self.df.drop(axis=0, columns=[vdir], inplace=True)
                self.livars.remove(vdir)
                #MC uncomment keep Wind_dir at FR-Hes
                # if 'wind_dir' in self.livars:
                #     self.ivars['wind_dir'] = 'Wind_dir'
                #MC keep Wind_dir at FR-Hes
            else:
                self.df['wind_n'] = self.df[vwind]                    # northerly
                self.df['wind_e'] = np.zeros_like(self.df['wind_n'])  # easterly
            if vwind not in ['wind_n', 'wind_e']:
                self.df.drop(axis=0, columns=[vwind], inplace=True)
                self.livars.remove(vwind)

        # #MC uncomment upcase wtdepth for mesocosm
        # if 'wtdepth' in self.livars:
        #     self.ivars['wtdepth'] = 'WTDEPTH'
        # #MC uncomment upcase wtdepth for mesocosm

        return

    #
    # Write output file
    #

    def write_netcdf(self, outputfile='', fill_value=None, ftimestep=0.5):
        """
        Write output netcdf file

        Parameters
        ----------
        ouputfile : string, optional
            Name of netcdf output file
            (default: suffix of infile replaced by .nc).
        fill_value : float, optional
            _FillValue for all variables in netcdf file
            (default: None, i.e. netCDF defaults).
        ftimestep : float, optional
            Fraction of time step at which time is given:
              0: beginning of time step

              0.5: middle of time step

              1: end of time step

            Note that MuSICA is using the middle of the time step, so
            output times will be transformed to the middle of the
            timestep. (Default: 0.5)

        """        
        if outputfile == '':
           outputfile = self.ifile[0:self.ifile.rfind('.')] + '.nc'
        self.outputfile = outputfile
        self.ftimestep = ftimestep

        # shift time to middle of time step
        if self.ftimestep != 0.5:
            dt = np.timedelta64(self.dt, 's')
            indx = self.df.index + (0.5 - ftimestep) * dt
            self.df.set_index(indx, inplace=True)

        # NETCDF
        if self.verbose:
            print('Create netcdf file ', outputfile)
        with nc.Dataset(outputfile, 'w', format='NETCDF3_64BIT') as fo:

            # Structure
            x    = fo.createDimension('x', 1)
            y    = fo.createDimension('y', 1)
            z    = fo.createDimension('z', 1)
            time = fo.createDimension('time', None)

            nav_lon            = fo.createVariable('nav_lon', 'f4', ('y', 'x',))
            nav_lon.long_name  = 'Longitude'
            nav_lon.units      = 'degrees_east'
            nav_lon.valid_min  = -180.
            nav_lon.valid_max  = 180.

            nav_lat            = fo.createVariable('nav_lat', 'f4', ('y', 'x',))
            nav_lat.long_name  = 'Latitude'
            nav_lat.units      = 'degrees_north'
            nav_lat.valid_min  = -90.
            nav_lat.valid_max  = 90.

            x                  = fo.createVariable('x', 'i4', ('x',))
            y                  = fo.createVariable('y', 'i4', ('y',))
            level              = fo.createVariable('level', 'f4', ('z',))

            vtime              = fo.createVariable('time', 'f8', ('time',))
            vtime.long_name    = 'Time'
            vtime.units        = f'hours since {self.df.index[0].year}-01-01 00:00:00'
            vtime.calendar     = 'gregorian'
            vtime.time_origin  = (f'{self.df.index[0].year}-JAN-01 00:00:00 (GMT+'
                                  f'{self.time2gmt})')

            vloc = dict()
            # Write all columns not only the MuSICA forcing columns
            # for iv in self.dvarnames.keys():
            for iv in self.df.columns:
                if iv in ['tair', 'qair', 'co2air', 'wind_n', 'wind_e']:
                    vloc[iv] = fo.createVariable(self.dvarnames[iv], 'f4',
                                                 ('time', 'z', 'y', 'x',),
                                                 fill_value=fill_value)
                else:
                    if iv in self.dvarnames:
                        vloc[iv] = fo.createVariable(self.dvarnames[iv], 'f4',
                                                     ('time', 'y', 'x',),
                                                     fill_value=fill_value)
                    elif iv in self.livars:
                        vloc[iv] = fo.createVariable(self.ivars[iv], 'f4',
                                                     ('time', 'y', 'x',),
                                                     fill_value=fill_value)
                    else:
                        vloc[iv] = fo.createVariable(iv, 'f4',
                                                     ('time', 'y', 'x',),
                                                     fill_value=fill_value)
                if iv in self.dvarnames:
                    if self.dvarunits[iv]:
                        vloc[iv].units = self.dvarunits[iv]
                    if self.dvarlongnames[iv]:
                        vloc[iv].long_name = self.dvarlongnames[iv]
                elif iv in self.livars:
                    if self.iunits[iv]:
                        vloc[iv].units = self.iunits[iv]

            fo.conventions = 'GDT 1.2'
            fo.file_name   = outputfile
            fo.production  = ('NetCDF file generated from ASCII file ' +
                              self.ifile + ' on ' +
                              ptime.strftime('%a %d %B %Y %H:%M:%S'))

            # calc netCDF time
            cdates = [ cf.datetime(dd.year, dd.month, dd.day, dd.hour,
                                   dd.minute, dd.second)
                       for dd in self.df.index ]
            dtime = cf.date2num(cdates,
                                f'hours since {self.df.index[0].year}-01-01 00:00:00',
                                calendar='gregorian')

            # Write variables
            nav_lon[:] = self.site_longitude
            nav_lat[:] = self.site_latitude
            x[:]       = 1
            y[:]       = 1
            level[:]   = self.forcing_height
            vtime[:]   = dtime

            # for iv in self.dvarnames.keys():
            for iv in self.df.columns:
                vloc[iv][:] = self.df[iv]

        return


# -------------------------------------------------------------------------
# Function using class
#

def ascii2musica(csvfile,
                 namelist='', site_latitude=None, site_longitude=None,
                 site_altitude=None, forcing_height=None, time2gmt=None,
                 startdate='', enddate='', interpolate=False,
                 outputfile='', fill_value=None, ftimestep=0.5,
                 verbose=False, **kwargs):
    """
    Write MuSICA netcdf forcing file from csv input file.

    See variable `__docstring__` or use --help on command line.

    Parameters
    ----------
    csvfile : string
        csv file with MuSICA's forcing variables
    namelist : str, optional
        musica.nml namelist file.
    site_latitude : float
        The site's latitude from -90 to 90.
    site_longitude : float
        The site's longitude from -180 to 180.
    site_altitude : float, optional
        The site's altitude (default: 0.).
    forcing_height : float, optional
        The height were temperature and humidity were measured in m
        (default: 2.).
    time2gmt : float, optional
        The time shift (in hours) between date/time in the csv file and
        GMT=UTC (from -12 to +12, default: 0).
    startdate : string, optional
        First possible date in netcdf output file in ISO8601 format.
        (Default: first date in input file)
    enddate : string, optional
        Last possible date in output netcdf file in ISO8601 format.
        (Default: last date in input file)
    interpolate : bool, optional
        Linearly interpolate missing data in input file given as NaN
        or empty cells.
    ouputfile : string, optional
        Name of netcdf output file
        (default: suffix of infile replaced by .nc).
    fill_value : float, optional
        _FillValue for all variables in netcdf file
        (default: None, i.e. netCDF defaults).
    ftimestep : float, optional
        Fraction of time step at which time is given:
          0: beginning of time step

          0.5: middle of time step

          1: end of time step

        Note that MuSICA is using the middle of the time step, so
        output times will be transformed to the middle of the
        timestep. (Default: 0.5)
    verbose : bool, optional
        Report progress if True.
    **kwargs : dict, optional
        All other keyword arguments will be ignored

    """
    musica = ascii2Musica(namelist=namelist,
                          site_latitude=site_latitude,
                          site_longitude=site_longitude,
                          site_altitude=site_altitude,
                          forcing_height=forcing_height,
                          time2gmt=time2gmt,
                          verbose=verbose)
    musica.read_data(csvfile,
                     startdate=startdate, enddate=enddate,
                     interpolate=interpolate)
    musica.convert_units_musica()
    musica.write_netcdf(outputfile, fill_value=fill_value,
                        ftimestep=ftimestep)

    return


# -------------------------------------------------------------------------
# Main
#

if __name__ == '__main__':

    import argparse

    startdate = ''
    enddate = ''
    outputfile = ''
    namelist = ''
    ftimestep = 0.5
    site_latitude = None
    site_longitude = None
    site_altitude = None
    time2gmt = None
    forcing_height = None
    fill_value = None
    interpolate = False
    verbose = False

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__docstring__)
    hstr = ('First possible date in netcdf output file (default: first date'
            ' in input file). Different formats are possible; if seconds or'
            ' time are omitted, 0 will be taken: DD.MM.YYYYTHH:MM:SS,'
            ' YYYY-MM-DDTHH:MM:SS, DD/MM/YYYYTHH:MM:SS.')
    parser.add_argument('-s', '--startdate', action='store',
                        default=startdate, dest='startdate',
                        metavar='start_date', help=hstr)
    hstr = ('Last possible date in output netcdf file (default: last date in'
            ' input file). Different formats are possible; if seconds or time'
            ' are omitted, 23 hours and 59 minutes and/or seconds will be'
            ' taken, respectively: DD.MM.YYYYTHH:MM:SS, YYYY-MM-DDTHH:MM:SS,'
            ' DD/MM/YYYYTHH:MM:SS.')
    parser.add_argument('-e', '--enddate', action='store', default=enddate,
                        dest='enddate', metavar='end_date', help=hstr)
    hstr = ('Name of netcdf output file (default: suffix of infile replaced'
            ' by .nc).')
    parser.add_argument('-o', '--ouputfile', action='store',
                        default=outputfile, dest='outputfile', metavar='outputfile',
                        help=hstr)
    parser.add_argument('-n', '--namelist', action='store',
                        default=namelist, dest='namelist', metavar='namelist',
                        help='MuSICA namelist file.')
    hstr = ('Fraction of time step at which time is given (default: 0.5):'
            ' 0=beginning of time step, 0.5=middle of time step, 1=end of'
            ' time step. Note that MuSICA is using the middle of the time'
            ' step so output times will be transformed to the middle of the'
            ' timestep.')
    parser.add_argument('-t', '--ftimestep', action='store', type=float,
                        default=ftimestep, dest='ftimestep',
                        metavar='fractional time step', help=hstr)
    parser.add_argument('--site_latitude', action='store', type=float,
                        default=site_latitude, dest='site_latitude',
                        metavar='latitude',
                        help="The site's latitude from -90 to 90.")
    parser.add_argument('--site_longitude', action='store', type=float,
                        default=site_longitude, dest='site_longitude',
                        metavar='longitude',
                        help="The site's longitude from -180 to 180.")
    parser.add_argument('--site_altitude', action='store', type=float,
                        default=site_altitude, dest='site_altitude',
                        metavar='altitude',
                        help="The site's altitude from -180 to 180.")
    hstr = ('Time difference between local time and GMT=UTC from -12 to 12.'
            ' (Default: 0.)')
    parser.add_argument('--time2gmt', action='store', type=float,
                        default=time2gmt, dest='time2gmt', metavar='time2gmt',
                        help=hstr)
    hstr = ('The height were the meteorological variables were measured in m.'
            ' (Default: 2.)')
    parser.add_argument('--forcing_height', action='store', type=float,
                        default=forcing_height, dest='forcing_height',
                        metavar='forcing_height', help=hstr)
    hstr = ('_FillValue for all variables in netcdf file. Negative fill'
            ' values have to be quoted. (default: None).')
    parser.add_argument('-f', '--fill_value', action='store', type=float,
                        default=fill_value, dest='fill_value',
                        metavar='fill_value', help=hstr)
    hstr = ('Linearly interpolate missing data in input file given as NaN'
            ' or empty cells.')
    parser.add_argument('--interpolate', action='store_true',
                        default=interpolate, dest='interpolate', help=hstr)
    parser.add_argument('-v', '--verbose', action='store_true',
                        default=verbose, dest='verbose',
                        help="Output progress on screen.")
    parser.add_argument('infile', nargs='?', default=None,
                        metavar='csv_forcing_file',
                        help="csv file with MuSICA's forcing variables.")

    args = parser.parse_args()

    assert args.infile is not None, (
        'csv file with MuSICA forcing variables must be given.')

    t1 = ptime.time()

    musica = ascii2Musica(namelist=args.namelist,
                          site_latitude=args.site_latitude,
                          site_longitude=args.site_longitude,
                          site_altitude=args.site_altitude,
                          forcing_height=args.forcing_height,
                          time2gmt=args.time2gmt,
                          verbose=True)
    musica.read_data(args.infile,
                     startdate=args.startdate, enddate=args.enddate,
                     interpolate=args.interpolate)
    musica.convert_units_musica()
    musica.write_netcdf(args.outputfile, fill_value=args.fill_value,
                        ftimestep=args.ftimestep)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.)
             if (t2 - t1) > 60. else '[s]: {:d}'.format(int(t2 - t1)))
    print('Time elapsed: ', strin)
