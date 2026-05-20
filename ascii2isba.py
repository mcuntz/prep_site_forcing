#!/usr/bin/env python
'''
Write ISBA netcdf forcing file from csv input file.

# Script can be called on the command line:
python ascii2isba.py --latitude=48.6742166667 --longitude=7.06461666667 --altitude=300.0 --outputfile=ISBA_in_FR-Hes_1997-era5-1.nc --ftimestep=1.0 --reference_height=22.5 --reference_height_wind=22.5 --slope=0.0 --aspect=0.0 --fill_value=-9999999.0 -v ISBA_in_FR-Hes_1997-era5.csv

# The wrapper function can be called from within Python:
from ascii2isba import ascii2isba
ascii2isba(csvfile,
           latitude=latitude, longitude=longitude, altitude=altitude,
           reference_height=reference_height,
           reference_height_wind=reference_height_wind,
           slope=slope, aspect=aspect,
           startdate=startdate, enddate=enddate,
           outputfile=outputfile, fill_value=fill_value,
           ftimestep=ftimestep, verbose=False)

# The class can be used from within Python:
from ascii2isba import ascii2Isba
isba = ascii2Isba(latitude=latitude, longitude=longitude,
                  altitude=altitude,
                  reference_height=reference_height,
                  reference_height_wind=reference_height_wind,
                  slope=slope, aspect=aspect,
                  verbose=False)
isba.read_data(csvfile,
               startdate=startdate, enddate=enddate)
isba.convert_units_isba()
isba.write_netcdf(outfile, fill_value=fill_value,
                  ftimestep=ftimestep)


History
-------
   * Written Jan 2023 by Matthias Cuntz, UMR Silva, INRAE
   * Removed cfunits package, Nov 2025, Matthias Cuntz
   * Rewritten as class, Nov 2025, Matthias Cuntz
   * Use pandas for reading ascii data, Dec 2025, Matthias Cuntz
   * Use class_ascii2netcdf, Dec 2025, Matthias Cuntz

'''
import time as ptime
import cftime as cf
import netCDF4 as nc
import numpy as np
from class_ascii2netcdf import ascii2Netcdf


__all__ = ['ascii2Isba', 'ascii2isba', '__docstring__']


"""
Development

py ascii2isba.py --latitude=48.6742166667 --longitude=7.06461666667 --altitude=300.0 --outputfile=ISBA_in_FR-Hes_1997-era5-1.nc --ftimestep=1.0 --reference_height=22.5 --reference_height_wind=22.5 --slope=0.0 --aspect=0.0 --fill_value=-9999999.0 -v ISBA_in_FR-Hes_1997-era5.csv

cdo diffn ISBA_in_FR-Hes_1997-era5.nc ISBA_in_FR-Hes_1997-era5-1.nc
cdo diffn,rellim=1e-09,abslim=1e-5 ISBA_in_FR-Hes_1997-era5.nc ISBA_in_FR-Hes_1997-era5-1.nc
"""


__docstring__ = """
Write ISBA netcdf forcing file from csv input file.

The script predicts the csv file structure from the first two lines, i.e.
variable names and units, the field delimiter and the date/time variables.

Apart from the date/time information
(year,jday,time or year,month,day,time or year,month,day,hour,minute),
the script makes one netcdf variable for each column known from its name.

Site latitude and longitude have to be given as optional inputs.

Header
- - -
Possible header structures have one or two lines:
  - 1 header line with only variables names that will be taken as netcdf
      variable names, e.g.
      year,jday,time,wind, ...
  - 1 header line with variables names followed by units with our without
      parentheses () or brackets [], e.g.
      year (years),jday (days),time (hours),wind (m/s), ...
  - 2 header lines, the first with the variables names and the second with
      the units with our without parentheses () or brackets [], e.g.
      year,jday,time,wind, ...
      years,days,hours,m/s, ...

Field Delimiter
- - - - - - - -
Commas, semicolons and whitespace are allowed as delimiters between
columns.
If whitespace is used then only the first and last header structures are
allowed, i.e. no 'varname (unit)' allowed in this case.

Date/Time
- - - - -
Date/time information can be given in three formats:
  - year,jday,time.
    jday is the doy-of-year starting at 1 and time is the fractional time
    starting at 0. doy can be used in the header line instead of jday.
    For example for date/time 18 August 2009 at 13:15 h
    2009,230,13.25
  - year,month,day,time.
    time is the fractional time starting at 0.
    For example for date/time 18 August 2009 at 13:15 h
    2009,8,18,13.25
  - year,month,day,hour,minute,second. Seconds can be omitted;
    assumed 0 in this case. That means for the 18 August 2009 at 13:15 h
    2009,8,18,13,15,0
    or
    2009,8,18,13,15
The script looks for month, and jday or doy (case in-sensitive),
respectively.
If both are given, month takes precedence if the script finds also day,
hour and minute, or if it finds time.

Mandatory Variables
- - - - - - - - - -
Apart from the date/time information
(year,jday,time  or  year,month,day,time  or  year,month,day,hour,minute),
the script makes one netcdf variable for each column known from its name.
ISBA needs the following forcing data:
       co2air | CO2 concentration                | kg/m3
       swdown | Incoming solar radiation         | W/m2
       lwdown | Incoming thermal radiation       | W/m2
        psurf | Atmospheric pressure             | hPa
         qair | Specific humidity                | kg/kg
              or
        rhair | Air relative humidity            | %
       precip | Precipitation                    | kg/m2/s
              or
       rainf  | Liquid precipitation             | kg/m2/s
              and
       snowf  | Solid precipitation              | kg/m2/s
   sca_swdown | Incoming diffuse solar radiation | W/m2
         tair | Air temperature                  | K
   wind_speed | Wind speed                       | m/s
and optionally
   sca_swdown | Incoming diffuse solar radiation | W/m2
     wind_dir | Wind direction                   | degree
If other units are given, the script tries to convert them.
If no units are given, the scripts assumes that the variables have already
the right units.
If precip is given, rain and snow will be split at 1 degree Celsius
air temperature.

"""


# -------------------------------------------------------------------------
# Class
#

class ascii2Isba(ascii2Netcdf):
    """
    Class for making ISBA netCDF forcing file from ASCII data

    Parameters
    ----------
    latitude : float
        The site's latitude from -90 to 90.
    longitude : float
        The site's longitude from -180 to 180.
    altitude : float, optional
        The site's elevation above sea level in m (default: 0.)
    reference_height : float
        The height were temperature and humidity were measured in m
        (default: 2.).
    reference_height_wind : float, optional
        The height were wind was measured in m (default: ``reference_height``).
    slope : float, optional
        Slope of terrain in degrees from horizontal (default: 0.).
    aspect : float, optional
        Aspect of the slope in degrees from north (default: 0.).

    """

    def __init__(self, *args,
                 latitude=None, longitude=None, altitude=0.,
                 reference_height=2., reference_height_wind=None,
                 slope=0., aspect=0., **kwargs):
        super().__init__(*args, **kwargs)

        if latitude is None:
            raise ValueError("Latitude must be given.")
        else:
            self.latitude = latitude

        if longitude is None:
            raise ValueError("Longitude must be given.")
        else:
            self.longitude = longitude

        self.altitude = altitude

        self.reference_height = reference_height

        if reference_height_wind is None:
            self.reference_height_wind = self.reference_height
        else:
            self.reference_height_wind = reference_height_wind

        self.slope = slope
        self.aspect = aspect

        return

    #
    # Output variables names and units (replace standard)
    #

    def output_variables(self):
        """
        Output variable names, longnames, and units

        """
        # name of downwelling shortwave radiation
        self.swdown = 'DIR_SWdown'
        self.lswdown = self.swdown.lower()
        # name of windspeed
        self.windspeed = 'Wind'
        self.lwindspeed = self.windspeed.lower()
        # output variables and units
        self.varnames = ['CO2air', 'DIR_SWdown', 'LWdown', 'PSurf',
                         'Qair', 'Rainf', 'SCA_SWdown', 'Snowf',
                         'Tair', 'Wind', 'Wind_DIR']
        self.varlownames = [ vv.lower() for vv in self.varnames ]
        self.varunits = ['kg/m3', 'W/m2', 'W/m2', 'Pa', 'kg/kg', 'kg/m2/s',
                         'W/m2', 'kg/m2/s', 'K', 'm/s', 'degree']
        self.varlongnames = ['Near Surface CO2 Concentration',
                             'Surface Incident Direct Shortwave Radiation',
                             'Surface Incident Longwave Radiation',
                             'Surface Pressure',
                             'Near Surface Specific Humidity',
                             'Rainfall Rate',
                             'Surface Incident Diffuse Shortwave Radiation',
                             'Snowfall Rate',
                             'Near Surface Air Temperature',
                             'Wind Speed',
                             'Wind Direction']

        # dictionary with lowercase variables as keys
        self.dvarunits = dict(zip(self.varlownames, self.varunits))
        self.dvarnames = dict(zip(self.varlownames, self.varnames))
        self.dvarlongnames = dict(zip(self.varlownames, self.varlongnames))

    #
    # Convert variable units
    #

    def convert_units_isba(self):
        """
        Convert variables to output units

        """
        # standard conversions
        self.convert_units()

        # SWdown
        if 'swdown' in self.livars:
            self.ivars['swdown'] = self.swdown
            self.df.rename(columns={'swdown': self.lswdown}, inplace=True)

        for vv in ['wind', 'windspeed', 'wind_speed']:
            if vv in self.livars:
                self.ivars[vv] = self.lwindspeed
                self.df.rename(columns={vv: self.lwindspeed}, inplace=True)

        # Diffuse radiation
        if 'sca_swdown' not in self.livars:
            self.df['sca_swdown'] = self.df['tair'].copy()
            self.df.loc[:, 'sca_swdown'] = 0.

        # Wind direction
        if 'wind_dir' not in self.livars:
            self.df['wind_dir'] = self.df['tair'].copy()
            self.df.loc[:, 'wind_dir'] = 0.

        return

    #
    # Write output file
    #

    def write_netcdf(self, outputfile='', fill_value=-9999999., ftimestep=1.0):
        """
        Write output netcdf file

        Parameters
        ----------
        ouputfile : string, optional
            Name of netcdf output file
            (default: suffix of infile replaced by .nc).
        fill_value : float, optional
            _FillValue for all variables in netcdf file (default: -9999999.).
        ftimestep : float, optional
            Fraction of time step at which time is given:
              0: beginning of time step

              0.5: middle of time step

              1: end of time step

            Note that ISBA is using the end of the time step, so
            output times will be transformed to the end of the
            timestep. (Default: 1.0)

        """
        if outputfile == '':
           outputfile = self.ifile[0:self.ifile.rfind(".")] + '.nc'
        self.outputfile = outputfile
        self.ftimestep = ftimestep

        # shift time to end of time step
        if self.ftimestep != 1.0:
            dt = np.timedelta64(self.dt, 's')
            indx = self.df.index + (1.0 - ftimestep) * dt
            self.df.set_index(indx, inplace=True)

        doupper = False
        # NETCDF
        if self.verbose:
            print('Create netcdf file ', outputfile)
        with nc.Dataset(outputfile, 'w', format='NETCDF3_64BIT') as fo:
            # Structure
            nb = fo.createDimension('Number_of_points', 1)
            time = fo.createDimension('time', None)

            lon = fo.createVariable('LON', 'f8', ('Number_of_points',),
                                    fill_value=fill_value)
            lon.long_name = "longitude"
            lon.units = "degrees_east"

            lat = fo.createVariable('LAT', 'f8', ('Number_of_points',),
                                    fill_value=fill_value)
            lat.long_name = "latitude"
            lat.units = "degrees_north"

            zs = fo.createVariable('ZS', 'f8', ('Number_of_points',),
                                   fill_value=fill_value)
            zs.long_name = "altitude"
            zs.units = "m"

            if doupper:
                vtime = fo.createVariable('TIME', 'f8', ('time',),
                                          fill_value=fill_value)
            else:
                vtime = fo.createVariable('time', 'f8', ('time',),
                                          fill_value=fill_value)
            vtime.long_name = "time"
            vtime.units = f"hours since {self.df.index[0].year}-01-01 00:00:00"

            uref = fo.createVariable('UREF', 'f8', ('Number_of_points',),
                                     fill_value=fill_value)
            uref.long_name = "Reference_Height_for_Wind"
            uref.units = "m"

            zref = fo.createVariable('ZREF', 'f8', ('Number_of_points',),
                                     fill_value=fill_value)
            zref.long_name = "Reference_Height"
            zref.units = "m"

            if doupper:
                vaspect = fo.createVariable('ASPECT', 'f8',
                                            ('Number_of_points',),
                                            fill_value=fill_value)
            else:
                vaspect = fo.createVariable('aspect', 'f8',
                                            ('Number_of_points',),
                                            fill_value=fill_value)
            vaspect.long_name = "slope aspect"
            vaspect.units = "degrees from north"

            if doupper:
                vslope = fo.createVariable('SLOPE', 'f8',
                                           ('Number_of_points',),
                                           fill_value=fill_value)
            else:
                vslope = fo.createVariable('slope', 'f8',
                                           ('Number_of_points',),
                                           fill_value=fill_value)
            vslope.long_name = "slope angle"
            vslope.units = "degrees from horizontal"

            frc_time_stp = fo.createVariable('FRC_TIME_STP', 'f8',
                                             fill_value=fill_value)
            frc_time_stp.long_name = "Forcing_Time_Step"
            frc_time_stp.units = "s"

            vloc = dict()
            for iv in self.dvarnames.keys():
                if doupper:
                    vloc[iv] = fo.createVariable(iv.upper(), 'f8',
                                                 ('time', 'Number_of_points',),
                                                 fill_value=fill_value)
                else:
                    vloc[iv] = fo.createVariable(self.dvarnames[iv], 'f8',
                                                 ('time', 'Number_of_points',),
                                                 fill_value=fill_value)
                vloc[iv].units = self.dvarunits[iv]
                vloc[iv].long_name = self.dvarlongnames[iv]

            fo.file_name   = outputfile
            fo.production  = ("NetCDF file generated from ASCII file " +
                              self.ifile + " on " +
                              ptime.strftime("%a %d %B %Y %H:%M:%S"))

            # Write variables

            # check radiation
            if 'sca_swdown' in self.livars:
                iidat = self.df['sca_swdown']
                ii = np.where(iidat == fill_value)[0]
                if ii.size > 0:
                    iidat[ii] = 0.

            # calc netCDF time
            cdates = [ cf.datetime(dd.year, dd.month, dd.day, dd.hour,
                                   dd.minute, dd.second)
                       for dd in self.df.index ]
            dtime = cf.date2num(
                cdates, f'hours since {self.df.index[0].year}-01-01 00:00:00',
                calendar='gregorian')

            lon[:]   = self.longitude
            lat[:]   = self.latitude
            zs[:]    = self.altitude
            vtime[:] = dtime
            uref[:]  = self.reference_height_wind
            zref[:]  = self.reference_height
            vaspect[:] = self.aspect
            vslope[:]  = self.slope
            frc_time_stp[:] = self.dt

            for iv in self.dvarnames.keys():
                vloc[iv][:] = self.df[iv]

        return


# -------------------------------------------------------------------------
# Function using class
#

def ascii2isba(csvfile,
               latitude=None, longitude=None, altitude=0.,
               reference_height=2., reference_height_wind=None,
               slope=0., aspect=0.,
               startdate='', enddate='', interpolate=False,
               outputfile='', fill_value=-9999999., ftimestep=1.0,
               verbose=False, **kwargs):
    """
    Write ISBA netcdf forcing file from csv input file.

    See variable `__docstring__` or use --help on command line.

    Parameters
    ----------
    csvfile : string
        csv file with ISBA's forcing variables
    latitude : float
        The site's latitude from -90 to 90.
    longitude : float
        The site's longitude from -180 to 180.
    altitude : float, optional
        The site's elevation above sea level in m (default: 0.)
    reference_height : float
        The height were temperature and humidity were measured in m
        (default: 2.).
    reference_height_wind : float, optional
        The height were wind was measured in m (default: ``reference_height``).
    slope : float, optional
        Slope of terrain in degrees from horizontal (default: 0.).
    aspect : float, optional
        Aspect of the slope in degrees from north (default: 0.).
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
        _FillValue for all variables in netcdf file (default: -9999999.).
    ftimestep : float, optional
        Fraction of time step at which time is given:
          0: beginning of time step

          0.5: middle of time step

          1: end of time step

        Note that ISBA is using the end of the time step, so
        output times will be transformed to the end of the
        timestep. (Default: 1.0)
    verbose : bool, optional
        Report progress if True.
    **kwargs : dict, optional
        All other keyword arguments will be ignored

    """
    isba = ascii2Isba(latitude=latitude, longitude=longitude,
                      altitude=altitude,
                      reference_height=reference_height,
                      reference_height_wind=reference_height_wind,
                      slope=slope, aspect=aspect,
                      verbose=verbose)
    isba.read_data(csvfile,
                   startdate=startdate, enddate=enddate,
                   interpolate=interpolate)
    isba.convert_units_isba()
    isba.write_netcdf(outputfile, fill_value=fill_value,
                      ftimestep=ftimestep)

    return


# -------------------------------------------------------------------------
# Command line
#

if __name__ == '__main__':

    import argparse

    latitude  = None
    longitude = None
    altitude  = 0.
    startdate = ''
    enddate   = ''
    outputfile = ''
    ftimestep = 1.0
    reference_height = 2.
    reference_height_wind = None
    slope = 0.
    aspect = 0.
    fill_value = -9999999.
    verbose = False

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__docstring__)
    parser.add_argument('--latitude', action='store',
                        default=latitude, dest='latitude',
                        metavar='latitude', type=float,
                        help="The site's latitude from -90 to 90.")
    parser.add_argument('--longitude', action='store',
                        default=longitude, dest='longitude',
                        metavar='longitude', type=float,
                        help="The site's longitude from -180 to 180.")
    parser.add_argument('--altitude', action='store',
                        default=altitude, dest='altitude',
                        metavar='altitude', type=float,
                        help="The site's elevation above sea level in m.")
    hstr = ('First possible date in netcdf output file in ISO8601 format'
            ' (default: first date in input file).')
    parser.add_argument('-s', '--startdate', action='store',
                        default=startdate, dest='startdate',
                        metavar='start_date', help=hstr)
    hstr = ('Last possible date in output netcdf file in ISO8601 format'
            ' (default: last date in input file).')
    parser.add_argument('-e', '--enddate', action='store', default=enddate,
                        dest='enddate', metavar='end_date', help=hstr)
    hstr = ('Name of netcdf output file'
            ' (default: suffix of infile replaced by .nc).')
    parser.add_argument('-o', '--outputfile', action='store',
                        default=outputfile, dest='outputfile',
                        metavar='outputfile', help=hstr)
    hstr = (f'Fraction of time step at which time is given'
            f' (default: {ftimestep}): 0=beginning of time step,'
            f' 0.5=middle of time step, 1=end of time step.'
            f' Note that ISBA is using the end of the time step so output'
            f' times will be transformed to the end of the timestep.')
    parser.add_argument('-t', '--ftimestep', action='store', type=float,
                        default=ftimestep, dest='ftimestep',
                        metavar='fractional time step', help=hstr)
    hstr = (f'The height of the temperature and humidity measurements in m'
            f' (default: {reference_height}).')
    parser.add_argument('--reference_height', action='store',
                        default=reference_height, dest='reference_height',
                        metavar='reference_height', type=float,
                        help=hstr)
    hstr = ('The height of the wind measurements in m'
            ' (default: same as reference_height).')
    parser.add_argument('--reference_height_wind', action='store',
                        default=reference_height_wind, type=float,
                        dest='reference_height_wind',
                        metavar='reference_height_wind',
                        help=hstr)
    hstr = (f'The slope of the terrain in degrees from horizontal'
            f' (default: {slope}).')
    parser.add_argument('--slope', action='store',
                        default=slope, dest='slope',
                        metavar='slope', type=float,
                        help=hstr)
    hstr = (f'The aspect of the slope in degrees from north'
            f' (default: {aspect}).')
    parser.add_argument('--aspect', action='store',
                        default=aspect, dest='aspect',
                        metavar='aspect', type=float,
                        help=hstr)
    hstr = (f'_FillValue for all variables in netcdf file. Negative fill'
            f' values have to be quoted. (default: {fill_value}).')
    parser.add_argument('-f', '--fill_value', action='store', type=float,
                        default=fill_value, dest='fill_value',
                        metavar='fill_value', help=hstr)
    parser.add_argument('-v', '--verbose', action='store_true',
                        default=verbose, dest='verbose',
                        help="Output progress on screen.")
    parser.add_argument('infile', nargs='?', default=None,
                        metavar='csv_forcing_file',
                        help="csv file with ISBA's forcing variables.")

    args = parser.parse_args()

    assert args.infile is not None, (
        'csv file with ISBA forcing variables must be given.')

    t1 = ptime.time()

    isba = ascii2Isba(latitude=args.latitude, longitude=args.longitude,
                      altitude=args.altitude,
                      reference_height=args.reference_height,
                      reference_height_wind=args.reference_height_wind,
                      slope=args.slope, aspect=args.aspect,
                      verbose=args.verbose)
    isba.read_data(args.infile,
                   startdate=args.startdate, enddate=args.enddate)
    isba.convert_units_isba()
    isba.write_netcdf(args.outputfile, fill_value=args.fill_value,
                      ftimestep=args.ftimestep)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.)
             if (t2 - t1) > 60. else '[s]: {:d}'.format(int(t2 - t1)))
    print('Time elapsed: ', strin)
