#!/usr/bin/env python
"""
Write ORCHIDEE netcdf forcing file from csv input file.

# Script can be called on the command line:
python ascii2orchidee.py -s 2005-01-01T00:30 -e 2006-01-01T00:30 -o ORCHIDEE_site_in_FR-Hes_2005-01.nc -t 1.0 --latitude=48.6742166667 --longitude=7.06461666667 --time2gmt=0 --forcing_height=22.5 -f=-9999 ORCHIDEE_site_in_FR-Hes_2005.csv

# The wrapper function can be called from within Python:
from ascii2orchidee import ascii2orchidee
ascii2orchidee(csvfile,
               latitude=latitude, longitude=longitude,
               forcing_height=forcing_height, time2gmt=time2gmt,
               startdate=startdate, enddate=enddate, interpolate=interpolate,
               outputfile=outputfile, fill_value=fill_value, ftimestep=ftimestep,
               verbose=False)

# The class can be used from within Python:
from ascii2orchidee import ascii2Orchidee
orchidee = ascii2Orchidee(latitude=latitude,
                          longitude=longitude,
                          forcing_height=forcing_height,
                          time2gmt=time2gmt,
                          verbose=False)
orchidee.read_data(csvfile,
                   startdate=startdate, enddate=enddate,
                   interpolate=interpolate)
orchidee.convert_units_orchidee()
orchidee.write_netcdf(outputfile, fill_value=fill_value,
                      ftimestep=ftimestep)


History
-------
   * Written Nov 2026 by Matthias Cuntz @ INRAE
     from ascii2musica.py

"""
import time as ptime
import cftime as cf
import netCDF4 as nc
import numpy as np
from class_ascii2netcdf import ascii2Netcdf


__all__ = ['ascii2Orchidee', 'ascii2orchidee', '__docstring__']


__docstring__ = """
Write ORCHIDEE netcdf forcing file from csv input file.

The script predicts the csv file structure from the first two lines, i.e.
variable names and units, field delimiter and date/time variables.

Apart from date/time information
(year,jday,time or year,month,day,hour,minute),
the script makes one netcdf variable for each column found in the csv file.

Additional information such as forcing height, site latitude, etc. can be given
as optional inputs.

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
However, ORCHIDEE needs the following forcing data (all at reference level) and
so these variables have to be present if sub-diurnal data:
       swdown | Incoming solar radiation   | W/m^2
       lwdown | Incoming thermal radiation | W/m^2
       rainf  | Liquid precipitation       | kg/m^2s
       snowf  | Solid precipitation        | kg/mv2s
         tair | Air temperature            | K
        psurf | Atmospheric pressure       | Pa
       co2air | CO2 mole fraction          | umolCO2 mol-1
         qair | air specific humidity      | kg/kg
   wind_speed | Wind speed                 | m/s
     wind_dir | Wind direction             | degreeN
              or
       wind_n | Northerly wind speed       | m/s
       wind_e | Easterly wind speed        | m/s
If other units are given, the script tries to convert them.
If no units are given, the scripts assumes that the variables have already
the right units.
If precip is given, rain and snow will be split at 1 degree Celsius
air temperature.

Additional Information
----------------------
Additional information about the station must be known. These are:
  - forcing_height: the reference height where meteorological variables and
    fluxes were measured (in meters)
  - latitude: the site latitude (in degrees, -90 to +90)
  - longitude: the site longitude (in degrees, -180 to +180)
  - time2gmt: the time shift (in hours) between date/time in the csv file and
    GMT=UTC (from -12 to +12)
These can be given on the command line by their name preceded by two dashes,
e.g. --latitude 17.52 or --latitude=17.52.

"""


class ascii2Orchidee(ascii2Netcdf):
    """
    Class for making ORCHIDEE netCDF forcing file from ASCII data

    Parameters
    ----------
    latitude : float
        The site's latitude from -90 to 90.
    longitude : float
        The site's longitude from -180 to 180.
    forcing_height : float, optional
        The height were temperature and humidity were measured in m
        (default: 2.).
    time2gmt : float, optional
        The time shift (in hours) between date/time in the csv file and
        GMT=UTC (from -12 to +12, default: 0).

    """

    def __init__(self, *args,
                 latitude=None, longitude=None,
                 forcing_height=2., time2gmt=0.,
                 **kwargs):
        super().__init__(*args, **kwargs)
        
        if latitude is None:
            raise ValueError("Latitude must be given.")
        else:
            self.latitude = latitude

        if longitude is None:
            raise ValueError("Longitude must be given.")
        else:
            self.longitude = longitude

        self.forcing_height = forcing_height

        self.time2gmt = time2gmt

        return

    #
    # Output variables names and units (replace standard)
    #

    def output_variables(self):
        """
        Output variable names, longnames, and units

        """
        # names of standard variables
        self.co2air = 'CO2_f'
        self.swdown = 'SWdown'
        self.lwdown = 'LWdown'
        self.psurf = 'PSurf'
        self.qair = 'Qair'
        self.rainf = 'Rainf'
        self.snowf = 'Snowf'
        self.tair = 'Tair'
        self.wind1 = 'Wind_N'
        self.wind2 = 'Wind_E'
        # output variables and units
        self.stdnames = ['co2air', 'swdown', 'lwdown',
                         'psurf', 'qair', 'rainf',
                         'snowf', 'tair', 'wind1', 'wind2']
        self.varnames = [self.co2air, self.swdown, self.lwdown,
                         self.psurf, self.qair, self.rainf,
                         self.snowf, self.tair, self.wind1, self.wind2]
        self.varunits = ['ppm', 'W/m2', 'W/m2', 'Pa',
                         'kg/kg', 'kg/m2/s', 'kg/m2/s', 'K',
                         'm/s', 'm/s']
        self.varlongnames = ['CO2 mole fraction',
                             'Surface incident shortwave radiation',
                             'Surface incident longwave radiation',
                             'Surface pressure',
                             'Near surface specific humidity',
                             'Rainfall rate',
                             'Snowfall rate',
                             'Near surface air temperature',
                             'Near surface northward wind component',
                             'Near surface eastward wind component']
        
        # dictionary with lowercase variables as keys
        self.varnames_lower = [ vv.lower() for vv in self.varnames ]
        self.dict_varunits = dict(zip(self.varnames_lower, self.varunits))
        self.dict_varnames = dict(zip(self.varnames_lower, self.varnames))
        self.dict_varlongnames = dict(zip(self.varnames_lower,
                                          self.varlongnames))
        self.dict_stdnames = dict(zip(self.stdnames, self.varnames_lower))

    #
    # Convert variable units
    #

    def convert_units_orchidee(self):
        """
        Convert variables to output units

        """
        # standard conversions
        self.convert_units()

        # rename variables
        if 'co2air' in self.inputvars_lower:
            self.inputvars_lower[self.inputvars_lower.index('co2air')] = (
                self.co2air.lower())
            self.inputvars[self.co2air] = self.inputvars.pop('co2air')
        if 'swdown' in self.inputvars_lower:
            self.inputvars_lower[self.inputvars_lower.index('swdown')] = (
                self.swdown.lower())
            self.inputvars[self.swdown] = self.inputvars.pop('swdown')
        if 'lwdown' in self.inputvars_lower:
            self.inputvars_lower[self.inputvars_lower.index('lwdown')] = (
                self.lwdown.lower())
            self.inputvars[self.lwdown] = self.inputvars.pop('lwdown')
        if 'psurf' in self.inputvars_lower:
            self.inputvars_lower[self.inputvars_lower.index('psurf')] = (
                self.psurf.lower())
            self.inputvars[self.psurf] = self.inputvars.pop('psurf')
        if 'qair' in self.inputvars_lower:
            self.inputvars_lower[self.inputvars_lower.index('qair')] = (
                self.qair.lower())
            self.inputvars[self.qair] = self.inputvars.pop('qair')
        if 'tair' in self.inputvars_lower:
            self.inputvars_lower[self.inputvars_lower.index('tair')] = (
                self.tair.lower())
            self.inputvars[self.tair] = self.inputvars.pop('tair')

        dren = {}
        if 'co2air' in self.df:
            dren.update({'co2air': self.co2air.lower()})
        if 'swdown' in self.df:
            dren.update({'swdown': self.swdown.lower()})
        if 'lwdown' in self.df:
            dren.update({'lwdown': self.lwdown.lower()})
        if 'psurf' in self.df:
            dren.update({'psurf': self.psurf.lower()})
        if 'qair' in self.df:
            dren.update({'qair': self.qair.lower()})
        if 'tair' in self.df:
            dren.update({'tair': self.tair.lower()})
        if len(dren) > 0:
            self.df = self.df.rename(columns=dren)

        # wind_n / wind_e
        dren = {}
        if self.wind1.lower() in self.inputvars_lower:  # wind_n
            if self.wind2.lower() not in self.inputvars_lower:  # wind_e
                self.df[self.wind2.lower()] = np.zeros_like(
                    self.df[self.wind1.lower()])
        else:
            for vv in ['wind', 'windspeed', 'wind_speed']:
                if vv in self.inputvars_lower:
                    vwind = vv
                    break
            vdir = ''
            for vv in ['wind_dir', 'winddir']:
                if vv in self.inputvars_lower:
                    vdir = vv
                    break
            if vdir:
                # northerly / easterly
                self.df[self.wind1.lower()] = (self.df[vwind] *
                                               np.sin(np.deg2rad(
                                                   self.df[vdir])))
                self.df[self.wind2.lower()] = (self.df[vwind] *
                                               np.cos(np.deg2rad(
                                                   self.df[vdir])))
                self.df.drop(axis=0, columns=[vdir], inplace=True)
                self.inputvars_lower.remove(vdir)
            else:
                self.df[self.wind1.lower()] = self.df[vwind]  # northerly
                self.df[self.wind2.lower()] = np.zeros_like(
                    self.df[self.wind2.lower()])  # easterly
            if vwind not in [self.wind1.lower(), self.wind2.lower()]:
                self.df.drop(axis=0, columns=[vwind], inplace=True)
                self.inputvars_lower.remove(vwind)

        return

    #
    # Write output file
    #

    def write_netcdf(self, outputfile='', fill_value=None, ftimestep=1.0):
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

            Note that ORCHIDEE is using the middle of the time step, so
            output times will be transformed to the middle of the
            timestep. (Default: 1.0)

        """        
        if outputfile == '':
           outputfile = self.inputfile[0:self.inputfile.rfind('.')] + '.nc'
        self.outputfile = outputfile
        self.ftimestep = ftimestep

        # shift time to middle of time step
        if self.ftimestep != 0.0:
            dt = np.timedelta64(self.dt, 's')
            indx = self.df.index + (0.0 - ftimestep) * dt
            self.df.set_index(indx, inplace=True)

        # NETCDF
        if self.verbose:
            print('Create netcdf file ', outputfile)
        with nc.Dataset(outputfile, 'w', format='NETCDF3_64BIT') as fo:

            # Structure
            lon  = fo.createDimension('lon', 1)
            lat  = fo.createDimension('lat', 1)
            lev  = fo.createDimension('lev', 1)
            tstep = fo.createDimension('tstep', None)
            bound = fo.createDimension('bound', 2)

            vlon           = fo.createVariable('lon', 'f8', ('lon',))
            vlon.units     = 'degrees_east'
            vlon.long_name = 'Longitude'

            vlat           = fo.createVariable('lat', 'f8', ('lat',))
            vlat.units     = 'degrees_north'
            vlat.long_name = 'Latitude'

            vlev           = fo.createVariable('lev', 'f8', ('lev',))
            vlev.units     = 'm'
            vlev.long_name = 'Vertical Levels'

            vtstep          = fo.createVariable('tstep', 'f8', ('tstep',))
            vtstep.units    = (f'seconds since'
                               f' {self.df.index[0].year}-01-01 00:00:00')
            vtstep.calendar = 'gregorian'
            vtstep.origin   = f'{self.df.index[0].year}-01-01 00:00:00'

            blon = fo.createVariable('bounds_lon', 'f4', ('lon', 'bound'))
            blat = fo.createVariable('bounds_lat', 'f4', ('lat', 'bound'))

            vloc = dict()
            # Write all columns not only the ORCHIDEE forcing columns
            # for iv in self.dict_varnames:
            for iv in self.df.columns:
                if iv in ['tair', 'qair', 'wind_n', 'wind_e']:
                    vloc[iv] = fo.createVariable(
                        self.dict_varnames[iv], 'f8',
                        ('tstep', 'lev', 'lat', 'lon',),
                        fill_value=fill_value)
                else:
                    if iv in self.dict_varnames:
                        vloc[iv] = fo.createVariable(
                            self.dict_varnames[iv], 'f8',
                            ('tstep', 'lat', 'lon',),
                            fill_value=fill_value)
                    elif iv in self.inputvars_lower:
                        vloc[iv] = fo.createVariable(
                            self.inputvars[iv], 'f8',
                            ('tstep', 'lat', 'lon',),
                            fill_value=fill_value)
                    else:
                        vloc[iv] = fo.createVariable(
                            iv, 'f8',
                            ('tstep', 'lat', 'lon',),
                            fill_value=fill_value)
                if iv in self.dict_varnames:
                    if self.dict_varunits[iv]:
                        vloc[iv].units = self.dict_varunits[iv]
                    if self.dict_varlongnames[iv]:
                        vloc[iv].long_name = self.dict_varlongnames[iv]
                elif iv in self.inputvars_lower:
                    if self.inputunits[iv]:
                        vloc[iv].units = self.inputunits[iv]

            fo.description = f'Half-hourly data'

            # calc netCDF time
            cdates = [ cf.datetime(dd.year, dd.month, dd.day, dd.hour,
                                   dd.minute, dd.second)
                       for dd in self.df.index ]
            dtime = cf.date2num(
                cdates,
                f'seconds since {self.df.index[0].year}-01-01 00:00:00',
                calendar='gregorian')
            # if time2gmt given, shift to UTC
            if self.time2gmt != 0:
                dtime = dtime + self.time2gmt * 3600.

            # Write variables
            vlon[:] = self.longitude
            vlat[:] = self.latitude
            vlev[:] = self.forcing_height
            vtstep[:] = dtime
            blon[:] = [self.longitude - 0.5, self.longitude + 0.5]
            blat[:] = [self.latitude - 0.5, self.latitude + 0.5]

            # for iv in self.dict_varnames:
            for iv in self.df.columns:
                vloc[iv][:] = self.df[iv]

        return


# -------------------------------------------------------------------------
# Function using class
#

def ascii2orchidee(csvfile,
                   latitude=None, longitude=None,
                   forcing_height=None, time2gmt=None,
                   startdate='', enddate='', interpolate=False,
                   outputfile='', fill_value=None, ftimestep=1.0,
                   verbose=False, **kwargs):
    """
    Write ORCHIDEE netcdf forcing file from csv input file.

    See variable `__docstring__` or use --help on command line.

    Parameters
    ----------
    csvfile : string
        csv file with ORCHIDEE's forcing variables
    latitude : float
        The site's latitude from -90 to 90.
    longitude : float
        The site's longitude from -180 to 180.
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

        Note that ORCHIDEE is using the middle of the time step, so
        output times will be transformed to the middle of the
        timestep. (Default: 1.0)
    verbose : bool, optional
        Report progress if True.
    **kwargs : dict, optional
        All other keyword arguments will be ignored

    """
    orchidee = ascii2Orchidee(latitude=latitude,
                              longitude=longitude,
                              forcing_height=forcing_height,
                              time2gmt=time2gmt,
                              verbose=verbose)
    orchidee.read_data(csvfile,
                       startdate=startdate, enddate=enddate,
                       interpolate=interpolate)
    orchidee.convert_units_orchidee()
    orchidee.write_netcdf(outputfile, fill_value=fill_value,
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
    ftimestep = 1.0
    latitude = None
    longitude = None
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
    hstr = ('Fraction of time step at which time is given (default: 1.0):'
            ' 0=beginning of time step, 0.5=middle of time step, 1=end of'
            ' time step. Note that ORCHIDEE is using the middle of the time'
            ' step so output times will be transformed to the middle of the'
            ' timestep.')
    parser.add_argument('-t', '--ftimestep', action='store', type=float,
                        default=ftimestep, dest='ftimestep',
                        metavar='fractional time step', help=hstr)
    parser.add_argument('--latitude', action='store', type=float,
                        default=latitude, dest='latitude',
                        metavar='latitude',
                        help="The site's latitude from -90 to 90.")
    parser.add_argument('--longitude', action='store', type=float,
                        default=longitude, dest='longitude',
                        metavar='longitude',
                        help="The site's longitude from -180 to 180.")
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
                        help="csv file with ORCHIDEE's forcing variables.")

    args = parser.parse_args()

    assert args.infile is not None, (
        'csv file with ORCHIDEE forcing variables must be given.')

    t1 = ptime.time()

    orchidee = ascii2Orchidee(latitude=args.latitude,
                              longitude=args.longitude,
                              forcing_height=args.forcing_height,
                              time2gmt=args.time2gmt,
                              verbose=True)
    orchidee.read_data(args.infile,
                       startdate=args.startdate, enddate=args.enddate,
                       interpolate=args.interpolate)
    orchidee.convert_units_orchidee()
    orchidee.write_netcdf(args.outputfile, fill_value=args.fill_value,
                          ftimestep=args.ftimestep)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.)
             if (t2 - t1) > 60. else '[s]: {:d}'.format(int(t2 - t1)))
    print('Time elapsed: ', strin)
