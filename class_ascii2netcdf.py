#!/usr/bin/env python
'''
Base class for conversion of ascii data into netcdf format


History
-------
   * Written Dec 2025 by Matthias Cuntz, UMR Silva, INRAE
     from script ascii2isba and ascii2musica
   * Add ORCHIDEE, Matthias Cuntz, Nov 2026
   * Class variables for all standard variables, Matthias Cuntz, Nov 2026
   * Mode explicit names for important class variables,
     Matthias Cuntz, Nov 2026

'''
import numpy as np
import pandas as pd
import pyjams as pj


__all__ = ['ascii2Netcdf']


class ascii2Netcdf(object):
    """
    Class for making netCDF forcing file from ASCII data

    Parameters
    ----------
    verbose : bool, optional
        Report progress if True.

    """

    #
    # init
    #

    def __init__(self, verbose=False):

        self.verbose = verbose
        self.co2air = ''
        self.swdown = ''
        self.lwdown = ''
        self.psurf = ''
        self.qair = ''
        self.rainf = ''
        self.snowf = ''
        self.tair = ''
        self.wind1 = ''
        self.wind2 = ''
        self.hsbl = ''
        self.varnames = []
        self.varunits = []
        self.varlongnames = []
        self.varnames_lower = []
        self.dict_varunits = {}
        self.dict_varnames = {}
        self.dict_varlongnames = {}
        self.sep = None
        self.skip = None
        self.inputvars_lower = []
        self.inputvars = {}
        self.inputunits = {}
        self.dtcols = []
        self.dtstr = ''
        self.isdaily = False
        self.inputfile = ''
        self.interpolate = False
        self.df = pd.DataFrame()
        self.dt = -1
        self.outputfile = ""
        self.ftimestep = -1

        # Sets
        #     varnames, varnames_lower, varunits, varlongnames
        self.output_variables()

        return


    #
    # Output variables names and units
    #

    def output_variables(self):
        """
        Output variable names, longnames, and units

        """
        # names of variables
        self.co2air = 'CO2air'
        self.swdown = 'SWdown'
        self.lwdown = 'LWdown'
        self.psurf = 'PSurf'
        self.qair = 'Qair'
        self.rainf = 'Rainf'
        self.snowf = 'Snowf'
        self.tair = 'Tair'
        self.wind1 = 'Wind_N'
        self.wind2 = 'Wind_E'
        self.hsbl = 'h_sbl'
        # output variables and units
        self.stdnames = ['co2air', 'swdown', 'lwdown',
                         'psurf', 'qair', 'rainf',
                         'snowf', 'tair', 'wind1', 'wind2',
                         'hsbl']
        self.varnames = [self.co2air, self.swdown, self.lwdown,
                         self.psurf, self.qair, self.rainf,
                         self.snowf, self.tair, self.wind1, self.wind2,
                         self.hsbl]
        self.varunits = ['ppmv', 'W/m2', 'W/m2',
                         'Pa', 'kg/kg', 'kg/m2/s',
                         'kg/m2/s', 'K', 'm/s', 'm/s',
                         'm']
        self.varlongnames = ['Near surface CO2 concentration',
                             'Surface incident shortwave radiation',
                             'Surface incident longwave radiation',
                             'Surface pressure',
                             'Near surface specific humidity',
                             'Rainfall rate',
                             'Snowfall Rate',
                             'Near surface air temperature',
                             'Near surface northerly wind speed',
                             'Near surface easterly wind speed',
                             'Surface boundary layer height']

        # dictionary with lowercase variables as keys
        self.varnames_lower = [ vv.lower() for vv in self.varnames ]
        self.dict_varunits = dict(zip(self.varnames_lower, self.varunits))
        self.dict_varnames = dict(zip(self.varnames_lower, self.varnames))
        self.dict_varlongnames = dict(zip(self.varnames_lower,
                                          self.varlongnames))
        self.dict_stdnames = dict(zip(self.stdnames, self.varnames_lower))

    #
    # Analyse header
    #

    def get_header_field_separator(self, line):
        """
        Determine separator symbol ',' or ';'

        Parameters
        ----------
        line : string
            First line of input file

        Returns
        -------
        str
            field separator

        """
        sep = None
        if ',' in line:
            sep = ','
        elif ';' in line:
            sep = ';'

        return sep


    def get_header_index_tair(self, arr):
        """
        Index of column header starting with 'tair' (case insensitive)

        Parameters
        ----------
        arr : array_like
            List of header strings

        Returns
        -------
        int
            index of header starting with 'tair'.
            Returns -1 if no element starts with 'tair'

        """
        lline = [ nn.lower() for nn in arr ]

        iitair = -1
        for ii, nn in enumerate(lline):
            if nn.startswith('tair'):
                iitair = ii
                break

        return iitair


    def get_header_num_lines(self, arr):
        """
        Determine number of header lines

        Parameters
        ----------
        arr : array_like
            List of header strings

        Returns
        -------
        int
            Number of header lines: 1 or 2.

        """
        iitair = self.get_header_index_tair(arr)
        try:               # 1 header line if second line are numbers
            # num   = [ float(nn) for nn in line2arr ]
            skip = 1
            _ = float(arr[iitair])
        except ValueError: # 2 header lines
            skip = 2

        return skip


    def get_header_vars_units(self, arr):
        """
        Get variables names and units if in form var (unit)

        Parameters
        ----------
        arr : array_like
            List of header strings

        Returns
        -------
        list, list
            First variable names, second corresponding units

        """
        ncol = len(arr)
        inputvars = list()
        inputunits = list()
        for ic in range(ncol):
            vv = arr[ic]
            if ' ' in vv:
                inputvars.append(vv[:vv.find(' ')])
                iu = vv[vv.find(' '):].strip()
                if iu.startswith('(') or iu.startswith('['):
                    iu = iu[1:]
                if iu.endswith(')') or iu.endswith(']'):
                    iu = iu[:-1]
                inputunits.append(iu.strip())
            else:
                inputvars.append(vv.strip())
                inputunits.append('')

        return inputvars, inputunits


    def get_header_datetime(self, inputvars_lower):
        """
        Get datetime columns and format string

        Parameters
        ----------
        inputvars_lower : array_like
            List of lowercase variable names

        Returns
        -------
        list, str
            datetime columns, format string

        """
        iiyear = -1
        iimonth = -1
        iiday = -1
        iihour = -1
        iiminute = -1
        iisecond = -1
        iitime = -1
        iidoy = -1

        if 'year' in inputvars_lower:
            iiyear = inputvars_lower.index('year')
        else:
            raise IOError("'year' must be in header line.")

        if 'month' in inputvars_lower:
            iimonth = inputvars_lower.index('month')
        if 'day' in inputvars_lower:
            iiday = inputvars_lower.index('day')
        if 'jday' in inputvars_lower:
            iidoy = inputvars_lower.index('jday')
        if 'doy' in inputvars_lower:
            iidoy = inputvars_lower.index('doy')
        if 'hour' in inputvars_lower:
            iihour = inputvars_lower.index('hour')
        if 'minute' in inputvars_lower:
            iiminute = inputvars_lower.index('minute')
        if 'second' in inputvars_lower:
            iisecond = inputvars_lower.index('second')
        if 'time' in inputvars_lower:
            iitime = inputvars_lower.index('time')

        if (iimonth < 0) and (iidoy < 0):
            raise IOError("Could not determine date/time structure:"
                          " month given but day is missing. Date/time"
                          " structure must be"
                          "  'year,doy,time'  or"
                          "  'year,doy,hour,minute'  or"
                          "  'year,month,day,time'  or"
                          "  'year,month,day,hour,minute'."
                          " The order can be arbitrary.")

        dtcols = [iiyear]
        dtstr = '%Y'
        if iimonth >= 0:
            dtcols += [iimonth]
            dtstr += ' %m'
        if iiday >= 0:
            dtcols += [iiday]
            dtstr += ' %d'
        if iidoy >= 0:
            dtcols += [iidoy]
            dtstr += ' %j'
        if iihour >= 0:
            dtcols += [iihour]
            dtstr += ' %H'
        if iiminute >= 0:
            dtcols += [iiminute]
            dtstr += ' %M'
        if iisecond >= 0:
            dtcols += [iisecond]
            dtstr += ' %S'
        if iitime >= 0:
            dtcols += [iitime]
            # decimal time such as 13.25 are not an datetime code
            # make hour, minute, second later
            dtstr += ' TT'

        return dtcols, dtstr


    def check_header_mandatory_vars(self, inputvars_lower):
        """
        Check for mandatory variables

        Parameters
        ----------
        inputvars_lower : array_like
            List of lowercase variable names

        """
        sinputvars_lower = set(inputvars_lower)

        vlists = [['co2air', self.co2air.lower()],
                  ['lwdown', self.lwdown.lower()],
                  ['psurf', self.psurf.lower()],
                  ['tair', self.tair.lower()],
                  ['swdown', self.swdown.lower()],
                  ['rhair', 'qair', self.qair.lower()],
                  ['precip', 'rainf', 'snowf', self.rainf.lower(),
                   self.snowf.lower()],
                  ['windspeed', 'wind_speed', 'wind_n', 'wind_e',
                   'wind_dir', 'winddir',
                   self.wind1.lower(), self.wind2.lower()]]
        for vlist in vlists:
            # [x for x in a if x in b]
            if len(sinputvars_lower.intersection(vlist)) == 0:
                raise IOError(f'Mandatory variable not in header line:'
                              f' {vlist}')

        return


    def analyse_header(self, inputfile=''):
        """
        Analyse csv file

        detect
            field separator,
            date and time variables,
            variable names and units,
            precipitation and wind variables

        Parameters
        ----------
        inputfile : string, optional
            csv file with ISBA's forcing variables (default: self.inputfile)

        """
        if inputfile == '':
            inputfile = self.inputfile

        with open(inputfile, 'r') as fi:
            line1 = fi.readline().strip()
            line2 = fi.readline().strip()

        # Analyse field delimiter
        self.sep = self.get_header_field_separator(line1)

        line1arr = line1.split(self.sep)
        line2arr = line2.split(self.sep)

        # Analyse variables and units
        # columns index of CO2 variable
        iitair = self.get_header_index_tair(line1arr)
        if iitair < 0:
            raise ValueError("One column has to start with string"
                             " 'tair' (case insensitive).")
        # number of header lines
        self.skip = self.get_header_num_lines(line2arr)
        # variable names and units
        if self.skip == 1:
            ncol = len(line1arr)
            line11 = line1arr[iitair].split()
            if len(line11) == 1:  # 1 header line without units
                inputvars  = [ iv.strip() for iv in line1arr ]
                inputunits = [''] * ncol
            else:                # 1 header line with units
                inputvars, inputunits = self.get_header_vars_units(line1arr)
        else:                    # 2 header lines: 1. names, 2. units
            inputvars = [ iv.strip() for iv in line1arr ]
            inputunits = [ iu.strip() for iu in line2arr ]
        # lowercase variable names
        self.inputvars_lower = [ iv.lower() for iv in inputvars ]

        # vars and units dictionary
        self.inputvars = dict(zip(self.inputvars_lower, inputvars))
        self.inputunits = dict(zip(self.inputvars_lower, inputunits))

        # analyse datetime
        self.dtcols, self.dtstr = self.get_header_datetime(self.inputvars_lower)
        if ('%H' in self.dtstr) or ('TT' in self.dtstr):
            self.isdaily = False
        else:
            self.isdaily = True

        return


    def make_hour_minute_second_from_time(self, df):
        """
        Make hour, minute, and second column from time column,
        correct columns list and datetime string

        Parameters
        ----------
        df : pandas.DataFrame
            data
        dtcol : list
            Date and time indexes
        dtstr : str
            Datetime string ending with TT

        Returns
        -------
        df
            DataFrame with additional hour, minute, second columns.
            Changes self.dtcols and self.dtstr accordingly.

        """
        # hour, minute, second columns
        df['hour'] = np.floor(df['time']).astype(int)
        df['minute'] = np.round((df['time'] - df['hour']) * 60.).astype(int)
        df['second'] = np.round(((df['time'] - df['hour']) * 60.
                                 - df['minute']) * 60.).astype(int)
        # date and time columns
        cc = list(df.columns)
        self.dtcols = self.dtcols[:-1]
        self.dtcols.extend([cc.index('hour'), cc.index('minute'),
                            cc.index('second')])

        # datetime string
        self.dtstr = self.dtstr.replace('TT', '%H %M %S')

        return df


    #
    # Read data and datetime
    #

    def read_data(self, inputfile='', startdate='', enddate='',
                  interpolate=False):
        """
        Read data, get datetime, select columns, select time_span

        Parameters
        ----------
        inputfile : string, optional
            csv file with ISBA's forcing variables (default: self.inputfile)
        startdate : string, optional
            First possible date in netcdf output file in ISO8601 format.
            (Default: first date in input file)
        enddate : string, optional
            Last possible date in output netcdf file in ISO8601 format.
            (Default: last date in input file)
        interpolate : bool, optional
            Linearly interpolate missing data in input file given as NaN
            or empty cells.

        """
        if inputfile == '':
            inputfile = self.inputfile
        else:
            self.inputfile = inputfile
        self.interpolate = interpolate

        if self.verbose:
            print('Read csv file ', inputfile)

        # Sets
        #     sep, skip, inputvars, inputvars_lower, inputunits, dtcols, dtstr
        self.analyse_header(inputfile)

        # check for mandatory variables
        self.check_header_mandatory_vars(self.inputvars_lower)

        # read data
        skiprows = [1,] if self.skip > 1 else None
        df = pd.read_csv(inputfile, sep=self.sep, header=0,
                         names=self.inputvars_lower,
                         na_values = ['NaN', 'NA', 'nan', 'NAN'],
                         skiprows=skiprows, skipinitialspace=True)
        # convert time column 
        itime = False
        if self.dtstr.endswith('TT'):
            # decimal time such as 13.25 are not a datetime code
            itime = True
            df = self.make_hour_minute_second_from_time(df)
            self.inputvars_lower.extend(['hour', 'minute', 'second'])
        # make datetime index
        indx = pd.to_datetime(df.iloc[:, self.dtcols], format=self.dtstr)
        df.set_index(indx, inplace=True)
        df.index.name = 'Datetime'
        # delete date and time columns
        for vv in list(df.columns[self.dtcols]):
            self.inputvars_lower.remove(vv)
        df.drop(axis=0, columns=df.columns[self.dtcols], inplace=True)
        if itime:  # cannot do it before without messing up the indexes
            self.inputvars_lower.remove('time')
            df.drop(axis=0, columns=['time'], inplace=True)

        # interpolate
        if self.interpolate:
            df.interpolate(method='linear', axis=0, inplace=True)

        # time step
        ddtime = df.index.diff()
        dt = ddtime.min()
        if not all(ddtime[1:] == dt):
            raise ValueError(f'Not all timesteps equally distributed min/max: '
                             f' {ddtime.min()} / {ddtime.max()}')

        # Start and end dates        
        if startdate == '':
            startdate = df.index[0]
        else:
            startdate = pd.to_datetime(startdate, format='ISO8601')
        if enddate == '':
            enddate = df.index[-1]
        else:
            enddate = pd.to_datetime(enddate, format='ISO8601')
        df = df[(df.index >= startdate) & (df.index <= enddate)]

        self.df = df
        self.dt = dt.seconds

        return

    #
    # Convert variable units
    #

    def convert_units(self):
        """
        Convert variables to output units

        """
        for vv in self.inputvars_lower:
            if (self.inputunits[vv] == '') and (vv in self.varnames_lower):
                self.inputunits[vv] = self.dict_varunits[
                    self.dict_stdnames[vv]]

        vpd2qair = False
        rhair2qair = False
        for vv in self.inputvars_lower:
            if vv in ['precip', 'rainf', 'snowf']:
                if self.inputunits[vv] not in ['mm', 'mm/dt', 'mm dt-1',
                                               'kg/m2/dt', 'kg m-2 dt-1',
                                               'mm/s', 'mm s-1',
                                               'kg/m2/s','kg m-2 s-1']:
                    raise ValueError(f'Only units mm (kg/m2/dt) or mm/s'
                                     f' (kg/m2/s) possible for precipitation.'
                                     f' Given: {self.inputunits[vv]}')
                if self.inputunits[vv] in ['mm', 'mm/dt', 'mm dt-1',
                                           'kg/m2/dt', 'kg m-2 dt-1']:
                    # mm -> kg/m2/s
                    self.df[vv] /= float(self.dt)
            elif vv in ['wind', 'windspeed', 'wind_speed', 'wind_n', 'wind_e']:
                if self.inputunits[vv] not in ['m/s', 'm s-1', '']:
                    raise ValueError(f'Only units m/s possible for wind speed.'
                                     f' Given: {self.inputunits[vv]}')
                continue  # nothing to do
            elif vv in ['wind_dir', 'winddir']:
                if self.inputunits[vv] not in ['degree', 'deg', '°', '']:
                    raise ValueError(f'Only units degree possible for wind'
                                     f' direction. Given:'
                                     f' {self.inputunits[vv]}')
                continue  # nothing to do
            elif vv == 'co2air':
                if self.inputunits[vv] not in ['ppm', 'ppmv',
                                               'µmol/mol', 'µmol mol-1',
                                               'kg/m3', 'kg m-3']:
                    raise ValueError(f'Only units ppm, ppmv, µmol/mol,'
                                     f' and kg/m3 possible for CO2air.'
                                     f' Given: {self.inputunits[vv]}')
                continue  # done below
            elif vv == 'vpdair':
                if self.inputunits[vv] in ['kPa', 'hPa', 'Pa']:
                    # ICOS var == VPD
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    if self.inputunits[vv] == 'kPa':
                        self.df[vv] *= 1000.
                    elif self.inputunits[vv] == 'hPa':
                        self.df[vv] *= 100.
                vpd2qair = True
            elif vv == 'rhair':
                rhair2qair = True
            elif vv == 'qair':
                if self.inputunits[vv] in ['mol/mol', 'mol mol-1']:
                    # mol/mol -> kg/kg
                    self.df[vv] *= pj.const.molmass_h2o / pj.const.molmass_air
                elif self.inputunits[vv] in ['mmol/mol', 'mmol mol-1']:
                    # mmol/mol -> kg/kg
                    self.df[vv] *= (1e-3 * pj.const.molmass_h2o
                                      / pj.const.molmass_air)
                elif self.inputunits[vv] in ['kPa', 'hPa', 'Pa']:
                    # ICOS var == VPD
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    if self.inputunits[vv] == 'kPa':
                        self.df[vv] *= 1000.
                    elif self.inputunits[vv] == 'hPa':
                        self.df[vv] *= 100.
                    vpd2qair = True
                elif self.inputunits[vv] in ['%', 'percent', '', '0-1']:
                    # observed var == relative humidity
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    rhair2qair = True
                elif self.inputunits[vv] == self.dict_varunits[
                        self.dict_stdnames[vv]]:
                    continue  # kg/kg
                else:
                    raise ValueError(
                        f'Cannot convert from qair unit'
                        f' {self.inputunits[vv]} to'
                        f' {self.dict_varunits[self.dict_stdnames[vv]]}.')
            elif vv == 'tair':
                if self.inputunits[vv] in ['C', 'degreeC', 'degree C',
                                           'degC', 'deg C', '°C']:
                    # degC -> K
                    self.df[vv] += 273.15
                elif self.inputunits[vv] == self.dict_varunits[
                        self.dict_stdnames[vv]]:
                    continue  # K
                else:
                    raise ValueError(
                        f'Cannot convert from tair unit'
                        f' {self.inputunits[vv]} to'
                        f' {self.dict_varunits[self.dict_stdnames[vv]]}.')
            elif vv == 'psurf':
                if self.inputunits[vv] in ['hPa', 'mbar']:
                    # hPa -> Pa
                    self.df[vv] *= 100.
                elif self.inputunits[vv] == 'kPa':
                    # kPa -> Pa
                    self.df[vv] *= 1000.
                elif self.inputunits[vv] == 'bar':
                    # bar -> Pa
                    self.df[vv] *= 100000.
                elif (self.inputunits[vv] ==
                      self.dict_varunits[self.dict_stdnames[vv]]):
                    continue  # Pa
                else:
                    raise ValueError(
                        f'Cannot convert from psurf unit'
                        f' {self.inputunits[vv]} to'
                        f' {self.dict_varunits[self.dict_stdnames[vv]]}.')
            elif vv == 'h_sbl':
                if self.inputunits[vv] not in ['m', '']:
                    raise ValueError(f'Only units m possible for surface'
                                     f' boundary layer. Given: '
                                     f'{self.inputunits[vv]}')
                continue  # nothing to do
            elif vv == 'o3_air_ref':
                if (not (self.inputunits[vv] == 'ppt')):
                    raise ValueError(f'Only units ppt possible for ozone.'
                                     f' Given: {self.inputunits[vv]}')
                # ppt -> mol/mol
                self.df[vv] *= 1.e-9
            elif vv in self.dict_varunits:
                if self.inputunits[vv] == self.dict_varunits[vv]:
                    continue  # nothing to do
                elif self.inputunits[vv] in ['W/m2', 'W m-2']:
                    if self.dict_varunits[vv] not in ['W/m2', 'W m-2']:
                        raise ValueError(
                            f'Cannot convert from unit {self.inputunits[vv]}'
                            f' to {self.dict_varunits[vv]} for variable {vv}.')
                else:
                    raise ValueError(f'Cannot convert from unit'
                                     f' {self.inputunits[vv]} to'
                                     f' {self.dict_varunits[vv]} for variable'
                                     f' {vv}.')

        if vpd2qair:  # vpd2shair
            if 'vpdair' in self.inputvars_lower:
                vair = 'vpdair'
            else:
                vair = 'qair'
            vpd = self.df[vair]   # Pa
            tk = self.df['tair']  # K
            p = self.df['psurf']  # Pa
            eair = pj.vpd2eair(vpd, tk)
            self.df['qair'] = pj.eair2shair(eair, p)
            if vair == 'vpdair':
                self.df.drop(axis=0, columns=['vpdair'], inplace=True)
                self.inputvars_lower.remove('vpdair')

        if rhair2qair:  # rhair2shair
            if 'rhair' in self.inputvars_lower:
                vair = 'rhair'
            else:
                vair = 'qair'
            rh = self.df[vair]   # 0-1
            if np.any(rh > 2.):  # %
                rh *= 0.01
            tk = self.df['tair']
            p = self.df['psurf']
            eair = rh * pj.esat(tk)
            self.df['qair'] = pj.eair2shair(eair, p)
            if vair == 'rhair':
                self.df.drop(axis=0, columns=['rhair'], inplace=True)
                self.inputvars_lower.remove('rhair')

        # CO2
        if ((self.inputunits['co2air'] in ['kg/m3', 'kg m-3']) and
            (self.dict_varunits[self.dict_stdnames['co2air']] in
             ['ppm', 'ppmv', 'µmol/mol', 'µmol mol-1'])):
            raise ValueError('Conversion of CO2air from kg/m3 to ppm'
                             ' not implemented yet')

        if ((self.inputunits['co2air'] in
             ['ppm', 'ppmv', 'µmol/mol', 'µmol mol-1']) and
            (self.dict_varunits[self.dict_stdnames['co2air']] in
             ['kg/m3', 'kg m-3'])):
            tk = self.df['tair']
            pres = self.df['psurf']
            qair = self.df['qair']
            # density kg(air)/m3
            rhoair = ( pres /
                       ( pj.const.Rair * tk *
                         (1. +
                          (pj.const.molmass_h2o / pj.const.molmass_air - 1.) *
                          qair) +
                         pj.const.gravity * self.reference_height ) )
            # ppm -> kg/m3
            self.df['co2air'] = (self.df['co2air'] * 1e-6 *
                                 pj.const.molmass_co2 / pj.const.molmass_air *
                                 rhoair)
            self.inputunits['co2air'] = 'kg/m3'

        # Direct radiation >= 0
        for vv in ['swdown', self.swdown.lower()]:
            if vv in self.inputvars_lower:
                self.df[vv] = np.maximum(self.df[vv], 0.)

        # Precip in rain and snow
        if 'precip' in self.inputvars_lower:
            tk = self.df['tair']
            precip = self.df['precip']
            self.df['rainf'] = np.where(tk > 274.15, precip, 0.)  # rain > 1 C
            self.df['snowf'] = np.where(tk > 274.15, 0., precip)  # snow < 1 C
            self.df.drop(axis=0, columns=['precip'], inplace=True)
            self.inputvars_lower.remove('precip')

        # daily input from Dietrich et al. (Annals of Forest Science, 2019)
        for dvar in ['grhds', 'grids', 'rrds', 'sddm',
                     'tadm', 'tadn', 'tadx', 'wsdm']:
            if dvar in self.inputvars_lower:
                # Input was given as value*100 as integer -> 2 digit precision
                self.df[dvar] *= 0.01

        return
