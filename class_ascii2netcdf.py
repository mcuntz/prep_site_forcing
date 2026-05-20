#!/usr/bin/env python
'''
Base class for conversion of ascii data into netcdf format


History
-------
   * Written Dec 2025 by Matthias Cuntz, UMR Silva, INRAE
     from script ascii2isba and ascii2musica

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
        self.swdown = ''
        self.lswdown = ''
        self.windspeed = ''
        self.lwindspeed = ''
        self.varnames = []
        self.varlownames = []
        self.varunits = []
        self.varlongnames = []
        self.dvarunits = {}
        self.dvarnames = {}
        self.dvarlongnames = {}
        self.sep = None
        self.skip = None
        self.livars = []
        self.ivars = {}
        self.iunits = {}
        self.dtcols = []
        self.dtstr = ''
        self.isdaily = False
        self.ifile = ''
        self.interpolate = False
        self.df = pd.DataFrame()
        self.dt = -1
        self.outputfile = ""
        self.ftimestep = -1

        # Sets
        #     varnames, varlownames, varunits, varlongnames
        self.output_variables()

        return


    #
    # Output variables names and units
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
        self.varnames = ['CO2air', 'SWdown', 'LWdown', 'PSurf', 'Qair', 'Rainf',
                         'Snowf', 'Tair', 'Wind_N', 'Wind_E', 'h_sbl']
        self.varlownames = [ vv.lower() for vv in self.varnames ]
        self.varunits = ['ppmv', 'W/m2', 'W/m2', 'Pa', 'kg/kg', 'kg/m2/s',
                         'kg/m2/s', 'K', 'm/s', 'm/s', 'm']
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
        self.dvarunits = dict(zip(self.varlownames, self.varunits))
        self.dvarnames = dict(zip(self.varlownames, self.varnames))
        self.dvarlongnames = dict(zip(self.varlownames, self.varlongnames))

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
        ivars = list()
        iunits = list()
        for ic in range(ncol):
            if ' ' in arr[ic]:
                iv, iu = arr[ic].split()
                ivars.append(iv.strip())
                iu = iu.strip()
                if iu.startswith('(') or iu.startswith('['):
                    iu = iu[1:]
                if iu.endswith(')') or iu.endswith(']'):
                    iu = iu[:-1]
                    iunits.append(iu.strip())
            else:
                ivars.append(arr[ic].strip())
                iunits.append('')

        return ivars, iunits


    def get_header_datetime(self, livars):
        """
        Get datetime columns and format string

        Parameters
        ----------
        livars : array_like
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

        if 'year' in livars:
            iiyear = livars.index('year')
        else:
            raise IOError("'year' must be in header line.")

        if 'month' in livars:
            iimonth = livars.index('month')
        if 'day' in livars:
            iiday = livars.index('day')
        if 'jday' in livars:
            iidoy = livars.index('jday')
        if 'doy' in livars:
            iidoy = livars.index('doy')
        if 'hour' in livars:
            iihour = livars.index('hour')
        if 'minute' in livars:
            iiminute = livars.index('minute')
        if 'second' in livars:
            iisecond = livars.index('second')
        if 'time' in livars:
            iitime = livars.index('time')

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


    def check_header_mandatory_vars(self, livars):
        """
        Check for mandatory variables

        Parameters
        ----------
        livars : array_like
            List of lowercase variable names

        """
        vlist = ['co2air', 'lwdown', 'psurf', 'tair']
        for iv in vlist:
            if iv.lower() not in livars:  # case in-sensitive
                raise IOError(f'Mandatory variable {iv} not in header line.')

        hasswdown = False
        vlist = ['swdown', self.lswdown]
        for iv in vlist:
            if iv.lower() in livars:
                hasswdown = True
        if not hasswdown:
            raise IOError(f'Mandatory variable swdown or'
                          f' {self.lswdown} not in header line.')

        hashum = False
        vlist = ['rhair', 'qair']
        for iv in vlist:
            if iv.lower() in livars:
                hashum = True
        if not hashum:
            raise IOError('Mandatory variable rhair or qair not in header line.')

        hasprecip = False
        vlist = ['precip', 'rainf', 'snowf']
        for iv in vlist:
            if iv.lower() in livars:
                hasprecip = True
        if not hasprecip:
            raise IOError('Mandatory variable precip, or rainf and snowf'
                          ' not in header line.')
        else:
            if 'precip' not in livars:
                if ('rainf' not in livars) or ('snowf' not in livars):
                    raise IOError('Both mandatory variables rainf and snowf'
                                  ' must be given. Otherwise provide total'
                                  ' precip.')

        hasspeed = False
        vlist = ['windspeed', 'wind_speed', 'wind_n', 'wind_e', self.lwindspeed]
        for iv in vlist:
            if iv.lower() in livars:
                hasspeed = True
        if not hasspeed:
            raise IOError('Mandatory variable wind, windspeed, wind_speed, wind_n,'
                          ' or wind_e not in header line.')

        return


    def analyse_header(self, ifile=''):
        """
        Analyse csv file

        detect
            field separator,
            date and time variables,
            variable names and units,
            precipitation and wind variables

        Parameters
        ----------
        ifile : string, optional
            csv file with ISBA's forcing variables (default: self.ifile)

        """
        if ifile == '':
            ifile = self.ifile

        with open(ifile, 'r') as fi:
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
                ivars  = [ iv.strip() for iv in line1arr ]
                iunits = [''] * ncol
            else:                # 1 header line with units
                ivars, iunits = self.get_header_vars_units(line1arr)
        else:                    # 2 header lines: 1. names, 2. units
            ivars = [ iv.strip() for iv in line1arr ]
            iunits = [ iu.strip() for iu in line2arr ]
        # lowercase variable names
        self.livars = [ iv.lower() for iv in ivars ]

        # vars and units dictionary
        self.ivars = dict(zip(self.livars, ivars))
        self.iunits = dict(zip(self.livars, iunits))

        # analyse datetime
        self.dtcols, self.dtstr = self.get_header_datetime(self.livars)
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

    def read_data(self, ifile='', startdate='', enddate='', interpolate=False):
        """
        Read data, get datetime, select columns, select time_span

        Parameters
        ----------
        ifile : string, optional
            csv file with ISBA's forcing variables (default: self.ifile)
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
        if ifile == '':
            ifile = self.ifile
        else:
            self.ifile = ifile
        self.interpolate = interpolate

        if self.verbose:
            print('Read csv file ', ifile)

        # Sets
        #     sep, skip, ivars, livars, iunits, dtcols, dtstr
        self.analyse_header(ifile)

        # check for mandatory variables
        self.check_header_mandatory_vars(self.livars)

        # read data
        skiprows = [1,] if self.skip > 1 else None
        df = pd.read_csv(ifile, sep=self.sep, header=0, names=self.livars,
                         na_values = ['NaN', 'NA', 'nan', 'NAN'],
                         skiprows=skiprows, skipinitialspace=True)
        # convert time column 
        itime = False
        if self.dtstr.endswith('TT'):
            # decimal time such as 13.25 are not a datetime code
            itime = True
            df = self.make_hour_minute_second_from_time(df)
            self.livars.extend(['hour', 'minute', 'second'])
        # make datetime index
        indx = pd.to_datetime(df.iloc[:, self.dtcols], format=self.dtstr)
        df.set_index(indx, inplace=True)
        df.index.name = 'Datetime'
        # delete date and time columns
        for vv in list(df.columns[self.dtcols]):
            self.livars.remove(vv)
        df.drop(axis=0, columns=df.columns[self.dtcols], inplace=True)
        if itime:  # cannot do it before without messing up the indexes
            self.livars.remove('time')
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
        for iv in self.livars:
            if (self.iunits[iv] == '') and (iv in self.varlownames):
                self.iunits[iv] = self.dvarunits[iv]

        vpd2qair = False
        rhair2qair = False
        for vv in self.livars:
            if vv in ['precip', 'rainf', 'snowf']:
                if self.iunits[vv] not in ['mm', 'mm/dt', 'mm dt-1',
                                           'kg/m2/dt', 'kg m-2 dt-1',
                                           'mm/s', 'mm s-1',
                                           'kg/m2/s','kg m-2 s-1']:
                    raise ValueError(f'Only units mm (kg/m2/dt) or mm/s'
                                     f' (kg/m2/s) possible for precipitation.'
                                     f' Given: {self.iunits[vv]}')
                if self.iunits[vv] in ['mm', 'mm/dt', 'mm dt-1',
                                       'kg/m2/dt', 'kg m-2 dt-1']:
                    # mm -> kg/m2/s
                    self.df[vv] /= float(self.dt)
            elif vv in ['wind', 'windspeed', 'wind_speed', 'wind_n', 'wind_e',
                        self.lwindspeed]:
                if self.iunits[vv] not in ['m/s', 'm s-1', '']:
                    raise ValueError(f'Only units m/s possible for wind speed.'
                                     f' Given: {self.iunits[vv]}')
                continue  # nothing to do
            elif vv in ['wind_dir', 'winddir']:
                if self.iunits[vv] not in ['degree', 'deg', '°', '']:
                    raise ValueError(f'Only units degree possible for wind'
                                     f' direction. Given: {self.iunits[vv]}')
                continue  # nothing to do
            elif vv == 'co2air':
                if self.iunits[vv] not in ['ppm', 'ppmv',
                                           'µmol/mol', 'µmol mol-1',
                                           'kg/m3', 'kg m-3']:
                    raise ValueError(f'Only units ppm, ppmv, µmol/mol,'
                                     f' and kg/m3 possible for CO2air.'
                                     f' Given: {self.iunits[vv]}')
                continue  # done below
            elif vv == 'vpdair':
                if self.iunits[vv] in ['kPa', 'hPa', 'Pa']:
                    # ICOS var == VPD
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    if self.iunits[vv] == 'kPa':
                        self.df[vv] *= 1000.
                    elif self.iunits[vv] == 'hPa':
                        self.df[vv] *= 100.
                vpd2qair = True
            elif vv == 'rhair':
                rhair2qair = True
            elif vv == 'qair':
                if self.iunits[vv] in ['mol/mol', 'mol mol-1']:
                    # mol/mol -> kg/kg
                    self.df[vv] *= pj.const.molmass_h2o / pj.const.molmass_air
                elif self.iunits[vv] in ['mmol/mol', 'mmol mol-1']:
                    # mmol/mol -> kg/kg
                    self.df[vv] *= (1e-3 * pj.const.molmass_h2o
                                      / pj.const.molmass_air)
                elif self.iunits[vv] in ['kPa', 'hPa', 'Pa']:
                    # ICOS var == VPD
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    if self.iunits[vv] == 'kPa':
                        self.df[vv] *= 1000.
                    elif self.iunits[vv] == 'hPa':
                        self.df[vv] *= 100.
                    vpd2qair = True
                elif self.iunits[vv] in ['%', 'percent', '', '0-1']:
                    # observed var == relative humidity
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    rhair2qair = True
                elif self.iunits[vv] == self.dvarunits[vv]:
                    continue  # kg/kg
                else:
                    raise ValueError(f'Cannot convert from qair unit'
                                     f' {self.iunits[vv]} to'
                                     f' {self.dvarunits[vv]}.')
            elif vv == 'tair':
                if self.iunits[vv] in ['C', 'degreeC', 'degree C',
                                       'degC', 'deg C', '°C']:
                    # degC -> K
                    self.df[vv] += 273.15
                elif self.iunits[vv] == self.dvarunits[vv]:
                    continue  # K
                else:
                    raise ValueError(f'Cannot convert from tair unit'
                                     f' {self.iunits[vv]} to'
                                     f' {self.dvarunits[vv]}.')
            elif vv == 'psurf':
                if self.iunits[vv] in ['hPa', 'mbar']:
                    # hPa -> Pa
                    self.df[vv] *= 100.
                elif self.iunits[vv] == 'kPa':
                    # kPa -> Pa
                    self.df[vv] *= 1000.
                elif self.iunits[vv] == 'bar':
                    # bar -> Pa
                    self.df[vv] *= 100000.
                elif self.iunits[vv] == self.dvarunits[vv]:
                    continue  # Pa
                else:
                    raise ValueError(f'Cannot convert from psurf unit'
                                     f' {self.iunits[vv]} to'
                                     f' {self.dvarunits[vv]}.')
            elif vv == 'h_sbl':
                if self.iunits[vv] not in ['m', '']:
                    raise ValueError(f'Only units m possible for surface'
                                     f' boundary layer. Given: '
                                     f'{self.iunits[vv]}')
                continue  # nothing to do
            elif vv == 'o3_air_ref':
                if (not (self.iunits[vv] == 'ppt')):
                    raise ValueError(f'Only units ppt possible for ozone.'
                                     f' Given: {self.iunits[vv]}')
                # ppt -> mol/mol
                self.df[vv] *= 1.e-9
            elif vv in self.dvarunits:
                if self.iunits[vv] == self.dvarunits[vv]:
                    continue  # nothing to do
                elif self.iunits[vv] in ['W/m2', 'W m-2']:
                    if self.dvarunits[vv] not in ['W/m2', 'W m-2']:
                        raise ValueError(
                            f'Cannot convert from unit {self.iunits[vv]}'
                            f' to {self.dvarunits[vv]} for variable {vv}.')
                else:
                    raise ValueError(f'Cannot convert from unit'
                                     f' {self.iunits[vv]} to'
                                     f' {self.dvarunits[vv]} for variable'
                                     f' {vv}.')

        if vpd2qair:  # vpd2shair
            if 'vpdair' in self.livars:
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
                self.livars.remove('vpdair')

        if rhair2qair:  # rhair2shair
            if 'rhair' in self.livars:
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
                self.livars.remove('rhair')

        # CO2
        if ((self.iunits['co2air'] in ['kg/m3', 'kg m-3']) and
            (self.dvarunits['co2air'] in
             ['ppm', 'ppmv', 'µmol/mol', 'µmol mol-1'])):
            raise ValueError('Conversion of CO2air from kg/m3 to ppm'
                             ' not implemented yet')

        if ((self.iunits['co2air'] in
             ['ppm', 'ppmv', 'µmol/mol', 'µmol mol-1']) and
            (self.dvarunits['co2air'] in ['kg/m3', 'kg m-3'])):
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
            self.iunits['co2air'] = 'kg/m3'

        # Direct radiation >= 0
        if 'swdown' in self.livars:
            self.df['swdown'] = np.maximum(self.df['swdown'], 0.)
        if self.lswdown in self.livars:
            self.df[self.lswdown] = np.maximum(self.df[self.lswdown], 0.)

        # Precip in rain and snow
        if 'precip' in self.livars:
            tk = self.df['tair']
            precip = self.df['precip']
            self.df['rainf'] = np.where(tk > 274.15, precip, 0.)  # rain > 1 degC
            self.df['snowf'] = np.where(tk > 274.15, 0., precip)  # snow < 1 degC
            self.df.drop(axis=0, columns=['precip'], inplace=True)
            self.livars.remove('precip')

        # daily input from Dietrich et al. (Annals of Forest Science, 2019)
        for dvar in ['grhds', 'grids', 'rrds', 'sddm',
                     'tadm', 'tadn', 'tadx', 'wsdm']:
            if dvar in self.livars:
                # Input was given as value*100 as integer -> 2 digit precision
                self.df[dvar] *= 0.01

        return
