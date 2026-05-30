#!/usr/bin/env python
'''
Class for preparation of meteo data to run ecosystem models

History
-------
   * Written Dec 2025 by Matthias Cuntz, UMR Silva, INRAE
     from script prep_site_forcing.py
   * grib, netcdf, and csv ERA5 files possible, Matthias Cuntz, Apr 2026
   * Allow no input file to get all ERA5 forcing, Matthias Cuntz, Apr 2026

'''
from collections.abc import Iterable
import configparser
from math import isfinite
import os
import tempfile
import warnings
import numpy as np
import pandas as pd
import pyjams as pj  # for mad and all functions about air humidity
import xarray as xr
from musica_qair import q_air_eair, e_air_sat #


__all__ = ['prepSiteForcing']


# Development options
use_musica_qair = False  # Use MuSICA code for qair/rhair/eair conversions
reproduce_old = False    # reproduce old version of prep_site_forcing

#
# Constants
#

# Molar mass of water (:math:`kg mol^{-1}`)
molmass_h2o = 18.01528e-3
# Molar mass of dry air of standard atmosphere (:math:`kg mol^{-1}`)
molmass_air = 28.9644e-3

#
# Functions
#

def str2bool(istr, default):
    """
    Return istr as bool if given, otherwise default

    """
    if istr:
        if istr == 'True':
            return True
        elif istr == 'False':
            return False
        else:
            raise ValueError(f'Boolean string not known: {istr}')
    else:
        return default

def str2float(istr, default):
    """
    Return istr as float if given, otherwise default

    """
    if istr:
        return float(istr)
    else:
        return default

def str2int(istr, default):
    """
    Return istr as int if given, otherwise default

    """
    if istr:
        return int(istr)
    else:
        return default

def parse_entry(text, default=''):
    """
    Convert text string to correct data type if given, otherwise default

    Parse an entry field to None, bool, int, float, datetime, list, dict

    Parameters
    ----------
    text : str
        String from entry field
    default :
        Default value to take if empty string

    Returns
    -------
    None, bool, int, float, datetime, list, dict

    Examples
    --------
    >>> parse_entry('7', -1)
    7
    >>> parse_entry('7,3', [0, 0])
    [7, 3]

    """
    if text is None:
        tt = None
    elif not text:
        tt = default
    elif ',' in text:
        # # list or str
        # try:
        #     tt = eval(f'[{text}]')
        # except SyntaxError:
        #     tt = text
        # parse each element
        stext = text.split(',')
        tt = [ parse_entry(ss, default) for ss in stext ]
    elif text == 'None':
        # None
        tt = None
    elif text == 'True':
        # bool True
        tt = True
    elif text == 'False':
        # bool False
        tt = False
    elif ':' in text:
        # dict, datetime, or str
        try:
            tt = eval(f'{{{text}}}')
        except SyntaxError:
            try:
                tt = np.datetime64(text)
            except ValueError:
                tt = text
    elif text.count('-') == 2:
        # datetime or str
        try:
            tt = np.datetime64(text)
        except ValueError:
            tt = text
    else:
        tt = text

    # if above gave str, check for scalars
    if tt == text:
        try:
            # int
            tt = int(text)
        except ValueError:
            try:
                # float
                tt = float(text)
            except ValueError:
                # str
                tt = text
            try:
                if not isfinite(tt):
                    # keep NaN and Inf string
                    tt = text
            except TypeError:
                pass
        except TypeError:
            # e.g. None
            pass
    return tt

# Line for curve_fit
def flin(x, a, b):
    """
    Linear equation to fit with curve_fit
    y = a + b * x

    """
    return a + b * x


# Line and cost function for fmin (reproduce_old)
def line_p(x, p):
    """
    Straight line: :math:`a + b*x`

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=2`)
          - `p[0]` = a
          - `p[1]` = b

    Returns
    -------
    float
        function value(s)
    """
    return p[0] + p[1] * x


def cost_line(p, x, y):
    """
    Sum of absolute deviations of obs and straight line: :math:`a + b*x`

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)
          - `p[0]` = a
          - `p[1]` = b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y - line_p(x, p)))

#
# Class
#

class prepSiteForcing(object):
    """
    Class for preparing meteo data to run ecosystem models

    Parameters
    ----------
    verbose : bool, optional
        Report progress if True.
    outtype : str, optional
        Plot fits between data and ERA5. Output type is pdf or png;
        everything else opens a screen window (default: no plotting).

    """

    def __init__(self, verbose=False):

        self.verbose = verbose

        # output_variables
        self.varnames = []
        self.varunits = []
        self.dvarunits = {}

        # plotting
        self.plotdict = None

        #
        # configfile
        #
        # Model
        self.ecomodel = ''
        # Options
        self.input = ''
        self.mad_z = 0
        self.imputation_method = 1
        self.make_netcdf = True
        self.keep_csv = False
        # Site
        self.site_name = ''
        self.latitude = -9999.
        self.longitude = -9999.
        self.altitude = 0.
        self.reference_height = 2.
        # ISBA
        self.reference_height_wind = self.reference_height
        self.slope = 0.
        self.aspect = 0.
        # MuSICA
        self.time2gmt = 0.
        self.rsl_yoyo = False
        # Input
        self.infile = ''
        self.sep = None
        self.header = 'infer'
        self.index_col = None
        self.usecols = None
        self.skiprows = None
        self.na_values = None
        self.parse_dates = None
        self.date_format = 'ISO8601'
        self.ftimestep = 1.0
        # Output
        self.outputfile = ''
        self.fill_value = ''
        self.startdate = ''
        self.enddate = ''
        self.timestep = ''
        # ICOS
        self.icos_product = ''
        self.icos_meteo = ''
        self.icos_qc = 2
        # ERA5
        self.era5path = '.'
        self.era5type = 'era5-land-ts'
        # CO2
        self.co2file = ''
        self.co2delimiter = ','
        self.co2date_column = 0
        self.co2co2_column = 1
        # VarNames
        self.dnames = {}
        # VarUnits
        self.dunits = {}
        # AlternativeVarNames
        self.anames = {}

        # Set
        #     varnames, varunits, pnames, punits, dvarunits
        self.output_variables()

        return

    #
    # Output variables names and units
    #

    def output_variables(self):
        """
        Output variable names and units

        Returns
        -------
        Sets attributes: 
            varnames, varunits,
            pnames, punits,
            dvarunits

        """
        # output variables and units
        self.varnames = ['co2air', 'swdown', 'lwdown', 'psurf', 'qair',
                         'tair', 'wind_speed', 'wind_dir']
        self.varunits = ['ppm', 'W/m2', 'W/m2', 'Pa', 'kg/kg',
                         'K', 'm/s', 'degree']
        # precipitation variables
        self.pnames = ['precip', 'rainf', 'snowf']
        self.punits = ['mm', 'mm', 'mm']

        # dictionary with lowercase variables as keys
        self.dvarunits = dict(zip(self.varnames, self.varunits))

        return

    #
    # Config file
    #

    def read_configfile(self, configfile=''):
        """
        Read config file

        Parameters
        ----------
        configfile : string, optional
            Configuration files (default: self.configfile)

        Returns
        -------
        Sets attributes from configfile:
            ecomodel,
            site_name, latitude, longitude, altitude
            reference_height, reference_height_wind,
            time2gmt, slope, aspect,
            infile, sep, header, index_col, usecols,
            skiprows, na_values, parse_dates, date_format,
            outputfile,
            era5path
            startdate, enddate, timestep,
            ftimestep,
            imputation_method, keep_csv, rsl_yoyo
            fill_value,
            varnames, varunits
            dnames, dunits, anames

        """
        if configfile == '':
            configfile = self.configfile

        if configfile == '':
            raise ValueError('Config file must be given')

        cfg = configparser.ConfigParser(interpolation=None)
        cfg.optionxform = str  # preserve case of keys
        cfg.read(configfile)

        # Check mandatory sections and options
        sect = ['Model', 'Site', 'VarNames']
        for ss in sect:
            if not cfg.has_section(ss):
                raise ValueError(f'Mandatory section "{ss}" missing in config file')
        odict = {'Options': 'input'}
        for ss in odict:
            if not cfg.has_option(ss, odict[ss]):
                raise ValueError(f'Mandatory option "{odict[ss]}" in section'
                                 f' "{ss}" missing in config file')

        # Model
        ss = 'Model'
        if not cfg.has_section(ss):
            raise ValueError(f'Mandatory section "{ss}" missing in config file')
        self.ecomodel = cfg['Model'].get('model', '')

        # Options
        if cfg.has_section('Options'):
            self.input = cfg['Options'].get('input', '')
            self.mad_z = str2float(cfg['Options'].get('mad_z', ''), 0.)
            self.imputation_method = str2int(cfg['Options'].get(
                'imputation_method', ''), 1)
            self.make_netcdf = str2bool(cfg['Options'].get(
                'make_netcdf', ''), True)
            self.keep_csv = str2bool(cfg['Options'].get('keep_csv', ''), False)

        # Site
        ss = 'Site'
        if not cfg.has_section(ss):
            raise ValueError(f'Mandatory section "{ss}" missing in config file')
        self.site_name = cfg['Site'].get('name', '')
        self.latitude = str2float(cfg['Site'].get('latitude', ''), -9999.)
        self.longitude = str2float(cfg['Site'].get('longitude', ''), -9999.)
        self.altitude = str2float(cfg['Site'].get('altitude', ''), 0.)
        self.reference_height = str2float(cfg['Site'].get(
            'reference_height', ''), 2.)

        # ISBA
        if cfg.has_section('ISBA'):
            self.reference_height_wind = str2float(cfg['ISBA'].get(
                'reference_height_wind', ''), self.reference_height)
            self.slope = str2float(cfg['ISBA'].get('slope', ''), 0.)
            self.aspect = str2float(cfg['ISBA'].get('aspect', ''), 0.)

        # MuSICA
        if cfg.has_section('MuSICA'):
            self.time2gmt = str2float(cfg['MuSICA'].get('time2gmt', ''), 0.)
            self.rsl_yoyo = str2bool(cfg['MuSICA'].get('rsl_yoyo', ''), False)
        if self.rsl_yoyo:
            self.varnames.append('h_sbl')
            self.varunits.append('m')
            self.dvarunits = dict(zip(self.varnames, self.varunits))

        # Input
        if cfg.has_section('Input'):
            self.infile = cfg['Input'].get('inputfile', '')
            self.sep = cfg['Input'].get('sep', None)
            self.header = parse_entry(cfg['Input'].get('header', 'infer'))
            self.index_col = parse_entry(cfg['Input'].get('index_col', None), None)
            self.usecols = parse_entry(cfg['Input'].get('usecols', None), None)
            self.skiprows = parse_entry(cfg['Input'].get('skiprows', None), None)
            self.na_values = parse_entry(cfg['Input'].get('na_values', None), None)
            self.parse_dates = str2bool(cfg['Input'].get('parse_dates', ''), True)
            self.date_format = cfg['Input'].get('date_format', 'ISO8601')
            self.ftimestep = str2float(cfg['Input'].get('ftimestep', ''), 1.0)

        # Output
        if cfg.has_section('Output'):
            self.outputfile = cfg['Output'].get('outputfile', '')
            self.fill_value = cfg['Output'].get('fill_value', '')
            self.startdate = cfg['Output'].get('startdate', '')
            self.enddate = cfg['Output'].get('enddate', '')
            self.timestep = cfg['Output'].get('timestep', '')  # if no inputfile
        if (self.outputfile == '') or (self.outputfile is None):
            if self.infile:
                self.outputfile = self.infile[0:self.infile.rfind(".")] + '.nc'
            else:
                self.outputfile = 'prep_site_forcing.nc'
        if self.fill_value == '':
            self.fill_value = None
        else:
            self.fill_value = float(self.fill_value)

        # ICOS
        if cfg.has_section('ICOS'):
            self.icos_product = cfg['ICOS'].get('icos_product', '')
            self.icos_meteo = cfg['ICOS'].get('icos_meteo', '')
            self.icos_qc = str2int(cfg['ICOS'].get('icos_qc', ''), 2)

        # ERA5
        if cfg.has_section('ERA5'):
            self.era5path = cfg['ERA5'].get('era', '.')
            self.era5path = cfg['ERA5'].get('era5path', self.era5path)
            self.era5type = cfg['ERA5'].get('era5type', 'era5-land-ts')

        # CO2
        if cfg.has_section('CO2'):
            self.co2file = cfg['CO2'].get('co2file', '')
            self.co2delimiter = cfg['CO2'].get('co2delimiter', ',')
            self.co2date_column = str2int(cfg['CO2'].get('co2date_column', ''), 0)
            self.co2co2_column = str2int(cfg['CO2'].get('co2co2_column', ''), 1)

        # VarNames
        ss = 'VarNames'
        if not cfg.has_section(ss):
            raise ValueError(f'Mandatory section "{ss}" missing in config file')
        # add appropriate precip name
        for iv, vv in enumerate(self.pnames):
            vname = cfg['VarNames'].get(f'name_{vv}', '')
            if vname:
                self.varnames.append(vv)
                self.varunits.append(self.punits[iv])
                self.dvarunits.update({vv: self.punits[iv]})
        # other variable names
        names = []
        for vv in self.varnames:
            names.append(cfg['VarNames'].get(f'name_{vv}', ''))
        self.dnames = dict(zip(self.varnames, names))

        # VarUnits
        if cfg.has_section('VarUnits'):
            units = []
            for vv in self.varnames:
                units.append(cfg['VarUnits'].get(f'unit_{vv}', ''))
        else:
            units = [''] * len(self.varnames)
        self.dunits = dict(zip(self.varnames, units))

        # AlternativeVarNames
        if cfg.has_section('AlternativeVarNames'):
            anames = []
            for vv in self.varnames:
                anames.append(cfg['AlternativeVarNames'].get(
                    f'aname_{vv}', ''))
        else:
            anames = [''] * len(self.varnames)
        self.anames = dict(zip(self.varnames, anames))

        # Extra variables
        if cfg.has_option('VarNames', 'extra_vars'):
            extra_vars = cfg['VarNames'].get('extra_vars', '')
            if extra_vars:
                extra_vars = [ vv.strip()
                               for vv in extra_vars.split(',') ]
            else:
                extra_vars = list()
        else:
            extra_vars = list()

        if cfg.has_option('VarNames', 'extra_names'):
            extra_names = cfg['VarNames'].get('extra_names', '')
            if extra_names:
                extra_names = [ vv.strip()
                                for vv in extra_names.split(',') ]
            else:
                extra_names = extra_vars[:]
        else:
            extra_names = extra_vars[:]
        for ii, ee in enumerate(extra_names):
            self.dnames.update({ee: extra_vars[ii]})
            self.dunits.update({ee: ''})

        # Check config
        leco = self.ecomodel.lower()
        if self.make_netcdf:
            a2n = f'ascii2{leco}.py'
            if not os.path.exists(a2n):
                raise ValueError(f'Cannot find {a2n} to make netcdf'
                                 f' file for ecosystem model {self.ecomodel}.')

        if self.input.lower() not in ['file', 'icos', 'era5']:
            raise ValueError('Input option no known.')

        if self.input.lower() == 'era5':
            if (self.startdate == '') or (self.enddate == ''):
                raise ValueError(
                    'startdate and enddate must be given if no input file'
                    ' nor ICOS data.')
            if self.imputation_method != 1:
                warnings.warn('\ninput=ERA5 implies imputation_method=1.'
                              ' Setting imputation_method=1.')
                self.imputation_method = 1

        if self.imputation_method > 1:
            raise ValueError(f'imputation_method = {self.imputation_method} not'
                             f' implemented.')

        if (not self.make_netcdf) and (not self.keep_csv):
            warnings.warn('\nkeep_csv must be True if make_netcdf is False.'
                          ' Setting keep_csv=True.')
            self.keep_csv = True

        if ((leco == 'musica') and self.rsl_yoyo and
            (self.era5type.lower() not in
             ['era5', 'era5-ts', 'era5-timeseries'])):
            warnings.warn(f'\nera5type must be era5 or era5-ts if rsl_yoyo'
                          f' for {self.ecomodel}. Taking era5-ts.')
            self.era5type = 'era5-ts'

        if ((leco == 'musica') and (self.fill_value is not None)):
            warnings.warn(f'\nfill_value should be empty for {self.ecomodel}')

        if ((leco == 'isba') and (self.fill_value != -9999999.)):
            warnings.warn(f'\nfill_value should be -9999999 for {self.ecomodel}')

        return

    #
    # Plotting
    #

    def setup_plot(self, outtype=''):
        """
        Setup output plots

        Parameters
        ----------
        outtype : str, optional
            Plot fits between data and ERA5. Output type is pdf or png;
            everything else opens a screen window (default: no plotting).

        Returns
        -------
        Sets the following attributes if outtype:
            nrow, ncol,
            hspace, vspace, right,
            textsize,
            dxabc, dyabc,
            lwidth, elwidth, alwidth,
            msize, mwidth,
            fgcolor, bgcolor,
            llxbbox, llybbox, llrspace, llcspace, llhtextpad,
            llhlength, frameon,
            dpi, transparent, bbox_inches, pad_inches

        """
        if outtype != '':
            import matplotlib as mpl
            import matplotlib.pyplot as plt

            figure = plt.Figure()
            fcb = figure.canvas.get_supported_filetypes()
            supported_file_types = list(fcb.keys())
            plt.close(figure)

            # dimensions
            self.nrow     = 4     # # of rows of subplots per figure
            self.ncol     = 2     # # of columns of subplots per figure
            self.hspace   = 0.30  # x-space between subplots
            self.vspace   = 0.20  # y-space between subplots
            self.gridspec = mpl.gridspec.GridSpec(
                self.nrow, self.ncol,
                wspace=self.hspace, hspace=self.vspace)

            # text
            self.textsize = 12    # standard text size

            # lines/markers
            self.lw  = 1.5  # linewidth
            self.alw = 1.0  # axis line width
            self.ms  = 1.5  # marker size
            self.mew = 1.0  # marker edge width

            # colors
            self.fgcolor  = 'black'
            self.bgcolor  = 'white'

            # png
            self.dpi         = 300
            self.transparent = False
            self.bbox_inches = 'tight'
            self.pad_inches  = 0.035

            # setup plot
            if outtype == 'pdf':
                mpl.use('PDF')
                # from matplotlib.backends.backend_pdf import PdfPages
                mpl.rc('ps', papersize='a4', usedistiller='xpdf')
                mpl.rc('figure', figsize=(8.27, 11.69))  # a4 portrait
            elif outtype in supported_file_types:
                mpl.use('Agg')
                mpl.rc('figure', figsize=(8.27, 11.69))
                mpl.rc('savefig', dpi=self.dpi, format='png')
            else:
                # 80% of a4 portrait
                mpl.rc('figure',
                       figsize=(10. / 12. * 8.27, 10. / 12. * 11.69))

            mpl.rcParams['font.family'] = 'sans-serif'
            mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'

            mpl.rc('axes', linewidth=self.alw,
                   edgecolor=self.fgcolor,
                   facecolor=self.bgcolor, labelcolor=self.fgcolor,
                   prop_cycle=mpl.rcsetup.cycler(
                       'color',
                       ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                        '#bcbd22', '#17becf']))
            mpl.rcParams['boxplot.boxprops.color'] = self.fgcolor
            mpl.rcParams['boxplot.capprops.color'] = self.fgcolor
            mpl.rcParams['boxplot.flierprops.color'] = self.fgcolor
            mpl.rcParams['boxplot.flierprops.markeredgecolor'] = self.fgcolor
            mpl.rcParams['boxplot.whiskerprops.color'] = self.fgcolor
            mpl.rc('figure', edgecolor=self.bgcolor, facecolor=self.bgcolor)
            mpl.rc('font', size=self.textsize)
            mpl.rc('grid', color=self.fgcolor)
            mpl.rc('lines', linewidth=self.lw, color=self.fgcolor)
            mpl.rc('patch', edgecolor=self.fgcolor)
            mpl.rc('path', simplify=False)  # do not remove
            mpl.rc('savefig', edgecolor=self.bgcolor, facecolor=self.bgcolor)
            mpl.rc('text', color=self.fgcolor)
            mpl.rc('xtick', color=self.fgcolor)
            mpl.rc('ytick', color=self.fgcolor)

        return

    
    def plot_filled(self, ds, plotdict, outtype='', plotname=''):
        """
        Setup output plots

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)
        plotdict : dict
            Dictionary with variables for scatter plots
        outtype : str, optional
            Plot fits between data and ERA5. Output type is pdf or png;
            everything else opens a screen window (default: no plotting).
        plotname : str, optional
            Name of output file for plots of types pdf and png.

        Returns
        -------
        Plots on screen window or into file if outtype.

        """
        if outtype != '':
            import matplotlib.pyplot as plt

            figure = plt.Figure()
            fcb = figure.canvas.get_supported_filetypes()
            supported_file_types = list(fcb.keys())
            plt.close(figure)

            self.setup_plot(outtype)

            # set plotting name
            if plotname == '':
                plotname = f'prep_site_forcing_era2obs.{outtype.lower()}'
            print(f'Plot {plotname}')

            # open plot
            if (outtype == 'pdf'):
                from matplotlib.backends.backend_pdf import PdfPages
                pdf_pages = PdfPages(plotname)

            fig = plt.figure()

            vlat, vlon = self.get_era5_latlon_names(ds)
            if vlat == 'latitude':
                emod = 'ERA5'
            else:
                emod = 'ERA5-Land'                
            tit = f'Obs vs. {emod}'                
            fig.text(0.5, 0.95, tit,
                     horizontalalignment='center',
                     fontweight='bold', fontsize='large')

            iplot = -1
            for dd in plotdict:
                if plotdict[dd] is not None:
                    pdict = plotdict[dd]
                    if pdict['fit_era'] is not None:
                        iplot += 1
                        sub = fig.add_subplot(self.gridspec[iplot])
                        sub.text(0.5, 0.9, dd,
                                 horizontalalignment='center',
                                 fontweight='bold', fontsize='large',
                                 transform=sub.transAxes)
                        # data
                        xx = pdict['fit_era']
                        yy = pdict['fit_obs']
                        mark1 = sub.plot(xx, yy, 'o', ms=self.ms)
                        # same x- and y-axis
                        xlim = sub.get_xlim()
                        ylim = sub.get_ylim()
                        lim = [min(xlim[0], ylim[0]), max(xlim[1], ylim[1])]
                        sub.set_xlim(lim)
                        sub.set_ylim(lim)
                        # 1:1 line
                        line0 = sub.plot(lim, lim, ':', lw=0.5*self.lw,
                                         color='grey')
                        # fit
                        popt = pdict['linear_params']
                        xmin, xmax = xx.min(), xx.max()
                        line1 = sub.plot(
                            [xmin, xmax],
                            [flin(xmin, *popt), flin(xmax, *popt)],
                            '-', lw=self.lw)
                        st = f'Obs = {popt[0]:.2f} + {popt[1]:.2f} {emod}'
                        sub.text(0.02, 0.8, st, fontsize='small',
                                 transform=sub.transAxes)

            if outtype == 'pdf':
                pdf_pages.savefig(fig)
                plt.close(fig)
                pdf_pages.close()
            elif outtype in supported_file_types:
                fig.savefig(plotname, transparent=self.transparent,
                            bbox_inches=self.bbox_inches,
                            pad_inches=self.pad_inches)
                plt.close(fig)
            else:
                plt.show()

        return

    #
    # Data
    #

    def read_csv_data(self, infile='', startdate=None, enddate=None):
        """
        Read data, get datetime, select columns, select time_span

        Parameters
        ----------
        infile : string, optional
            csv file with forcing variables (default: self.infile)
        startdate : string, optional
            First possible date in netcdf output file in ISO8601 format.
            (Default: first date in input file)
        enddate : string, optional
            Last possible date in output netcdf file in ISO8601 format.
            (Default: last date in input file)

        Returns
        -------
        df : pandas.DataFrame
           File read into pandas.DataFrame with datetime index

        """
        if infile == '':
            infile = self.infile
        if startdate is None:
            startdate = self.startdate
        if enddate is None:
            enddate = self.enddate

        print(f'Read input file {infile}')

        navalues = ['NaN', 'NA', 'nan', 'NAN']
        if self.na_values is not None:
            if isinstance(self.na_values, Iterable):
                navalues.extend(self.na_values)
            else:
                navalues.append(self.na_values)
        rcargs = {'sep': self.sep,
                  'header': self.header,
                  'index_col': self.index_col,
                  'usecols': self.usecols,
                  'skiprows': self.skiprows,
                  'na_values': self.na_values,
                  'parse_dates': self.parse_dates,
                  'date_format': self.date_format}
        df = pd.read_csv(infile, **rcargs)

        if not np.issubdtype(df.index.dtype, np.datetime64):
            raise ValueError(
                f'DataFrame after reading file must have datetime index.'
                f' Currently it is: {df.index.dtype}.\n'
                f'You might want to set "index_col" to the column index of'
                f' the dates in the input file. Current options to'
                f' pandas.read_csv are:\n'
                f'{rcargs}'
                f'which gives the DataFrame:\n'
                f'{df}')
        try:
            _ = [ float(cc) for cc in df.columns ]
            warnings.warn(f'\nAll column names are numbers in input file:\n'
                          f'{list(df.columns)}.\n'
                          f'You might have set "header" or "skiprows" to wrong'
                          f' numbers. Currently they are: header {self.header},'
                          f' skiprows {self.skiprows}.')
        except ValueError:
            pass

        in_columns = list(df.columns.copy())
        wanted_columns = list(self.dnames.values())

        # remove units from columns names
        cols = [ hh.split('(')[0].rstrip() if '(' in hh
                 else hh.split('[')[0].rstrip()
                 for hh in df.columns ]
        df.rename(columns=dict(zip(df.columns, cols)), inplace=True)

        # aggregate multiple columns
        print('Aggregate variables')
        dnames = self.dnames.copy()
        ll = []
        for dd in self.dnames:   # standard and extra vars
            if self.dnames[dd]:
                da = df.filter(regex=self.dnames[dd], axis=1)
                dvars = list(da.columns)
                if len(dvars) == 0:
                    raise ValueError(f'No column found for {self.dnames[dd]}.\n'
                                     f'Available columns are:\n'
                                     f'{in_columns}')
                vkeep = dvars[0]
                dnames[dd] = vkeep
                ll.append(vkeep)
                if len(dvars) > 1:
                    print(f'    {dd}: {dvars}')
                    da = da.mean(axis=1)
                    df[vkeep] = da
                    vdel = dvars[1:]
                    df = df.drop(columns=vdel)
        if len(ll) == 0:
            raise ValueError(f'No columns left after aggregation of variables.\n'
                             f'Available columns were:\n'
                             f'{in_columns}\n'
                             f'Variables demanded were:\n'
                             f'{wanted_columns}')
        df = df[ll]
        self.dnames = dnames

        # fill standard vars with alternative vars
        for dd in self.anames:
            aa = self.anames[dd]     # alternative var
            if aa:
                da = df.filter(regex=aa, axis=1).mean(axis=1)
                vv = self.dnames[dd] # standard var
                df[vv] = df[vv].where(df[vv].notna(), other=da)

        # select variables
        ll = []
        for dd in self.dnames:   # standard and extra vars
            if self.dnames[dd]:
                ll.append(self.dnames[dd])
        if len(ll) == 0:
            raise ValueError(f'No columns left after selection of variables.\n'
                             f'Available columns were:\n'
                             f'{in_columns}\n'
                             f'Variables demanded were:\n'
                             f'{wanted_columns}')
        df = df[ll]

        # rename vars to standard names,
        # making non-existent variables
        for dd in self.dnames:
            if self.dnames[dd]:
                df.rename(columns={self.dnames[dd]: dd}, inplace=True)
            else:
                df[dd] = np.nan

        # start and end dates        
        if startdate == '':
            startdate = df.index[0]
        else:
            startdate = pd.to_datetime(startdate, format='ISO8601')
        if enddate == '':
            enddate = df.index[-1]
        else:
            enddate = pd.to_datetime(enddate, format='ISO8601')
        was_start = df.index[0]
        was_end = df.index[-1]
        df = df[(df.index >= startdate) & (df.index <= enddate)]
        if len(df) == 0:
            raise ValueError(f'No timesteps left after selecting between'
                             f' startdate {startdate} and enddate {enddate}.\n'
                             f'Available dates were between {was_start} and'
                             f' {was_end}.')

        self.df = df

        return df


    def convert_units(self, df=''):
        """
        Convert variables to output units (self.varunits)

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)

        Returns
        -------
        df : pandas.DataFrame
           Input DataFrame with units converted to output units

        """
        if not isinstance(df, str):
            idf = df
        else:
            idf = self.df.copy()

        for dd in self.dnames:
            if (self.dunits[dd] == '') and (dd in self.varnames):
                self.dunits[dd] = self.dvarunits[dd]

        vpd2qair = False
        rhair2qair = False
        for dd in self.dnames:
            if dd in self.pnames:
                if self.dunits[dd] not in ['mm', 'mm/dt', 'mm dt-1',
                                           'kg/m2/dt', 'kg m-2 dt-1']:
                    raise ValueError(
                        f'Only units mm (kg/m2/dt) possible for precipitation.'
                        f' Given: {self.dunits[dd]}')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'wind_speed':
                if self.dunits[dd] not in ['m/s', 'm s-1', '']:
                    raise ValueError(
                        f'Only units m/s possible for wind speed.'
                        f' Given: {self.dunits[dd]}')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'wind_dir':
                if self.dunits[dd] not in ['degree', 'deg', '°', '']:
                    raise ValueError(f'Only units degree possible for wind'
                                     f' direction. Given: {self.dunits[dd]}')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'co2air':
                if self.dunits[dd] not in ['ppm', 'ppmv',
                                           'µmol/mol', 'µmol mol-1']:
                    raise ValueError(f'Only units ppm, ppmv, and µmol/mol'
                                     f' possible for CO2air.'
                                     f' Given: {self.dunits[dd]}')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'qair':
                if self.dunits[dd] in ['mol/mol', 'mol mol-1']:
                    # mol/mol -> kg/kg
                    idf[dd] *= molmass_h2o / molmass_air
                elif self.dunits[dd] in ['mmol/mol', 'mmol mol-1']:
                    # mmol/mol -> kg/kg
                    idf[dd] *= 1e-3 * molmass_h2o / molmass_air
                elif self.dunits[dd] in ['kPa', 'hPa', 'Pa']:
                    # ICOS var == VPD
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    if self.dunits[dd] == 'kPa':
                        idf[dd] *= 1000.
                    elif self.dunits[dd] == 'hPa':
                        idf[dd] *= 100.
                    vpd2qair = True
                elif self.dunits[dd] in ['%', 'percent', '', '0-1']:
                    # observed var == relative humidity
                    # -> calc water vapour mixing ratio after change of
                    # units of tair and psurf
                    rhair2qair = True
                elif self.dunits[dd] == self.dvarunits[dd]:
                    continue  # kg/kg
                else:
                    raise ValueError(
                        f'Cannot convert from qair unit'
                        f' {self.dunits[dd]} to {self.dvarunits[dd]}.')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'tair':
                if (self.dunits[dd] in ['C', 'degreeC', 'degree C',
                                        'degC', 'deg C', '°C']):
                    # degC -> K
                    idf[dd] += 273.15
                elif self.dunits[dd] == self.dvarunits[dd]:
                    continue  # K
                else:
                    raise ValueError(f'Cannot convert from tair unit'
                                     f' {self.dunits[dd]} to'
                                     f' {self.dvarunits[dd]}.')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'psurf':
                if self.dunits[dd] in ['hPa', 'mbar']:
                    # hPa -> Pa
                    idf[dd] *= 100.
                elif self.dunits[dd] == 'kPa':
                    # kPa -> Pa
                    idf[dd] *= 1000.
                elif self.dunits[dd] == 'bar':
                    # bar -> Pa
                    idf[dd] *= 100000.
                elif self.dunits[dd] == self.dvarunits[dd]:
                    continue  # Pa
                else:
                    raise ValueError(f'Cannot convert from psurf unit'
                                     f' {self.dunits[dd]} to'
                                     f' {self.dvarunits[dd]}.')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd == 'h_sbl':
                if self.dunits[dd] not in ['m', '']:
                    raise ValueError(
                        f'Only units m possible for surface'
                        f' boundary layer. Given: {self.dunits[dd]}')
                self.dunits[dd] = self.dvarunits[dd]
            elif dd in self.dvarunits:
                if self.dunits[dd] == self.dvarunits[dd]:
                    continue  # nothing to do
                elif self.dunits[dd] in ['W/m2', 'W m-2']:
                    if self.dvarunits[dd] not in ['W/m2', 'W m-2']:
                        raise ValueError(
                            f'Cannot convert from unit {self.dunits[dd]}'
                            f' to {self.dvarunits[dd]} for variable {dd}.')
                else:
                    raise ValueError(
                        f'Cannot convert from unit {self.dunits[dd]}'
                        f' to {self.dvarunits[dd]} for variable {dd}.')

        if vpd2qair:  # vpd2shair
            vpd = idf['qair']  # Pa
            tk = idf['tair']   # K
            p = idf['psurf']   # Pa
            if use_musica_qair:
                eair = e_air_sat(tk) - vpd
                idf['qair'] = q_air_eair(eair, p)
            else:
                eair = pj.vpd2eair(vpd, tk)
                idf['qair'] = pj.eair2shair(eair, p)

        if rhair2qair:  # rhair2shair
            rh = idf['qair']     # 0-1
            if np.any(rh > 2.):  # %
                rh *= 0.01
            tk = idf['tair']
            p = idf['psurf']
            if use_musica_qair:
                # % -> kg/kg
                eair = rh * e_air_sat(tk)
                idf['qair'] = q_air_eair(eair, p)
            else:
                eair = rh * pj.esat(tk)
                idf['qair'] = pj.eair2shair(eair, p)

        if not isinstance(df, str):
            return idf
        else:
            self.df = idf
            return idf


    def make_empty_data(self, startdate=None, enddate=None, timestep=''):
        """
        Make DataFrame with output variables that are all NaN

        Parameters
        ----------
        startdate : string, optional
            First date in netcdf output file in ISO8601 format.
            (Default: startdate from config file)
        enddate : string, optional
            Last date in output netcdf file in ISO8601 format.
            (Default: enddate from config file)
        timestep : str, optional
            Timestep of output dates. 
            (Default: timestep from config file otherwise 3600s)

            Notation from:
            https://pandas.pydata.org/docs/user_guide/timeseries.html#timeseries-offset-aliases

        Returns
        -------
        df : pandas.DataFrame

        """
        if startdate is None:
            startdate = self.startdate
        if enddate is None:
            enddate = self.enddate
        if timestep == '':
            if self.timestep == '':
                timestep = '3600s'
            else:
                timestep = self.timestep

        print('Make DataFrame with NaNs')

        dates = pd.date_range(start=startdate, end=enddate,
                              freq=f'{timestep}',
                              name='datetime', inclusive='both')
        columns = list(self.dvarunits.keys())

        if 'precip' in columns:
            columns.remove('precip')
            columns.extend(['rainf', 'snowf'])
            self.dvarunits.update({'rainf': 'mm', 'snowf': 'mm'})

        self.dnames = dict(zip(columns, columns))
        self.dunits = self.dvarunits.copy()
        
        df = pd.DataFrame(index=dates, columns=columns)

        self.df = df

        return df


    def read_icos_data(self, station='', product='', meteo='',
                       startdate=None, enddate=None):
        """
        Read data, get datetime, select columns, select time_span

        Parameters
        ----------
        station : str, optional
            ICOS Station name (default: self.site_name)
        product : str, optional
            ICOS-CP data product (default: self.icos_product)
              'NRT' : near-real-time data
        
              'L2' : ICOS L2 data
        
              'Fluxnet' : europe-fluxdata.eu data
        meteo : str, optional
            NRT : 'Meteo', 'Meteosens'

            L2 : 'Meteo', 'Meteosens', 'Fluxnet'
        
            Fluxnet : ignored

            (default: self.icos_meteo)
        startdate : string, optional
            First possible date in netcdf output file in ISO8601 format.
            (Default: first date in input file)
        enddate : string, optional
            Last possible date in output netcdf file in ISO8601 format.
            (Default: last date in input file)

        Returns
        -------
        df : pandas.DataFrame
           File read into pandas.DataFrame with datetime index

        """
        from icos import read_icos

        if station == '':
            station = self.site_name
        if product == '':
            product = self.icos_product
        if meteo == '':
            meteo = self.icos_meteo
        if startdate is None:
            startdate = self.startdate
        if enddate is None:
            enddate = self.enddate

        print(f'Read ICOS product {product} with meteo {meteo}')

        df, dfunit = read_icos(station, product=product, meteo=meteo,
                               units=True, concat=True)
        in_columns = list(df.columns.copy())
        wanted_columns = list(self.dnames.values())
        icos_units = dict(zip(list(df.columns), dfunit))

        # filter quality flags
        if self.icos_qc < 2:
            print('Filter quality flags')
            for cc in df.columns:
                cc_qc = f'{cc}_QC'
                if cc_qc in df.columns:
                    mask = df[cc_qc] <= self.icos_qc
                    print(f'    {cc}: {df.shape[0] - mask.sum()}')
                    df[cc] = df[cc].where(mask, other=np.nan)

        # select variables
        print('Aggregate variables')
        dnames = self.dnames.copy()
        ll = []
        for dd in self.dnames:   # standard and extra vars
            if self.dnames[dd]:
                da = df.filter(regex=self.dnames[dd], axis=1)
                dvars = list(da.columns)
                if len(dvars) == 0:
                    raise ValueError(f'No column found for {self.dnames[dd]}.\n'
                                     f'Available columns are:\n'
                                     f'{in_columns}')
                vkeep = dvars[0]
                dnames[dd] = vkeep
                ll.append(vkeep)
                if len(dvars) > 1:
                    print(f'    {dd}: {dvars}')
                    da = da.mean(axis=1)
                    df[vkeep] = da
                    vdel = dvars[1:]
                    df = df.drop(columns=vdel)
        if len(ll) == 0:
            raise ValueError(f'No columns left after aggregation of variables.\n'
                             f'Available columns were:\n'
                             f'{in_columns}\n'
                             f'Variables demanded were:\n'
                             f'{wanted_columns}')
        df = df[ll]
        self.dnames = dnames

        # fill standard vars with alternative vars
        for dd in self.anames:
            aa = self.anames[dd]     # alternative var
            if aa:
                da = df.filter(regex=aa, axis=1).mean(axis=1)
                vv = self.dnames[dd] # standard var
                df[vv] = df[vv].where(df[vv].notna(), other=da)

        # select variables
        ll = []
        for dd in self.dnames:   # standard and extra vars
            if self.dnames[dd]:
                ll.append(self.dnames[dd])
        if len(ll) == 0:
            raise ValueError(f'No columns left after selection of variables.\n'
                             f'Available columns were:\n'
                             f'{in_columns}\n'
                             f'Variables demanded were:\n'
                             f'{wanted_columns}')
        df = df[ll]

        # rename vars to standard names,
        # making non-existent variables
        for dd in self.dnames:
            if self.dnames[dd]:
                df.rename(columns={self.dnames[dd]: dd}, inplace=True)
            else:
                df[dd] = np.nan

        # update unit dictionary
        for dd in self.dnames:   # standard and extra vars
            if self.dnames[dd]:
                self.dunits.update({dd: icos_units[self.dnames[dd]]})

        # start and end dates        
        if startdate == '':
            startdate = df.index[0]
        else:
            startdate = pd.to_datetime(startdate, format='ISO8601')
        if enddate == '':
            enddate = df.index[-1]
        else:
            enddate = pd.to_datetime(enddate, format='ISO8601')
        was_start = df.index[0]
        was_end = df.index[-1]
        df = df[(df.index >= startdate) & (df.index <= enddate)]
        if len(df) == 0:
            raise ValueError(f'No timesteps left after selecting between'
                             f' startdate {startdate} and enddate {enddate}.\n'
                             f'Available dates were between {was_start} and'
                             f' {was_end}.')

        self.df = df

        return df


    def mad_data(self, df='', z=7):
        """
        Apply Mean Absolute Deviation (MAD) filter on data
        setting everything to NaN, which deviates more than
        *z* standard deviations from median.

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)
        z : float, optional
            Input is allowed to deviate maximum *z* (estimators of) standard
            deviations from the median (default: 7)

        Returns
        -------
        df : pandas.DataFrame
           Input DataFrame with filtered data set to NaN

        """
        import matplotlib.pyplot as plt
        if not isinstance(df, str):
            idf = df
        else:
            idf = self.df.copy()

        print(f'Filter columns with MAD z={z}')

        for cc in idf.columns:
            if ((cc not in self.pnames) and
                (cc != 'swdown') and
                idf[cc].notna().any()):
                mask = pj.mad(idf[cc].values, z=z, nozero=True)
                idf.loc[mask, cc] = np.nan
                print(f'    {cc}: {mask.sum()}')

        return idf


    def get_timestep_seconds(self, df=''):
        """
        Get timestep in seconds

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)

        Returns
        -------
        dt : float
           Time step in seconds

        """
        if isinstance(df, str):
            df = self.df

        dt = (df.index[1] - df.index[0]).seconds

        return dt

    #
    # ERA5
    #

    def open_era5_files(self, ifiles):
        """
        Open ERA5 netCDF files or read csv files

        Parameters
        ----------
        ifiles : str or list
            Filename(s)

        Returns
        -------
        pandas.DataFrame or xarray.Dataset

        """
        try:
            ds = xr.open_mfdataset(ifiles)
        except:
            ropt = {'sep': ',',
                    'header': 'infer',
                    'skiprows': None,
                    'index_col': 0,
                    'parse_dates': True}
            dsl = []
            for ff in ifiles:
                dsl.append(pd.read_csv(ff, **ropt))
            # try to combine files with different variables and
            # files with different time steps
            ds = dsl[0]
            for ds1 in dsl[1:]:
                ds, ds1 = ds.align(ds1, join="outer")
                ds1 = ds1.fillna(ds)
                ds = ds.fillna(ds1)

        return ds


    def get_era5_files(self, df=''):
        """
        Get ERA5 filenames

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)

        Returns
        -------
        era5files : list
           List of filenames of ERA5 files

        """
        from get_era5 import get_era5

        if isinstance(df, str):
            df = self.df

        area = f'{self.latitude},{self.longitude}'
        # always take one day before to get rid of
        # cumulative precip of first day
        tmin = df.index.min() - np.timedelta64(1, 'D')
        tmax = df.index.max()
        isoform = '%Y-%m-%d'
        date = tmin.strftime(isoform) + '/' + tmax.strftime(isoform)

        print(f'Calling get_era5.py -a {area} -d {date},'
              f' -p {self.era5path} -r {self.era5type}')
        era5files = get_era5(area=area, date=date,
                             path=self.era5path,
                             reanalysis_model=self.era5type)

        return era5files


    def get_era5_latlon_names(self, ds):
        """
        Get ERA5 names of latitude and longitude

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)

        Returns
        -------
        lat, lon : str
           Names of latitude and longitude dimensions and variables

        """
        if isinstance(ds, (str, list, tuple)):
            ids = self.open_era5_files(ds)
        else:
            ids = ds

        if isinstance(ids, pd.DataFrame):
            if 'latitude' in ids.columns:
                vlat = 'latitude'
            elif 'lat' in ids.columns:
                vlat = 'lat'
            else:
                raise ValueError("Neither 'lat' nor 'latitude' column found")

            if 'longitude' in ids.columns:
                vlon = 'longitude'
            elif 'lon' in ids.columns:
                vlon = 'lon'
            else:
                raise ValueError("Neither 'lon' nor 'longitude' column found")

        elif isinstance(ids, xr.DataArray):
            if 'latitude' in ids.dims:  # era5
                if 'latitude' in ids.variables:
                    vlat = 'latitude'
                else:
                    ids.close()
                    raise ValueError("'latitude' dimension but no variable found")
            elif 'lat' in ids.dims:     # era5-land
                if 'lat' in ids.variables:
                    vlat = 'lat'
                else:
                    ids.close()
                    raise ValueError("'lat' dimension but no variable found")
            else:
                ids.close()
                raise ValueError("Neither 'lat' nor 'latitude' dimension found")

            if 'longitude' in ids.dims:  # era5
                if 'longitude' in ids.variables:
                    vlon = 'longitude'
                else:
                    ids.close()
                    raise ValueError("'longitude' dimension but no variable found")
            elif 'lon' in ids.dims:      # era5-land
                if 'lon' in ids.variables:
                    vlon = 'lon'
                else:
                    ids.close()
                    raise ValueError("'lon' dimension but no variable found")
            else:
                ids.close()
                raise ValueError("Neither 'lon' nor 'longitude' dimension found")

        else:
            try:
                ids.close()
            except:
                pass
            raise ValueError(f"Could not get lat/lon from ERA5 files: {ds}")

        if isinstance(ds, str) and isinstance(ids, xr.DataArray):
            ids.close()

        return vlat, vlon


    def get_era5_time_name(self, ds):
        """
        Get ERA5 name of time variable

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)

        Returns
        -------
        time : str
           Name of time variable

        """
        if isinstance(ds, (str, list, tuple)):
            ids = self.open_era5_files(ds)
        else:
            ids = ds

        if isinstance(ids, (pd.DataFrame, pd.Series)):
            if ids.index.name == 'time':
                vtime = 'time'
            elif ids.index.name == 'valid_time':
                vtime = 'valid_time'
            else:
                raise ValueError("Neither 'time' nor 'valid_time' found")

        elif isinstance(ids, xr.DataArray):
            if 'time' in ids.dims:          # era5
                if 'time' in ids.variables:
                    vtime = 'time'
                else:
                    ids.close()
                    raise ValueError("'time' dimension but no variable found")
            elif 'valid_time' in ids.dims:  # era5-land
                if 'valid_time' in ids.variables:
                    vtime = 'valid_time'
                else:
                    ids.close()
                    raise ValueError("'valid_time' dimension but no variable"
                                     " found")
            else:
                ids.close()
                raise ValueError("Neither 'time' nor 'valid_time' dimension found")

        else:
            try:
                ids.close()
            except:
                pass
            raise ValueError(f"Could not get time from ERA5 files: {ds}")

        if isinstance(ds, str) and isinstance(ids, xr.DataArray):
            ids.close()

        return vtime


    def get_era5_windspeed_names(self, ds):
        """
        Get ERA5 names of windspeed variables

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)

        Returns
        -------
        vu10, vv10 : str
           Name of windspeed variables

        """
        if isinstance(ds, (str, list, tuple)):
            ids = self.open_era5_files(ds)
        else:
            ids = ds

        if isinstance(ids, pd.DataFrame):
            if 'u10' in ids.columns:
                vu10 = 'u10'
            elif 'u10m' in ids.columns:
                vu10 = 'u10m'
            else:
                raise ValueError("Neither 'u10' nor 'u10m' column found")

            if 'v10' in ids.columns:
                vv10 = 'v10'
            elif 'v10m' in ids.columns:
                vv10 = 'v10m'
            else:
                raise ValueError("Neither 'v10' nor 'v10m' column found")

        elif isinstance(ids, xr.DataArray):
            if 'u10' in ids.variables:     # era5
                vu10 = 'u10'
            elif 'u10m' in ids.variables:  # era5-land
                vu10 = 'u10m'
            else:
                ids.close()
                raise ValueError("Neither 'u10' nor 'u10m' variable found")

            if 'v10' in ids.variables:     # era5
                vv10 = 'v10'
            elif 'v10m' in ids.variables:  # era5-land
                vv10 = 'v10m'
            else:
                ids.close()
                raise ValueError("Neither 'v10' nor 'v10m' variable found")

        else:
            try:
                ids.close()
            except:
                pass
            raise ValueError(f"Could not get wind from ERA5 files: {ds}")

        if isinstance(ds, str) and isinstance(ids, xr.DataArray):
            ids.close()

        return vu10, vv10


    def get_era5_timestep_seconds(self, ds):
        """
        Get ERA5 timestep in seconds

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)

        Returns
        -------
        dt : float
           Time step in seconds

        """
        if isinstance(ds, (str, list, tuple)):
            ids = self.open_era5_files(ds)
        else:
            ids = ds

        if isinstance(ids, (pd.DataFrame, pd.Series)):
            dt = (ids.index[1] - ids.index[0]) / np.timedelta64(1, 's')
        elif isinstance(ids, xr.DataArray):
            vtime = self.get_era5_time_name(ids)
            dt = ( (ids[vtime][1] - ids[vtime][0])
                   / np.timedelta64(1, 's') ).values

        if isinstance(ds, str) and isinstance(ids, xr.DataArray):
            ids.close()

        return dt


    def check_forecast(self, ds):
        """
        Check if surface solar radiation downwards (ssrd) is cumulated variable

        Check that surface solar radiation downwards (ssrd) is always below
        1412 * 3600 ~ 5e6 J/m^2, otherwise some variables are cumulated
        values after start of forecast (00:00 UTC).

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)

        Returns
        -------
        bool
           True if forecast variable detected

        """
        if isinstance(ds, (str, list, tuple)):
            ids = self.open_era5_files(ds)
        else:
            ids = ds

        if isinstance(ids, pd.DataFrame):
            ssrd = ids['ssrd']
        elif isinstance(ids, xr.DataArray):
            ssrd = ids['ssrd'].compute()
        smax = ssrd.max()

        if isinstance(ds, str) and isinstance(ids, xr.DataArray):
            ids.close()
        
        if smax > 5e6:
            return True
        else:
            return False


    def check_era5_variables(self, ds):
        """
        Check that all variables needed are in ERA5 files

        Parameters
        ----------
        ds : string or xarray.Dataset, optional
            If string then xarray.open_mfdataset(ds)

        Returns
        -------
        raise ValueError if variable is missing

        """
        if isinstance(ds, (str, list, tuple)):
            ids = self.open_era5_files(ds)
        else:
            ids = ds

        # ERA5 variable names needed for gap filling
        evars = ['ssrd', 'strd', 'sp', 'd2m', 'sp', 'tp', 't2m']
        evars.extend(self.get_era5_windspeed_names(ids))
        if self.rsl_yoyo:
            evars.append('blh')

        nn = []
        if isinstance(ids, pd.DataFrame):
            for vv in evars:
                if vv not in ids.columns:
                    nn.append(vv)
        elif isinstance(ids, xr.DataArray):
            for vv in evars:
                if vv not in ids.variables:
                    nn.append(vv)

        if len(nn) > 0:
            try:
                ids.close()
            except:
                pass
            raise ValueError(f'{nn} not in ERA5 files')

        if isinstance(ds, str) and isinstance(ids, xr.DataArray):
            ids.close()

        return


    def era5_cumulative_to_hourly(self, dsvar):
        """
        Convert ERA5-Land forecast variable from cumulated values after
        start of forecast (00:00 UTC) to single hour cumulated values.

        Parameters
        ----------
        dsvar : xarray.DataArray
            DataArray of ERA5 forecast variable

        Returns
        -------
        numpy.ndarray
            var with hourly values

        """
        # cumulated of hour is diff to timestep before
        val = dsvar.values
        out = np.diff(val, axis=0, prepend=0)
        # except in the first hour after forecast (01:00 UTC)
        ii = np.where(dsvar['time'].dt.hour == 1)[0]
        out[ii, ...] = val[ii, ...]

        return out

    #
    # Fill data
    #

    def linear_fit_data_vs_era5(self, dfvar, dsvar, daily=False):
        """
        Convert ERA5-Land forecast variable from cumulated values after
        start of forecast (00:00 UTC) to single hour cumulated values.

        Parameters
        ----------
        dfvar : pandas.DataSeries
            DataSeries of observed data
        dsvar : xarray.DataArray
            DataArray of ERA5 forecast variable
        daily : bool, optional
            If True, aggregate to daily values before linear fit,
            forward fill the daily values to whole day

        Returns
        -------
        numpy.ndarray
            dfvar but with fitted values

        """
        import scipy.optimize as opt

        pdict = {'linear_params': None,
                 'fit_era': None, 'fit_obs': None,
                 'fitted_era': None, 'fitted_obs': None,
                 'interpol_era': None, 'interpol_obs': None}
        if daily:
            df = dfvar.resample('1D').mean()
            if isinstance(dsvar, (pd.DataFrame, pd.Series)):
                ds = dsvar.resample('1D').mean()
                ds = ds.loc[(ds.index >= df.index[0]) &
                            (ds.index <= df.index[-1])]
            elif isinstance(dsvar, xr.DataArray):
                ds = dsvar.resample(time='1D').mean()
                ds = ds.loc[(ds['time'] >= df.index[0]) &
                            (ds['time'] <= df.index[-1])]
            
            if any(df.notna()):
                ivar = df[df.notna()].values
                if isinstance(dsvar, (pd.DataFrame, pd.Series)):
                    evar = ds[df.notna()].values
                elif isinstance(dsvar, xr.DataArray):
                    evar = ds[df.notna(), :].values.squeeze()
                if reproduce_old:
                    popt = opt.fmin(cost_line, np.array([0., 1.]),
                                    args=(evar, ivar), disp=False)
                else:
                    popt, pcov = opt.curve_fit(flin, evar, ivar, p0=[0, 1],
                                               nan_policy='omit', maxfev=5000)
                pdict.update({'linear_params': popt,
                              'fit_era': evar, 'fit_obs': ivar})
                out = flin(ds, *popt)
                print(f'        daily linear bias: {popt}')
            else:
                out = ds
                print('        all daily ERA5')
            pdict.update({'fitted_era': ds, 'fitted_obs': out})

            # interpolate to original dfvar using the same value for
            # one hour
            dt = self.get_timestep_seconds(dfvar)
            x = dfvar.index.values.astype('int64') / 1e9
            if isinstance(dsvar, (pd.DataFrame, pd.Series)):
                iout = out.resample(f'{dt}s').ffill()
                iout.index = iout.index + pd.Timedelta(f'{dt/2}s')
                # last day missing thus use interp to fill
                xp = iout.index.values.astype('int64') / 1e9
                fp = iout.values
            elif isinstance(dsvar, xr.DataArray):
                iout = out.resample(time=f'{dt}s').ffill()
                iout['time'] = iout['time'] + pd.Timedelta(f'{dt/2}s')
                # last day missing thus use interp to fill
                xp = iout['time'].values.astype('int64') / 1e9
                fp = iout.values.squeeze()
            iout = np.interp(x, xp, fp)
        else:
            # same time step as ERA5
            # will also put index at beginning of the hour
            df = dfvar.resample('1h').mean()
            if isinstance(dsvar, (pd.DataFrame, pd.Series)):
                ds = dsvar.loc[(dsvar.index >= df.index[0]) &
                               (dsvar.index <= df.index[-1])]
            elif isinstance(dsvar, xr.DataArray):
                ds = dsvar.loc[(dsvar['time'] >= df.index[0]) &
                               (dsvar['time'] <= df.index[-1])]

            if any(df.notna()):
                ivar = df[df.notna()].values
                if isinstance(dsvar, (pd.DataFrame, pd.Series)):
                    evar = ds[df.notna()].values
                elif isinstance(dsvar, xr.DataArray):
                    evar = ds[df.notna(), :].values.squeeze()
                if reproduce_old:
                    popt = opt.fmin(cost_line, np.array([0., 1.]),
                                    args=(evar, ivar), disp=False)
                else:
                    popt, pcov = opt.curve_fit(flin, evar, ivar, p0=[0, 1])
                out = flin(dsvar, *popt)
                print(f'        linear bias: {popt}')
                pdict.update({'linear_params': popt,
                              'fit_era': evar, 'fit_obs': ivar})
            else:
                out = dsvar
                print('        all ERA5')
            pdict.update({'fitted_era': dsvar, 'fitted_obs': out})

            # interpolate to original dfvar
            # dfvar.index is at 15 / 45 min.
            # dsvar['time'] is at 0 min but 30 min should be correct
            x = dfvar.index.values.astype('int64') / 1e9
            if isinstance(dsvar, (pd.DataFrame, pd.Series)):
                if reproduce_old:
                    xp = dsvar.index.values
                else:
                    xp = dsvar.index.values + np.timedelta64(30, 'm')
            elif isinstance(dsvar, xr.DataArray):
                if reproduce_old:
                    xp = dsvar['time'].values
                else:
                    xp = dsvar['time'].values + np.timedelta64(30, 'm')
            xp = xp.astype('int64') / 1e9
            fp = out.values.squeeze()
            iout = np.interp(x, xp, fp)
        pdict.update({'interpol_era': dfvar, 'interpol_obs': iout})

        return iout, pdict


    def impute_data(self, df, dsvar=None, minimum=None,
                    imputation_method=None):
        """
        Bias correct data

        Parameters
        ----------
        df : pandas.DataSeries
            DataSeries of observed data
        dsvar : xarray.DataArray, optional
            DataArray of ERA5 forecast variable; only needed
            if `imputation_method == True`
        minimum : float, optional
            If not None, limit output to *minimum*
        imputation_method : int, optional
            Method for imputation of missing data (gap filling):
            0: linear interpolation between data,
            1: linearly corrected ERA5 data
            (default: self.imputation_method)

        Returns
        -------
        pandas.DataSeries, dict
            df with NaN values filled, dict with plotting variables

        """
        pdict = {'linear_params': None,
                 'fit_era': None, 'fit_obs': None,
                 'fitted_era': None, 'fitted_obs': None,
                 'interpol_era': None, 'interpol_obs': None}
        if imputation_method is None:
            imputation_method = self.imputation_method

        if imputation_method == 0:
            ii = df.notna()
            if ii.sum() == 0:
                warnings.warn(f'\nNo valid data point in {df.name}')
                return df, pdict
            ivar = np.interp(df.index,
                             df.index[ii],
                             df.loc[ii])
            pdict.update({'interpol_obs': ivar})
        elif imputation_method == 1:
            ivar, pdict = self.linear_fit_data_vs_era5(df, dsvar)

        out = df.where(df.notna(), other=ivar)

        if minimum is not None:
            out = out.where(out > minimum, other=minimum)

        return out, pdict


    def fill_data(self, df='', imputation_method=None, outtype='', plotname=''):
        """
        Fill missing data

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame of observed data
        imputation_method : int, optional
            Method for imputation of missing data (gap filling):
            0: linear interpolation between data,
            1: linearly corrected ERA5 data
            (default: self.imputation_method)
        outtype : str, optional
            Plot fits between data and ERA5. Output type is pdf or png;
            everything else opens a screen window (default: no plotting).
        plotname : str, optional
            Name of output file for plots of types pdf and png.

        Returns
        -------
        pandas.DataFrame
            df with NaN values filled with ERA5 values
        Plots on screen window or into file if outtype.

        """
        if isinstance(df, str):
            df = self.df
        dt = self.get_timestep_seconds(df)

        if imputation_method is None:
            imputation_method = self.imputation_method

        if self.imputation_method == 1:
            era5files = self.get_era5_files(df)
            print(f'Open ERA5 files {era5files}')
            ds = self.open_era5_files(era5files)
            era5dt = self.get_era5_timestep_seconds(ds)
            iforecast = self.check_forecast(ds)

        print('Fill')
        plotdict = dict()
        for dd in self.dnames:
            if (dd == 'precip') or (dd == 'rainf'):
                pkey = 'precip'
            else:
                pkey = dd
            plotdict.update({pkey: None})
            if any(df[dd].isna()):
                print(f'    {dd}')
                if dd == 'swdown':
                    if self.imputation_method == 1:
                        evar = ds['ssrd']
                        if iforecast:
                            evar.values = self.era5_cumulative_to_hourly(evar)
                        evar = evar / era5dt  # mean W/m2
                    else:
                        evar = None
                    plotdict.update({pkey: {'data': df[dd], 'era': evar}})
                    df[dd], pdict = self.impute_data(
                        df[dd], evar, minimum=0.,
                        imputation_method=imputation_method)
                    plotdict[pkey].update(pdict)
                elif dd == 'lwdown':
                    if self.imputation_method == 1:
                        evar = ds['strd']
                        if iforecast:
                            evar.values = self.era5_cumulative_to_hourly(evar)
                        evar = evar / era5dt  # mean W/m2
                    else:
                        evar = None
                    plotdict.update({pkey: {'data': df[dd], 'era': evar}})
                    df[dd], pdict = self.impute_data(
                        df[dd], evar, minimum=0.,
                        imputation_method=imputation_method)
                    plotdict[pkey].update(pdict)
                elif dd == 'psurf':
                    if self.imputation_method == 1:
                        evar = ds['sp']
                    else:
                        evar = None
                    plotdict.update({pkey: {'data': df[dd], 'era': evar}})
                    df[dd], pdict = self.impute_data(
                        df[dd], evar, minimum=0.,
                        imputation_method=imputation_method)
                    plotdict[pkey].update(pdict)
                elif dd == 'qair':
                    if self.imputation_method == 1:
                        d2m  = ds['d2m']
                        pres = ds['sp']
                        evar = pj.eair2shair(pj.esat(d2m), pres)  # kg/kg
                    else:
                        evar = None
                    plotdict.update({pkey: {'data': df[dd], 'era': evar}})
                    df[dd], pdict = self.impute_data(
                        df[dd], evar, imputation_method=imputation_method)
                    plotdict[pkey].update(pdict)
                elif dd == 'tair':
                    if self.imputation_method == 1:
                        evar = ds['t2m']
                    else:
                        evar = None
                    plotdict.update({pkey: {'data': df[dd], 'era': evar}})
                    df[dd], pdict = self.impute_data(
                        df[dd], evar, minimum=0.,
                        imputation_method=imputation_method)
                    plotdict[pkey].update(pdict)
                elif dd == 'wind_speed':
                    if self.imputation_method == 1:
                        vu10, vv10 = self.get_era5_windspeed_names(ds)
                        u10 = ds[vu10]
                        v10 = ds[vv10]
                        evar = np.sqrt(u10**2 + v10**2)
                    else:
                        evar = None
                    plotdict.update({pkey: {'data': df[dd], 'era': evar}})
                    df[dd], pdict = self.impute_data(
                        df[dd], evar, minimum=0.,
                        imputation_method=imputation_method)
                    plotdict[pkey].update(pdict)
                elif dd == 'h_sbl':
                    if self.imputation_method == 1:
                        sbl = ds['blh']
                        evar = sbl * 0.1
                        if isinstance(evar, (pd.DataFrame, pd.Series)):
                            ivar = np.interp(df.index, evar.index, evar)
                        elif isinstance(evar, xr.DataArray):
                            vtime = self.get_era5_time_name(evar)
                            ivar = np.interp(df.index, evar[vtime], evar)
                        df[dd] = df[dd].where(df[dd].notna(), other=ivar)
                elif (dd == 'precip') or (dd == 'rainf'):
                    # data - rain, snow, and total precip
                    if dd == 'precip':
                        rainf = df[dd]
                        snowf = rainf.copy()
                        snowf[:] = 0.
                    else:
                        rainf = df[dd]
                        snowf = df['snowf']
                    dvar = (rainf + snowf) / dt
                    if imputation_method == 0:
                        df[dd] = df[dd].where(df[dd].notna(), other=0.)
                        if dd != 'precip':
                            df['snowf'] = df['snowf'].where(
                                df['snowf'].notna(), other=0.)
                    elif imputation_method == 1:
                        # era5 - total precip
                        evar = ds['tp']  # m
                        # there is always noise in the ERA5 data ~ 1e-18
                        evar = evar.where(evar > np.finfo(float).eps, other=0.)
                        if iforecast:
                            evar.values = self.era5_cumulative_to_hourly(evar)
                        evar = evar * 1000. / era5dt  # mean mm/s = kg/m2/s
                        plotdict.update({pkey: {'data': dvar, 'era': evar}})
                        ivar, pdict = self.linear_fit_data_vs_era5(
                            dvar, evar, daily=True)
                        plotdict[pkey].update(pdict)
                        ivar = np.where(ivar > np.finfo(float).eps, ivar, 0.)
                        ivar *= dt
                        if dd == 'precip':
                            df[dd] = df[dd].where(df[dd].notna(), other=ivar)
                        else:
                            tair = ds['t2m']
                            if isinstance(tair, (pd.DataFrame, pd.Series)):
                                tair = tair.loc[(tair.index >= df.index[0]) &
                                                (tair.index <= df.index[-1])]
                                ttair = tair.index
                            elif isinstance(tair, xr.DataArray):
                                if 'time' in tair:
                                    tair = tair.loc[
                                        (tair['time'] >= df.index[0]) &
                                        (tair['time'] <= df.index[-1])]
                                    ttair = tair['time']
                                elif 'valid_time' in tair:
                                    tair = tair.loc[
                                        (tair['valid_time'] >= df.index[0]) &
                                        (tair['valid_time'] <= df.index[-1])]
                                    ttair = tair['valid_time']
                                else:
                                    ValueError(
                                        'time variable not found in ERA5 data')
                            tair = np.interp(df.index, ttair, tair)
                            rainf = np.where(tair >= 274.15, ivar, 0.)
                            snowf = np.where(tair < 274.15, ivar, 0.)
                            df[dd] = df[dd].where(df[dd].notna(), other=rainf)
                            df['snowf'] = df['snowf'].where(
                                df['snowf'].notna(), other=snowf)
                    else:
                        continue  # do not fill anything
                elif dd == 'snowf':
                    continue  # done with rainf
                else:
                    plotdict.update({pkey: None})
                    continue  # co2air and wind_dir are done separately
            else:
                plotdict.update({pkey: None})

        if self.imputation_method == 1:
            self.plot_filled(ds, plotdict,
                             outtype=outtype, plotname=plotname)

            if isinstance(ds, xr.DataArray):
                ds.close()

        return df


    def fill_co2(self, df='', co2file=''):
        """
        Fill missing CO2air

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)
        co2file : str, optional
            Filename with CO2 concentrations (default: self.co2file)
            If (not co2file) and (not self.co2file), gaps will be
            filled with the mean value of the time series.

        Returns
        -------
        df : pandas.DataFrame
           Input DataFrame with missing co2air values filled

        """
        if not isinstance(df, str):
            idf = df
        else:
            idf = self.df.copy()
        if co2file == '':
            co2file = self.co2file

        print('Fill CO2')

        co2 = idf['co2air']
        if any(co2.isna()):
            if co2file == '':
                ivar = co2.mean()
                print(f'    mean value {ivar}')
            else:
                print(f'    using {co2file}')
                rcargs = {'sep': self.co2delimiter,
                          'header': 0,
                          'parse_dates': False,
                          'skipinitialspace': True}
                dfc = pd.read_csv(co2file, **rcargs)
                aco2 = dfc.iloc[:, self.co2co2_column]
                date = dfc.iloc[:, self.co2date_column].values
                year = (date // 1).astype(int).astype(str)
                doy = ((date % 1.) * 365.25).astype(int).astype(str)
                aco2.index = pd.to_datetime(year + '.' + doy, format='%Y.%j')
                ivar = np.interp(co2.index, aco2.index, aco2)
            idf['co2air'] = co2.where(co2.notna(), other=ivar)

        # all NaN
        if all(idf['co2air'].isna()):
            warnings.warn('\nNo CO2 concentration given. Fill in csv'
                          ' manually and run ascii2netcdf again.')

        if not isinstance(df, str):
            return idf
        else:
            self.df = idf
            return idf


    def fill_wind_dir(self, df=''):
        """
        Fill missing wind_dir with median or zero

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)

        Returns
        -------
        df : pandas.DataFrame
           Input DataFrame with missing wind_dir values filled with median

        """
        if not isinstance(df, str):
            idf = df
        else:
            idf = self.df.copy()

        wdir = idf['wind_dir']
        if any(wdir.isna()):
            idf['wind_dir'] = wdir.where(wdir.notna(), other=wdir.median())

        # all NaN
        if all(idf['wind_dir'].isna()):
            idf['wind_dir'] = 0.

        if not isinstance(df, str):
            return idf
        else:
            self.df = idf
            return idf

    #
    # Check
    #

    def check_vars_filled(self, df=''):
        """
        Check that no NaN left in vars

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)

        Returns
        -------
        raise ValueError if still NaN in any variable

        """
        if isinstance(df, str):
            df = self.df.copy()

        nn = []
        for vv in df.columns:
            if vv != 'co2air':
                if any(df[vv].isna()):
                    nn.append(vv)

        if len(nn) > 0:
            raise ValueError(f'{nn} not completely filled.')

        return

    #
    # Write
    #

    def write_csv(self, df='', csvfile='', keep_csv=None):
        """
        Write csv file

        Parameters
        ----------
        df : pandas.DataFrame, optional
            DataFrame with input data (default: self.df)
        csvfile : str, optional
            Name of intermediate csv file
            (default: outputfile with .csv file ending)
        keep_csv : bool, optional
            Keep intermediate csv file
            (default: self.keep_csv)

        Returns
        -------
        csvfile : str
            Name of csvfile written

        """
        print('Write csv')
        if isinstance(df, str):
            idf = self.df.copy()
        else:
            idf = df.copy()

        if keep_csv is None:
            keep_csv = self.keep_csv

        columns = ['year', 'month', 'day', 'hour', 'minute', 'second']
        columns.extend(idf.columns)

        # make datetime columns
        idf['year'] = idf.index.year
        idf['month'] = idf.index.month
        idf['day'] = idf.index.day
        idf['hour'] = idf.index.hour
        idf['minute'] = idf.index.minute
        idf['second'] = idf.index.second

        # sort datetime columns first
        idf = idf[columns]
        
        # include units in column names
        ocol = {}
        for cc in idf.columns:
            if cc in self.dvarunits:
                ocol.update({cc: f'{cc} ({self.dvarunits[cc]})'})
        if len(ocol) > 0:
            idf.rename(columns=ocol, inplace=True)

        if isinstance(self.na_values, Iterable):
            undef = self.na_values[0]
        elif np.isnan(self.na_values) or (self.na_values is None):
            undef = 'NaN'
        else:
            undef = self.na_values

        # write file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as ff:
            tmpfile = ff.name
            idf.to_csv(ff, na_rep=undef, index=False)

        if keep_csv or csvfile:
            if csvfile:
                self.csvfile = csvfile
            else:
                outbase = self.outputfile[0:self.outputfile.rfind(".")]
                self.csvfile = f'{outbase}.csv'
                if self.csvfile == self.infile:
                    self.csvfile = f'{outbase}.tmp.csv'
            os.rename(tmpfile, self.csvfile)
            print(f'    {self.csvfile}')
        else:
            self.csvfile = tmpfile

        return self.csvfile


    def write_netcdf(self, csvfile='', keep_csv=None):
        """
        Write netCDF file

        Parameters
        ----------
        csvfile : str, optional
            Name of intermediate csv file
            (default: outputfile with .csv file ending)
        keep_csv : bool, optional
            Keep intermediate csv file
            (default: self.keep_csv)

        Returns
        -------
        netcdf file on disk

        """
        if csvfile:
            csvfile = self.csvfile

        if keep_csv is None:
            keep_csv = self.keep_csv

        leco = self.ecomodel.lower()  # musica, isba
        exec(f'from ascii2{leco} import ascii2{leco}', globals())
        print(f'Write netcdf {self.outputfile}')

        kwargs = {
            # all
            'outputfile': self.outputfile,
            'ftimestep': self.ftimestep,
            'enddate': self.enddate,
            'startdate': self.startdate,
            'fill_value': self.fill_value,
            'verbose': self.verbose}
        aopt = (f'-o {self.outputfile}'
                f' -t {self.ftimestep}')
        if self.enddate:
            aopt += f' -e {self.enddate}'
        if self.startdate:
            aopt += f' -s {self.startdate}'
        if self.fill_value:
            aopt += f' -f {self.fill_value}'
        if self.verbose:
            aopt += ' -v'
        if leco == 'musica':
            kwargs.update({
                'site_latitude': self.latitude,
                'site_longitude': self.longitude,
                'altitude': self.altitude,
                'time2gmt': self.time2gmt,
                'forcing_height': self.reference_height})
            aopt += (f' --site_latitude {self.latitude}'
                     f' --site_longitude {self.longitude}'
                     f' --site_altitude {self.altitude}'
                     f' --time2gmt {self.time2gmt}'
                     f' --forcing_height {self.reference_height}')
            print(f'Calling python ascii2{leco}.py {aopt} {csvfile}')
        elif leco == 'isba':
            kwargs.update({
                'latitude': self.latitude,
                'longitude': self.longitude,
                'altitude': self.altitude,
                'reference_height': self.reference_height,
                'reference_height_wind': self.reference_height_wind,
                'slope': self.slope,
                'aspect': self.aspect})
            aopt += (f' --latitude {self.latitude}'
                     f' --longitude {self.longitude}'
                     f' --altitude {self.altitude}'
                     f' --reference_height {self.reference_height}'
                     f' --reference_height_wind {self.reference_height_wind}'
                     f' --slope {self.slope}'
                     f' --aspect {self.aspect}')
            print(f'Calling python ascii2{leco}.py {aopt} {csvfile}')

        exec(f'ascii2{leco}(csvfile, **kwargs)')

        if not keep_csv:
            os.remove(csvfile)

        return
