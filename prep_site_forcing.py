#!/usr/bin/env python
'''
Preparation of meteo data to run ecosystem models

History
-------
   * Written Dec 2025 by Matthias Cuntz, UMR Silva, INRAE
     from script prep_site_forcing.py
   * grib, netcdf, and csv ERA5 files possible, Matthias Cuntz, Apr 2026
   * Allow no input file to get all ERA5 forcing, Matthias Cuntz, Apr 2026

'''
import argparse
import warnings
import matplotlib.pyplot as plt
from class_prep_site_forcing import prepSiteForcing


# Development options
noera5 = False  # Test without ERA5. Output is not be filled then.


if __name__ == '__main__':

    plotname = ''
    outtype  = ''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            'Prepare input data for the ecosystem models MuSICA or ISBA,'
            ' reading the configuration given on the command line.\n\n'
            'Default values will be taken in most cases if options are empty,'
            ' missing or commented out.\n\n'
            'Lots of information is the same as needed by ascii2musica and'
            ' ascii2isba. It is given directly in the config file, e.g.'
            ' longitude = 7.06.'
            '\n\n'
            'Other information is related to the input file such as names of'
            ' variables in the header line, etc.'))
    default = __file__[0:__file__.rfind(".")] + '_'
    hstr = (f'Name of output file for plots of types pdf or png'
            f' (default: {default}).')
    parser.add_argument('-p', '--plotname', action='store', default=plotname,
                        dest='plotname', metavar='plotname', help=hstr)
    hstr = ('Output type is pdf or png. Everything else opens screen'
            ' window (default: no plotting).')
    parser.add_argument('-t', '--type', action='store', default=outtype,
                        dest='outtype', metavar='outtype', help=hstr)
    hstr = 'Configuration file for station in ConfigParser format.'
    parser.add_argument('configfile', nargs='?', default=None,
                        metavar='ConfigurationFile', help=hstr)

    args = parser.parse_args()
    plotname = args.plotname
    outtype  = args.outtype
    configfile = args.configfile

    del parser, args

    # Check that config file given
    if configfile is None:
        raise IOError('Configuration file must be given.')
    cfgfile = configfile

    # Determine outtype from plotname
    if (plotname != '') and (outtype == ''):
        outtype = plotname[plotname.rfind(".") + 1:]
        figure = plt.Figure()
        fcb = figure.canvas.get_supported_filetypes()
        supported_file_types = list(fcb.keys())
        plt.close(figure)
        if outtype not in supported_file_types:
            warnings.warn(f'\nplot outtype {outtype} not in supported file types:\n'
                          f'{supported_file_types}')
            plotname = ''
            outtype = 'X'

    import time as ptime
    t1 = ptime.time()

    prep = prepSiteForcing()
    prep.read_configfile(configfile)

    if prep.input.lower() == 'era5':
        df = prep.make_empty_data()
    elif prep.input.lower() == 'icos':
        df = prep.read_icos_data(product=prep.icos_product,
                                 meteo=prep.icos_meteo)
        df = prep.convert_units(df)
    elif prep.input.lower() == 'file':
        df = prep.read_csv_data(prep.infile)
        df = prep.convert_units(df)
    else:
        raise ValueError(f'Option input={prep.input} unknown.')

    if prep.mad_z > 0:
        df = prep.mad_data(df, z=prep.mad_z)

    if any(df.isna()):
        imputation_method = prep.imputation_method
        if noera5:
            imputation_method = 0
        df = prep.fill_data(df, imputation_method=imputation_method,
                            outtype=outtype, plotname=plotname)
    df = prep.fill_co2(df)
    df = prep.fill_wind_dir(df)
    prep.check_vars_filled(df)

    csvfile = prep.write_csv(df)
    if prep.make_netcdf:
        prep.write_netcdf(csvfile, keep_csv=prep.keep_csv)

    t2    = ptime.time()
    strin = ('[m]: {:.1f}'.format((t2 - t1) / 60.)
             if (t2 - t1) > 60. else '[s]: {:d}'.format(int(t2 - t1)))
    print(f'Time elapsed: {strin}')
