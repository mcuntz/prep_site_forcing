#!/usr/bin/env python
'''
Read data from ICOS Data Portal

History
-------
   * Written 2025 by Matthias Cuntz
     Routines: read_icos
   * keywords product and meteo case-insensitive, Matthias Cuntz, Apr 2026
   * Add Fluxnet product, Matthias Cuntz, Apr 2026
   * Added keyword concat, Matthias Cuntz, Apr 2026
   * Added keyword units, Matthias Cuntz, Apr 2026
   * Make '-' NaN in _QC columns, Matthias Cuntz, Apr 2026
   * Add write_icos, Matthias Cuntz, Apr 2026
   * Add info_icos, Matthias Cuntz, Apr 2026
   * Print variable names (info_icos), Matthias Cuntz, Apr 2026
   * Return Dobj if station and full product name (info_icos),
     Matthias Cuntz, Apr 2026
   * Remove Dobj from icoscp, Matthias Cuntz, Jun 2026
   * Add command line interface, Matthias Cuntz, Jun 2026
   * Remove meteo keyword: use full ICOS product names,
     Matthias Cuntz, Jun 2026
   * Return units as dictionary, Matthias Cuntz, Jun 2026

'''
import os
import numpy as np
import pandas as pd
from icoscp_core.icos import meta, data, ECO_STATION
# To initialize, run the following code interactively once for every machine
# from icoscp_core.icos import auth
# auth.init_config_file()


__all__ = ['info_icos', 'read_icos', 'write_icos']



def _get_header_vars_units(arr):
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
            sarr = arr[ic].split()
            iv = sarr[0]
            iu = ' '.join(sarr[1:])
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

    dunits = dict(zip(ivars, iunits))

    return ivars, dunits


# [ dd.label for dd in meta.list_datatypes() ]
# dtypes = [ dd for dd in meta.list_datatypes()
#            if (dd.label.startswith('ETC') or dd.label.startswith('Fluxnet')) ]
# known_products = [ dd.label for dd in dtypes if dd.has_data_access ]
known_products = ['ETC L2 Fluxes', 'ETC L2 Fluxnet (half-hourly)',
                  'ETC L2 Meteo', 'ETC L2 Meteosens',
                  'ETC NRT Fluxes',
                  'ETC NRT Meteo', 'ETC NRT Meteosens',
                  'Fluxnet Product']


def info_icos(station='', product=''):
    '''
    Info on available ICOS products of ecosystem stations

    Parameters
    ----------
    station : str, optional
        ICOS ecosystem station id (case-sensitive).
        Print available products for station if given,
        otherwise print available ICOS stations and products.
    product : str, optional
        ICOS-CP data product (case insensitive) (default: '')
            'ETC L2 Fluxes' : L2 ecosystem fluxes
            'ETC L2 Fluxnet (half-hourly)' : Fluxnet product for ICOS L2 data
            'ETC L2 Meteo' : L2 aggregated meteorological variables
            'ETC L2 Meteosens' : L2 individual meteorological sensors
            'ETC NRT Fluxes' : NRT ecosystem fluxes
            'ETC NRT Meteo' : NRT aggregated meteorological variables
            'ETC NRT Meteosens' : NRT individual meteorological sensors
            'Fluxnet Product' : Fluxnet Shuttle product for station

        `product` can be comma-separated list of ICOS-CP data products.

    Returns
    -------
    icoscp_core.metaclient.DataObject
        if station and product, otherwise None

    '''
    eco_stations = meta.list_stations(ECO_STATION)
    eco_stations_ids = [ ss.id for ss in eco_stations ]

    if station == '':
        print('Known ecosystem stations:')
        eco_stations_ids.sort()
        print(eco_stations_ids)
        
        print('')
        print('Known products:')
        print(known_products)
        
        return

    if station not in eco_stations_ids:
        raise ValueError(f'Station {station} not known.'
                         f' Known ecosystem stations:\n'
                         f'{eco_stations_ids}')

    istation = [ ee for ee in eco_stations if ee.id == station ][0]
    i_have_product = []
    for pp in known_products:
        dtype = [ dd for dd in meta.list_datatypes() if dd.label == pp ][0]
        # meta.list_data_objects is empty is dtype does not exist
        i_have_product.append(len(
            meta.list_data_objects(dtype, station=istation)))
    products_avail = [ known_products[i] for i in range(len(known_products))
                       if i_have_product[i] != 0 ]
    products_avail = sorted(set(products_avail))
    lproducts_avail = [ pp.lower() for pp in products_avail ]

    print(f'Station {station}')

    if product == '':
        print('')
        print('Available products:')
        print('   ', ', '.join(products_avail))
    else:
        if ',' in product:
            products = product.split(',')
        else:
            products = [product]

        for pp in products:
            if pp.lower() not in lproducts_avail:
                raise ValueError(f'Product "{pp}" not known for station'
                                 f' {station}.')

        dobj_metas = []
        for pp in products:
            dtype = [ dd for dd in meta.list_datatypes()
                      if dd.label == pp ][0]
            srelease = meta.list_data_objects(dtype, station=istation)[0]

            dobj_meta = meta.get_dobj_meta(srelease.uri)
            sdata = data.get_columns_as_arrays(dobj_meta, length=1)
            idf = pd.DataFrame(sdata)

            pvars = list(idf.columns)
            pvars.sort()
            print('')
            print(f'Available variables for product "{pp}":')
            print(', '.join(pvars))

            dobj_metas.append(dobj_meta)

        if len(dobj_metas) == 1:
            return dobj_metas[0]
        else:
            return dobj_metas

    return


def read_icos(station, product='ETC L2 Meteosens',
              units=False, concat=False):
    '''
    List of pandas DataFrame of ICOS-CP data

    Parameters
    ----------
    station : str
        ICOS ecosystem station code (case-sensitive)
    product : str, optional
        ICOS-CP data product (default: 'ETC L2 Meteosens')
            'ETC L2 Fluxes' : L2 ecosystem fluxes
            'ETC L2 Fluxnet (half-hourly)' : Fluxnet product for ICOS L2 data
            'ETC L2 Meteo' : L2 aggregated meteorological variables
            'ETC L2 Meteosens' : L2 individual meteorological sensors
            'ETC NRT Fluxes' : NRT ecosystem fluxes
            'ETC NRT Meteo' : NRT aggregated meteorological variables
            'ETC NRT Meteosens' : NRT individual meteorological sensors
            'Fluxnet Product' : Fluxnet Shuttle product for station

        `product` can be comma-separated list of ICOS-CP data products.

        `product` can also be the filename of a file with a structure
        as written, for example, by `write_icos` (-9999 as NaN allowed),
        for example:
    
          TIMESTAMP_END,CO2 (µmol mol-1),FC (µmol m-2 s-1),...
          2025-12-31 23:00:00,432.47,-9999,...

    units : bool, optional
        Return also dictionary with units for variables
        (default: False)
    concat : bool, optional
        Concat different data streams into one pandas.DataFrame if True
        (default: False)

    Returns
    -------
    (List of) pandas.DataFrame with ICOS-CP data products[, dict of units]

    '''
    # check if product is a file
    if os.path.exists(product):
        df = pd.read_csv(product, index_col=0, parse_dates=True,
                         date_format='ISO8601', na_values='-9999')

        dvars, dunits = _get_header_vars_units(df.columns)
        df.rename(columns=dict(zip(df.columns, dvars)), inplace=True)

        if units:
            return df, dunits
        else:
            return df

    # ICOS Carbon Portal
    eco_stations = meta.list_stations(ECO_STATION)
    eco_stations_ids = [ ss.id for ss in eco_stations ]
    if station not in eco_stations_ids:
        raise ValueError(f'Station {station} not known.'
                         f' Known ecosystem stations:\n'
                         f'{eco_stations_ids}')

    # station
    istation = [ ss for ss in eco_stations if ss.id == station ][0]
    # smeta = meta.get_station_meta(istation)

    # available products for station
    i_have_product = []
    for pp in known_products:
        dtype = [ dd for dd in meta.list_datatypes() if dd.label == pp ][0]
        # meta.list_data_objects is empty is dtype does not exist
        i_have_product.append(len(
            meta.list_data_objects(dtype, station=istation)))
    products_avail = [ known_products[i] for i in range(len(known_products))
                       if i_have_product[i] != 0 ]
    products_avail = sorted(set(products_avail))
    lproducts_avail = [ pp.lower() for pp in products_avail ]

    if ',' in product:
        products = product.split(',')
    else:
        products = [product]

    # check that product(s) exist at station
    for pp in products:
        if pp.lower() not in lproducts_avail:
            raise ValueError(f'Product "{pp}" not known for station'
                             f' {station}.')

    # get products
    df = []
    unit = []
    for pp in products:
        dtype = [ dd for dd in meta.list_datatypes() if dd.label == pp ]
        dtype = dtype[0]

        srelease = meta.list_data_objects(dtype, station=istation)[0]
        dobj_meta = meta.get_dobj_meta(srelease.uri)
        sdata = data.get_columns_as_arrays(dobj_meta)

        # _QC column is '-' if data column is NaN,
        # e.g. in product Fluxnet -> set to NaN so that dtype float
        # in pandas.DataFrame
        for dd in sdata:
            ids = sdata[dd]
            if isinstance(ids[0], str):
                sdata[dd] = np.where(
                    ids == '-', 'NaN', ids).astype(np.float32)

        idf = pd.DataFrame(sdata)

        idf.set_index('TIMESTAMP', inplace=True)
        if 'TIMESTAMP_END' in idf.columns:
            idf.drop(columns='TIMESTAMP_END', inplace=True)
        idf.index.name = 'Date Time'

        df.append(idf)

        if units:
            # units of all variables
            si = dobj_meta.specificInfo
            vunits = {}
            for vv in si.columns:
                vt = vv.valueType
                if hasattr(vt, 'unit'):
                    vunits.update({vv.label: vt.unit})
            # units of selected variables (should be the same to vunits)
            iunit = {}
            for cc in idf.columns:
                iunit.update({cc: vunits[cc]})
            # units of all products
            unit.append(iunit)

    if concat:
        dfc = df[0]
        for idf in df[1:]:
            # assumes same time steps
            dfc = pd.concat([dfc, idf], axis=1)
        df = dfc

        if units:
            unitc = unit[0]
            for iunit in unit[1:]:
                unitc.update(iunit)
            unit = unitc

    if units:
        return df, unit
    else:
        return df


def write_icos(station, outfile, product='ETC L2 Meteosens',
               undef=-9999., verbose=False):
    '''
    Write ICOS-CP data to file

    Parameters
    ----------
    station : str
        ICOS ecosystem station code (case-sensitive)
    outfile : str
        Name of csv output file. If empty, it will be:
        {station}_{product.replace(" ", "_").replace(",", "-")}.csv
    product : str, optional
        ICOS-CP data product (default: 'ETC L2 Meteosens')
            'ETC L2 Fluxes' : L2 ecosystem fluxes
            'ETC L2 Fluxnet (half-hourly)' : Fluxnet product for ICOS L2 data
            'ETC L2 Meteo' : L2 aggregated meteorological variables
            'ETC L2 Meteosens' : L2 individual meteorological sensors
            'ETC NRT Fluxes' : NRT ecosystem fluxes
            'ETC NRT Meteo' : NRT aggregated meteorological variables
            'ETC NRT Meteosens' : NRT individual meteorological sensors
            'Fluxnet Product' : Fluxnet Shuttle product for station

        `product` can be comma-separated list of ICOS-CP data products.

    undef : float, optional
        NaN will be undef in output file
        (default: -9999.)
    verbose : bool, optional
        Report progress if True

    Returns
    -------
    csv file on disk

    '''
    # read icos data
    if verbose:
        print(f'  Get ICOS product "{product}" for station "{station}".')
    df, dfunit = read_icos(station, product=product, units=True, concat=True)

    # include units in column names
    ocol = {}
    for cc in df.columns:
        if cc in dfunit:
            ocol.update({cc: f'{cc} ({dfunit[cc]})'})
    if len(ocol) > 0:
        df.rename(columns=ocol, inplace=True)

    # write csv file
    if verbose:
        print(f'  Write into csv file "{outfile}".')
    df.to_csv(outfile, na_rep=undef, index_label='TIMESTAMP_END')

    return


if __name__ == '__main__':

    import argparse
    import os
    import sys

    info = False
    output = ''
    product = ''
    station = ''

    desc = "Info and retrieval of ICOS ecosystem data."
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=desc)

    hstr = ('If True, print info about ICOS stations, products,'
            ' and variables. Otherwise download product at station.'
            ' (default: False)')
    parser.add_argument('-i', '--info', action='store_true',
                        default=info, dest='info',
                        help=hstr)

    hstr = ('Output filename. If not given, filename will be'
            ' {station}_{product.replace(" ", "_").replace(",", "-")}'
            '.csv. Compression is inferred from the following extensions:'
            ' .gz, .bz2, .zip, .xz, .zst, .tar, .tar.gz, .tar.xz, or'
            ' .tar.bz2 (otherwise no compression).')
    parser.add_argument('-o', '--output', action='store',
                        default=output, dest='output',
                        metavar='output_filename', help=hstr)

    hstr = ('Full product name (e.g "ETC L2 Meteosens")'
            ' (case insensitive); product name must be in quotes.'
            ' Print available variable names for product if info (-i),'
            ' otherwise download product at station. Can be'
            ' comma-separated list of products. (default: "")')
    parser.add_argument('-p', '--product', action='store',
                        default=product, dest='product',
                        metavar='product', help=hstr)

    hstr = ('ICOS station ID such as FR-Hes (case-sensitive).'
            ' Prints available products for station (-i)'
            ' or downloads product at station.'
            ' If empty, prints available ICOS stations and products.'
            ' (default: "").')
    parser.add_argument('station', nargs='?', default=station,
                        metavar='ICOS_ID', help=hstr)

    args = parser.parse_args()

    info = args.info
    output = args.output
    product = args.product
    station = args.station

    # info, or station or product missing
    if info or (station == '') or (product == ''):
        info_icos(station=station, product=product)
    else:
        eco_stations = meta.list_stations(ECO_STATION)
        eco_stations_ids = [ ss.id for ss in eco_stations ]

        # valid station
        if station not in eco_stations_ids:
            print(f'Station {station} not known.'
                  f' Known ecosystem stations:\n'
                  f'{eco_stations_ids}')
        else:
            istation = [ ee for ee in eco_stations if ee.id == station ][0]
            i_have_product = []
            for pp in known_products:
                dtype = [ dd for dd in meta.list_datatypes() if dd.label == pp ][0]
                # meta.list_data_objects is empty is dtype does not exist
                i_have_product.append(len(
                    meta.list_data_objects(dtype, station=istation)))
            products_avail = [ known_products[i] for i in range(len(known_products))
                               if i_have_product[i] != 0 ]
            products_avail = sorted(set(products_avail))
            lproducts_avail = [ pp.lower() for pp in products_avail ]

            if ',' in product:
                products = product.split(',')
            else:
                products = [product]

            # check that product(s) exist at station
            allgood = True
            for pp in products:
                if pp.lower() not in lproducts_avail:
                    allgood = False
                    print(f'Product "{pp}" not known for station {station}.')

            if allgood:
                if info:
                    info_icos(station=station, product=product)
                else:
                    if output == '':
                        output = (
                            f'{station}_'
                            f'{product.replace(" ", "_").replace(",", "-")}'
                            f'.csv')
                    print(f'Write {output} for product "{product}"'
                          f' at {station}')
                    write_icos(station, output, product=product,
                               undef=-9999., verbose=False)
            else:
                info_icos(station=station, product='')
