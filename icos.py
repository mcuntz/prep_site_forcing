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

'''
import os
import numpy as np
import pandas as pd
from icoscp_core.icos import meta, data, ECO_STATION
from icoscp.dobj import Dobj
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

    return ivars, iunits


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
        Print available products and meteo for station if given,
        otherwise print available ICOS stations and products.
    product : str, optional
        Print available meteo for product of station if short product name,
        and available variable names if full product name (case insensitive).

    Returns
    -------
    icoscp.dobj.Dobj
        if station and full product name, otherwise None
    Prints info on ICOS ecosystem stations

    '''
    eco_stations = meta.list_stations(ECO_STATION)
    eco_stations_ids = [ ss.id for ss in eco_stations ]

    if station == '':
        print('Known ecosystem stations:')
        eco_stations_ids.sort()
        print(eco_stations_ids)
        
        print('')
        print('Known full product names:')
        print(known_products)
        
        return

    if station not in eco_stations_ids:
        raise ValueError(f'Station {station} not known.'
                         f' Known ecosystem stations: {eco_stations_ids}')

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
        print('Available full product names:')
        print('   ', ', '.join(products_avail))

        products = [ pp.split()[1] if pp.split()[0] == 'ETC' else pp.split()[0]
                     for pp in products_avail ]
        products = sorted(set(products))
        print('')
        print('Available short product names:')
        print('   ', ', '.join(products))
    else:
        if product.lower() in lproducts_avail:
            pp = products_avail[lproducts_avail.index(product.lower())]
            dtype = [ dd for dd in meta.list_datatypes() if dd.label == pp ][0]
            srelease = meta.list_data_objects(dtype, station=istation)
            dobj = Dobj(srelease[0].uri)

            # pvars = dobj.colNames  # only preview-able variables
            pvars = dobj.variables.name.values
            pvars.sort()
            print('')
            print(f'Available variables for product "{pp}":')
            print('   ', ', '.join(pvars))

            return dobj
        else:
            ipa = []
            for pp in products_avail:
                p = pp.split()[1]
                if p == 'Product':
                    p = 'Fluxnet'
                if p.lower() == product.lower():
                    ipa.append(pp)
            if len(ipa) == 0:
                raise ValueError(f'Product {product} not known for station'
                                 f' {station}.')
            products_avail = ipa

    pm = []
    for pp in products_avail:
        if pp.startswith('ETC'):
            ipp = pp[4:]
        elif pp.startswith('Fluxnet'):
            ipp = pp[:7]
        else:
            raise ValueError(f'Unknown product: {pp}')
        if ipp.endswith('(half-hourly)'):
            ipp = ipp[:-14]
        if not ipp.endswith('Fluxes'):
            pm.append(ipp)
    pm.sort()
    print('')
    print('Available product meteo combinations:')
    print('   ', ', '.join(pm))

    return


def read_icos(station, product='NRT', meteo='Meteosens',
              units=False, concat=False):
    '''
    List of pandas DataFrame of ICOS-CP data

    Parameters
    ----------
    station : str
        ICOS ecosystem station code (case-sensitive)
    product : str, optional
        ICOS-CP data product (default: 'NRT')
           'NRT' : near-real-time data
           'L2' : ICOS L2 data
           'Fluxnet' : europe-fluxdata.eu data

        Can also be file with structure (-9999 as NaN allowed),
        as written, for example, by `write_icos`:
    
          TIMESTAMP_END,CO2 (µmol mol-1),FC (µmol m-2 s-1),...
          2025-12-31 23:00:00,432.47,-9999,...
    
    meteo : str, optional
        NRT : 'Meteo', 'Meteosens'
        L2 : 'Meteo', 'Meteosens', 'Fluxnet'
        Fluxnet : ignored
    units : bool, optional
        Return also list of units for columns
        (default: False)
    concat : bool, optional
        Concat different data streams into one pandas.DataFrame if True
        (default: False)

    Returns
    -------
    (List of) pandas.DataFrame with ICOS-CP data products[, list of units]

    '''
    eco_stations = meta.list_stations(ECO_STATION)
    eco_stations_ids = [ ss.id for ss in eco_stations ]
    if station not in eco_stations_ids:
        raise ValueError(f'Station {station} not known.'
                         f' Known ecosystem stations: {eco_stations_ids}')

    isstation = [ ss for ss in eco_stations if ss.id == station ]
    isstation = isstation[0]
    # smeta = meta.get_station_meta(isstation)

    if product.lower() == 'l2':
        if (not meteo) or (meteo.lower() == 'meteosens'):
            products = ['ETC L2 Fluxes', 'ETC L2 Meteosens']
        elif meteo.lower() == 'meteo':
            products = ['ETC L2 Fluxes', 'ETC L2 Meteo']
        elif meteo.lower() == 'fluxnet':
            products = ['ETC L2 Fluxnet (half-hourly)']
        else:
            raise ValueError(f'Meteo not known {meteo}')
    elif product.lower() == 'nrt':
        if (not meteo) or (meteo.lower() == 'meteosens'):
            products = ['ETC NRT Fluxes', 'ETC NRT Meteosens']
        elif meteo.lower() == 'meteo':
            products = ['ETC NRT Fluxes', 'ETC NRT Meteo']
        elif meteo.lower() == 'fluxnet':
            products = ['ETC NRT Fluxnet (half-hourly)']
        else:
            raise ValueError(f'Meteo not known: {meteo}')
    elif product.lower() == 'fluxnet':
        products = ['Fluxnet Product']
    else:
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
        else:
            raise ValueError(f'Product not known: {product}')

    df = []
    unit = []
    for pp in products:
        if pp not in known_products:
            raise ValueError(f'Product/meteo combination not known: {pp}.')

        dtype = [ dd for dd in meta.list_datatypes() if dd.label == pp ]
        dtype = dtype[0]

        srelease = meta.list_data_objects(dtype, station=isstation)
        sdata = data.batch_get_columns_as_arrays(srelease)

        idf = [ arrs for dobj, arrs in sdata ]
        idf = idf[0]
        # _QC column is '-' if data column is NaN,
        # e.g. in product Fluxnet -> set to NaN so that dtype float
        # in pandas.DataFrame
        for dd in idf:
            ids = idf[dd]
            if isinstance(ids[0], str):
                idf[dd] = np.where(ids == '-', 'NaN', ids).astype(np.float32)
        idf = pd.DataFrame(idf)

        idf.set_index('TIMESTAMP', inplace=True)
        if 'TIMESTAMP_END' in idf.columns:
            idf.drop(columns='TIMESTAMP_END', inplace=True)
        idf.index.name = 'Date Time'

        df.append(idf)

        if units:
            dobj = Dobj(srelease[0].uri)
            iunit = []
            for cc in idf.columns:
                iunit.append(dobj.variables[
                    dobj.variables.name == cc].unit.values[0])
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
                unitc.extend(iunit)
            unit = unitc

    if units:
        return df, unit
    else:
        return df


def write_icos(station, outfile,
               product='NRT', meteo='Meteosens',
               undef=-9999., verbose=False):
    '''
    Write ICOS-CP data to file

    Parameters
    ----------
    station : str
        ICOS ecosystem station code (case-sensitive)
    outfile : str
        Name of csv output file
    product : str, optional
        ICOS-CP data product (default: 'NRT')
           'NRT' : near-real-time data
           'L2' : ICOS L2 data
           'Fluxnet' : europe-fluxdata.eu data
    meteo : str, optional
        NRT : 'Meteo', 'Meteosens'
        L2 : 'Meteo', 'Meteosens', 'Fluxnet'
        Fluxnet : ignored
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
        print(f'  Get ICOS product "{product}" with meteo "{meteo}"'
              f' for station "{station}".')
    df, dfunit = read_icos(station, product=product, meteo=meteo,
                           units=True, concat=True)
    icos_units = dict(zip(list(df.columns), dfunit))

    # include units in column names
    ocol = {}
    for cc in df.columns:
        if cc in icos_units:
            ocol.update({cc: f'{cc} ({icos_units[cc]})'})
    if len(ocol) > 0:
        df.rename(columns=ocol, inplace=True)

    # write csv file
    if verbose:
        print(f'  Write into csv file "{outfile}".')
    df.to_csv(outfile, na_rep=undef, index_label='TIMESTAMP_END')

    return
