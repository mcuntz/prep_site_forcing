#!/usr/bin/env python
'''
Read data from the Fluxnet Shuttle



History
-------
   * Written Jun 2026 by Matthias Cuntz

'''
import datetime as dt
import glob
import os
import zipfile
import pandas as pd
import fluxnet_shuttle as fs
# Need to install Fluxnet shuttle
# python -m pip install git+https://github.com/fluxnet/shuttle.git


__all__ = ['check_fluxnet_snapshot_today', 'info_fluxnet', 'read_fluxnet']


def _download_fluxnet_sites(site, snapshot, download_dir='.'):
    '''
    Download Fluxnet file

    Parameters
    ----------
    site : str
        Fluxnet site id (case-sensitive).
    snapshot : str
        Name of the snapshot file
    download_dir : str, optional
        Directory for downloaded Fluxnet files (default: '.').

    Returns
    -------
    str
        Name of the downloaded file

    '''
    # AmeriFlux intended use codes:
    #   1 = Synthesis
    #   2 = Model
    #   3 = Remote Sensing
    #   4 = Other Research
    #   5 = Education
    #   6 = Other
    user_info = {
        'ameriflux': {
            'user_name': '',    # Optional: your full name
            'user_email': '',   # Optional: your email address
            'intended_use': 2,  # Optional: intended use code (see below)
            'description': 'Using local meteorology to run ecosystem model',
        }
    }

    iamone = False
    if isinstance(site, str):
        iamone = True
        site = [site]

    dfilenames = fs.download(
        site_ids=site,
        snapshot_file=snapshot,
        output_dir=download_dir,
        user_info=user_info)

    if iamone:
        return dfilenames[0]
    else:
        return dfilenames


def check_fluxnet_snapshot_today(snapshot_dir='.'):
    '''
    Returns newest snapshot file of today

    Parameters
    ----------
    snapshot_dir : str, optional
        Directory to look for snapshot file(s) (default: '.').

    Returns
    -------
    str
        Name of the newest snapshot file of today, otherwise empty

    '''
    snapshot = ''
    sfiles = glob.glob(f'{snapshot_dir}/fluxnet_shuttle_snapshot_*.csv')
    if len(sfiles) > 0:
        today = dt.datetime.today().date()
        sdate = today.strftime('%Y%m%d')
        for ss in sfiles:
            bb = os.path.basename(ss)
            if bb.startswith(f'fluxnet_shuttle_snapshot_{sdate}'):
                snapshot = max(snapshot, ss)

    return snapshot


def info_fluxnet(site='', snapshot_dir='.', force_new=False,
                 variables=False, download_dir='.'):
    '''
    Info on available Fluxnet site

    Parameters
    ----------
    site : str, optional
        Fluxnet site id (case-sensitive).
        Print available information on site if given,
        otherwise print available Fluxnet sites.
    snapshot_dir : str, optional
        Directory with snapshot file. A new file will be downloaded
        if no snapshot file is found of the current day (default: '.').
    force_new : bool, optional
        If True, download new snapshot file independent of existsing snapshot
        file.
    variables : bool, optional
        If True, also give info on variables available for site.
        Data will be downloaded if not yet available locally.
    download_dir : str, optional
        Directory for downloaded Fluxnet files (default: '.').

    Returns
    -------
    None

    '''
    # check existing snapshot files
    snapshot = check_fluxnet_snapshot_today(snapshot_dir)
    if (snapshot == '') or force_new:
        # get new snapshot file
        snapshot = fs.listall(output_dir=snapshot_dir)

    snap = pd.read_csv(snapshot)

    if site == '':
        print('Known Fluxnet sites:')
        print(sorted(snap['site_id'].values))

        return
    else:
        print(f'Current snapshot entries for {site}:')
        ss = snap[snap['site_id'] == site]
        for cc in ss.columns:
            print(f'  {cc}: {ss[cc].values[0]}')

        if variables:
            ff = ss['fluxnet_product_name'].values[0]
            dfile = glob.glob(f'{download_dir}/{ff}')
            if len(dfile) == 0:
                dfile = _download_fluxnet_sites(
                    site, snapshot, download_dir=download_dir)
            else:
                dfile = dfile[0]

            with zipfile.ZipFile(dfile) as zfile:
                efile = ''
                for nn in zfile.namelist():
                    if '_FLUXMET_HH_' in nn:
                        if efile == '':
                            efile = nn
                        else:
                            ValueError(f'Several FLUXMET_HH files in {dfile}:\n'
                                       f'{zfile.namelist()}')
                if efile == '':
                    ValueError(f'No FLUXMET_HH file found in {dfile}:\n'
                               f'{zfile.namelist()}')

                with zfile.open(efile) as hfile:
                    df = pd.read_csv(
                        hfile,
                        parse_dates=['TIMESTAMP_START', 'TIMESTAMP_END'],
                        na_values=('-9999', '-9999.0', '-9999.9') )
                
                pvars = list(df.columns)
                pvars.sort()
                print('')
                print(f'Available variables in {dfile}/{efile}:')
                print(', '.join(pvars))

    return


def read_fluxnet(site, snapshot, download_dir='.', units=False):
    '''
    pandas.DataFrame for data of Fluxnet site

    Parameters
    ----------
    site : str
        Fluxnet site id (case-sensitive).
    snapshot : str
        Snapshot file.
    download_dir : str, optional
        Directory for downloaded Fluxnet files (default: '.').
    units : bool, optional
        Return also dictionary with units for variables
        (default: False)

    Returns
    -------
    pandas.DataFrame with Fluxnet data product[, dict of units]

    '''
    # check existing snapshot files
    if not os.path.exists(snapshot):
        raise ValueError(f'Snapshot file {snapshot} must exist.')

    snap = pd.read_csv(snapshot)

    ss = snap[snap['site_id'] == site]
    if len(ss) ==0:
        raise ValueError(f'Site {site} not in snapshot file {snapshot}.')
    
    ff = ss['fluxnet_product_name'].values[0]
    dfile = glob.glob(f'{download_dir}/{ff}')
    if len(dfile) == 0:
        dfile = _download_fluxnet_sites(
            site, snapshot, download_dir=download_dir)
    else:
        dfile = dfile[0]

    with zipfile.ZipFile(dfile) as zfile:
        efile = ''
        for nn in zfile.namelist():
            if '_FLUXMET_HH_' in nn:
                if efile == '':
                    efile = nn
                else:
                    ValueError(f'Several FLUXMET_HH files in {dfile}:\n'
                               f'{zfile.namelist()}')
        if efile == '':
            ValueError(f'No FLUXMET_HH file found in {dfile}:\n'
                       f'{zfile.namelist()}')

        with zfile.open(efile) as hfile:
            df = pd.read_csv(
                hfile,
                parse_dates=['TIMESTAMP_START', 'TIMESTAMP_END'],
                na_values=('-9999', '-9999.0', '-9999.9') )
            df.set_index('TIMESTAMP_END', inplace=True)
            df.index.name = 'Date Time'

        if units:
            bfile = efile.replace('FLUXMET', 'BIFVARINFO')
            with zfile.open(bfile) as bbfile:
                meta = pd.read_csv(bbfile)
            gids = meta['GROUP_ID'].unique()
            vunits = {}
            for gid in gids:
                vv = meta[meta['GROUP_ID'] == gid]
                ivar = vv.loc[vv['VARIABLE'] == 'VAR_INFO_VARNAME', 'DATAVALUE']
                iunit = vv.loc[vv['VARIABLE'] == 'VAR_INFO_UNIT', 'DATAVALUE']
                if (len(ivar) > 0) and (len(iunit) > 0):
                    vvar = ivar.values[0]
                    vunit = iunit.values[0]
                    vunit = vunit.replace('�', 'µ')  # strange encoding of BIF
                    if vvar in vunits:
                        if vunit != vunits[vvar]:
                            raise ValueError(
                                f'Units do not match for {vvar}:'
                                f' {vunits[vvar]} != {vunit}')
                    vunits.update({vvar: vunit})

    if units:
        return df, vunits
    else:
        return df


if __name__ == '__main__':

    import argparse
    import os

    download_dir = '.'
    force_new_snapshot = False
    info = False
    snapshot_dir = '.'
    site = ''
    variables = False

    desc = 'Info and retrieval of Fluxnet data.'
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=desc)

    hstr = 'Directory for downloaded Fluxnet files (default: ".")'
    parser.add_argument('-d', '--download_dir', action='store',
                        default=download_dir, dest='download_dir',
                        metavar='directory', help=hstr)

    hstr = 'Force download of a new snapshot.'
    parser.add_argument('-f', '--force_new_snapshot', action='store_true',
                        default=force_new_snapshot, dest='force_new_snapshot',
                        help=hstr)

    hstr = 'Print info in snapshot about Fluxnet site.'
    parser.add_argument('-i', '--info', action='store_true',
                        default=info, dest='info',
                        help=hstr)

    hstr = 'Directory with snapshot file(s) (default: ".")'
    parser.add_argument('-s', '--snapshot_dir', action='store',
                        default=snapshot_dir, dest='snapshot_dir',
                        metavar='directory', help=hstr)

    hstr = ('Fluxnet site ID such as FR-Hes (case-sensitive).'
            ' Prints info for site (-i, -v) or downloads site data.'
            ' If not given, prints available Fluxnet sites.')
    parser.add_argument('site', nargs='?', default=site,
                        metavar='site_id', help=hstr)

    hstr = ('Print info about available variables at Fluxnet site.'
            ' Data will be downloaded if not yet on disk.')
    parser.add_argument('-v', '--variables', action='store_true',
                        default=variables, dest='variables',
                        help=hstr)

    args = parser.parse_args()

    download_dir = args.download_dir
    force_new_snapshot = args.force_new_snapshot
    info = args.info
    snapshot_dir = args.snapshot_dir
    site = args.site
    variables = args.variables

    # info, variables, or no site id
    if info or variables or (site == ''):
        info_fluxnet(site=site,
                     snapshot_dir=snapshot_dir,
                     force_new=force_new_snapshot,
                     variables=variables,
                     download_dir=download_dir)
    else:
        snapshot = check_fluxnet_snapshot_today(snapshot_dir)
        if (snapshot == '') or force_new_snapshot:
            # get new snapshot file
            snapshot = fs.listall(output_dir=snapshot_dir)

        snap = pd.read_csv(snapshot)

        ss = snap[snap['site_id'] == site]
        if len(ss) ==0:
            raise ValueError(f'Site {site} not in snapshot file {snapshot}.')
        
        ff = ss['fluxnet_product_name'].values[0]
        dfile = glob.glob(f'{download_dir}/{ff}')
        if len(dfile) == 0:
            dfile = _download_fluxnet_sites(
                site, snapshot, download_dir=download_dir)
        else:
            dfile = dfile[0]

        print(f'Fluxnet data file for {site}: {dfile}')
