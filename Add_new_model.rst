Add new model
=============

Example of how to add a new ecosystem model.


Content
-------

Preparation of forcing files uses the classes ``prepSiteForcing`` and
``ascii2Netcdf``. Both need to be enhanced. The easiest is to start
with ``ascii2musica.py`` for the latter. The former can be enhanced
pragmatically but there is very little model specific code in
``class_prep_site_forcing.py`` so the easiest is to change the file
directly. We will be happy to receive Pull Requests (PR) that include
new models.

Steps are:

* `Make ascii2netcdf script`_
* `Enhance prepSiteForcing`_
* `Adapt config file`_

Here we add the land surface model ORCHIDEE_ as an example.


Make ascii2netcdf script
------------------------

The netcdf forcing file for ORCHIDEE has the following structure: ::

   netcdf FR-Fon_2005-2020 {
   dimensions:
           tstep = UNLIMITED ; // (280512 currently)
           lev = 1 ;
           bound = 2 ;
           lon = 1 ;
           lat = 1 ;
   variables:
           double tstep(tstep) ;
                   tstep:origin = "2005-01-01 00:00:00" ;
                   tstep:units = "seconds since 2005-01-01 00:00:00" ;
                   tstep:calendar = "gregorian" ;
           float bounds_lon(lon, bound) ;
           float bounds_lat(lat, bound) ;
           double lat(lat) ;
                   lat:units = "degrees_north" ;
                   lat:long_name = "Latitude" ;
           double lon(lon) ;
                   lon:units = "degrees_east" ;
                   lon:long_name = "Longitude" ;
           double lev(lev) ;
                   lev:units = "m" ;
                   lev:long_name = "Vertical levels" ;
           double SWdown(tstep, lat, lon) ;
                   SWdown:units = "W/m^2" ;
                   SWdown:long_name = "Surface incident shortwave radiation" ;
           double LWdown(tstep, lat, lon) ;
                   LWdown:units = "W/m^2" ;
                   LWdown:long_name = "Surface incident longwave radiation" ;
           double Rainf(tstep, lat, lon) ;
                   Rainf:units = "kg/m^2s" ;
                   Rainf:long_name = "Rainfall rate" ;
           double Snowf(tstep, lat, lon) ;
                   Snowf:units = "kg/m^2s" ;
                   Snowf:long_name = "Snowfall rate" ;
           double PSurf(tstep, lat, lon) ;
                   PSurf:units = "Pa" ;
                   PSurf:long_name = "Surface pressure" ;
           double Qair(tstep, lev, lat, lon) ;
                   Qair:units = "kg/kg" ;
                   Qair:long_name = "Near Surface specific humidity" ;
           double Tair(tstep, lev, lat, lon) ;
                   Tair:units = "K" ;
                   Tair:long_name = "Near surface air temperature" ;
           double Wind_N(tstep, lev, lat, lon) ;
                   Wind_N:units = "m/s" ;
                   Wind_N:long_name = "Near surface northward wind component" ;
           double Wind_E(tstep, lev, lat, lon) ;
                   Wind_E:units = "m/s" ;
                   Wind_E:long_name = "Near surface eastward wind component" ;
           double CO2_f(tstep, lat, lon) ;
                   CO2_f:units = "umolCO2 mol-1" ;
                   CO2_f:long_name = "CO2 mole fraction" ;

   // global attributes:
                   :description = "Half-hourly data for site FR-Fon" ;

   data:

    tstep = 0, 1800, ...

    bounds_lon =
     2.2801, 3.2801 ;

    bounds_lat =
     47.9764, 48.9764 ;

    lat = 48.4764 ;

    lon = 2.7801 ;

    lev = 2 ;

    SWdown = ...

    ...
   }

One needs a function ``ascii2orchidee`` within a file
``ascii2orchidee.py`` that can be called from within the method
``write_netcdf`` of the class ``prepSiteForcing``. It is convenient
that this script can also be called from the command line with command
line options (see below).

We enhance the class ``ascii2Netcdf``. We have to set the name of the
output variables, convert the variables to the correct units, and
write the netcdf file.

* We copy the file ``ascii2musica.py`` to ``ascii2orchidee.py`` for
  a start.
* In the file, replace all instances of *musica* by *orchidee* (I
  would follow the same uppercase/lowercase conventions as used in
  ``ascii2musica``).
* One can delete all code that uses *site_altitude* and
  *namelist*. That's in ``__init__``, ``ascii2orchidee``, the command
  line options and the docstrings.
* Replace all *site_* with nothing ''. That changes all
  *site_latitude* to *latitude* and all *site_longitude* to
  *longitude*.
* The timestep in ORCHIDEE seems to indicate the start of the
  timestep. Change the two 0.5 in the block ``if self.ftimestep !=
  1.0:`` in ``write_netcdf`` to 0.0. This will shift all model forcing
  time steps to the beginning of the timestep.

The ``__init__`` method should now look like this:

.. code-block:: python

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

Then we change the method ``output_variables`` for the correct names and
units:

* Change *CO2air* to *CO2_f* in ``self.co2air``.
* Change *Near surface CO2 concentration* to *CO2 mole fraction*,
  *Near surface northerly wind speed* to *Near surface northward
  wind component*, and *Near surface easterly wind speed* to *Near
  surface eastward wind component* in ``self.varlongnames``.
* We will not change the units to the non-standard units of ORCHIDEE.
     
The ``output_variables`` method now looks like this:

.. code-block:: python

   def output_variables(self):
       """
       Output variable names, longnames, and units

       """
       # name of downwelling shortwave radiation 
       # names of standard variables
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
       # output variables and units
       self.stdnames = ['co2air', 'swdown', 'lwdown',
                        'psurf', 'qair', 'rainf',
                        'snowf', 'tair', 'wind1', 'wind2']
       self.varnames = [self.co2air, self.swdown, self.lwdown,
                        self.psurf, self.qair, self.rainf,
                        self.snowf, self.tair, self.wind1, self.wind2]  # ,
                        # 'h_sbl']
       self.varnames_lower = [ vv.lower() for vv in self.varnames ]
       self.varunits = ['ppmv', 'W/m2', 'W/m2', 'Pa',
                        'kg/kg', 'kg/m2/s', 'kg/m2/s', 'K',
                        'm/s', 'm/s']  # ,
                        # 'm']
       self.varlongnames = ['Near surface CO2 concentration',
                            'Surface incident shortwave radiation',
                            'Surface incident longwave radiation',
                            'Surface pressure',
                            'Near surface specific humidity',
                            'Rainfall rate',
                            'Snowfall rate',
                            'Near surface air temperature',
                            'Near surface northerly wind speed',
                            'Near surface easterly wind speed']  # ,
                            # 'Surface boundary layer height']

       # dictionary with lowercase variables as keys
       self.varnames_lower = [ vv.lower() for vv in self.varnames ]
       self.dict_varunits = dict(zip(self.varnames_lower, self.varunits))
       self.dict_varnames = dict(zip(self.varnames_lower, self.varnames))
       self.dict_varlongnames = dict(zip(self.varnames_lower,
                                         self.varlongnames))
       self.dict_stdnames = dict(zip(self.stdnames, self.varnames_lower))

The ``convert_units_orchidee`` does already exactly what we want. It
sets the correct output variable names, and deals with the wind
components *Wind_N* and *Wind_E* if wind speed and wind direction are
given.

.. code-block:: python

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

Last, we adapt the method ``write_netcdf``. First we set and write the
dimensions:

.. code-block:: python

   lon  = fo.createDimension('lon', 1)
   lat  = fo.createDimension('lat', 1)
   lev  = fo.createDimension('lev', 1)
   tstep = fo.createDimension('tstep', None)
   bound = fo.createDimension('bound', 2)

Then we define the dimension variables:

.. code-block:: python

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

Then we define the variables. Some forcing variables have a lev
dimension, such Tair, and some not. ORCHIDEE allows to have other variables in the file such as heat fluxes for model-data comparison, for example. So we just take the code from ``ascii2musica``:

.. code-block:: python

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

Last, we write the values into the variables. Times can either be
given in local time or in UTC in ORCHIDEE. We shift the times to
UTZ if the option time2gmt is different to 0.

.. code-block:: python

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

The final ``write_netcdf`` method looks like this:

.. code-block:: python

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


Enhance prepSiteForcing
-----------------------

Basically, one has to specify the options for the call to
``ascii2orchidee`` in the method ``write_netcdf`` in the class
``prepSiteForcing``. That's pretty much all, so we write it just into
``class_prep_site_forcing.py`` and do not enhance the class via
another class. You can search for 'isba' in the file to see the few
places the specific model is used. There is only a check for
*_FillValue* apart from the options in ``write_netcdf``.

We just copy the section ``if leco == 'musica':`` at the end of the if
block and change it to ``elif leco == 'orchidee':``. In this block, we
remove again the two lines with *altitude* and *self_altitude*, and
remove *site_* from *longitude* and *latitude*. This gives the block:

.. code-block:: python

   elif leco == 'orchidee':
       kwargs.update({
           'latitude': self.latitude,
           'longitude': self.longitude,
           'time2gmt': self.time2gmt,
           'forcing_height': self.reference_height})
       aopt += (f' --latitude {self.latitude}'
                f' --longitude {self.longitude}'
                f' --time2gmt {self.time2gmt}'
                f' --forcing_height {self.reference_height}')
       print(f'Calling python ascii2{leco}.py {aopt} {csvfile}')

The *kwargs* are for the call to ``ascii2orchidee``. The *aopt* are only
in the print. It prints the command line that one has type to make the same conversion of the ascii file to the netcdf file. This is practical, for example, if one wants impute long time series but later have shorter forcing files. Or one wants to add extra columns to the imputed ascii file after prep_site_forcing and run only the conversion ascii2netcdf again.


Adapt config file
-----------------

If you added a new variable or option for a new model, you can add it
to the section [Model Options] of the config file. Do not forget to
read it in the method ``read_configfile`` in the class
``prepSiteForcing``. Here we just "recycled" the option time2gmt for
ORCHIDEE that was already present for MuSICA.


.. _ORCHIDEE: https://orchidee.ipsl.fr
