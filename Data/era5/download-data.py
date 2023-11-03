##Script para solicitud MARS de ERA-5
#LINK:   https://apps.ecmwf.int/data-catalogues/era5/?type=an&class=ea&stream=oper&expver=1
# Select the variables, conditions, etc, ask the MARS request.  Copy the code of the CDS API script details, and reeplace the code below
# OR, if you understand the code:
#  Modify the code referring to your interest: grid, format and area
#  'area'    : '80/-50/-25/0',          # North, West, South, East. If yo want global, you need to put 90/-180/-90/180
# Include the directory path in the ouput filename.
# 
###############################  
###############################   Surface e Invariant
#!/usr/bin/env python


import cdsapi

c = cdsapi.Client()

c.retrieve("reanalysis-era5-complete", {
    "class": "ea",
    "dataset": "era5",
    "date": "2022-08-07/2022-08-08",
    "expver": "1",
    "levtype": "sfc",
    "param": "129.128/159.128/165.128/166.128/167.128/172.128",
    "stream": "oper",
    "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
    "type": "an"
    'area'    : '80/-50/-25/0',          # North, West, South, East. Default: global
    'grid'    : '1.0/1.0',               # Latitude/longitude. Default: spherical harmonics or reduced Gaussian grid
    'format'  : 'netcdf',                # Output needs to be regular lat-lon, so only works in combination with 'grid'!
}, "output.nc")
#
#
###############################   Pressure levels
#!/usr/bin/env python
import cdsapi

c = cdsapi.Client()

c.retrieve("reanalysis-era5-complete", {
    "class": "ea",
    "dataset": "era5",
    "date": "2022-08-07/2022-08-08",
    "expver": "1",
    "levelist": "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    "levtype": "pl",
    "param": "129.128/130.128/131/132/135.128/157.128",
    "stream": "oper",
    "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
    "type": "an"
    'area'    : '80/-50/-25/0',          # North, West, South, East. Default: global
    'grid'    : '1.0/1.0',               # Latitude/longitude. Default: spherical harmonics or reduced Gaussian grid
    'format'  : 'netcdf',                # Output needs to be regular lat-lon, so only works in combination with 'grid'!
}, "outputPress.nc")
