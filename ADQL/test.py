#!/usr/bin/env python

from adql import ADQL
from spatial_index import SpatialIndex
  
 
# Set up the ADQL statement
 
adql_string = []

# adql_string.append("select TOP 100 psc.ra, psc.dec distance(point('icrs', psc.ra, psc.dec), point('GALACTIC', 234.56, 34.567)) as dist from iraspsc as psc where contains(point('icrs', psc.ra, psc.dec), circle('GALACTIC', 234.56, 34.567, 0.006)) = 1 and psc.glat > 34.567 order by dec desc")
adql_string.append("select ra, dec from iraspsc where contains(point('icrs', ra, dec), box('ECLIPTIC', 233.56, 34.567, 1., 2.)) = 1 order by dec desc")


for i in range(len(adql_string)):

    print('------------------------------------')
    print('')
    print('BEFORE: ', adql_string[i])
    print('')

    adql = ADQL(level=20, debugfile='test.debug',
                racol='ra', deccol='dec',
                xcol = 'x', ycol='y', zcol='z', indxcol='htm20',
                mode=SpatialIndex.HTM, encoding=SpatialIndex.BASE10)

    try:
        sql_string = adql.sql(adql_string[i])

        print('AFTER:  ', sql_string)
        print('')

    except Exception as adql_error:
        print('ERROR: ', adql_error)
