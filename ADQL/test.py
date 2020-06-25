#!/usr/bin/env python

from adql import ADQL
from spatial_index import SpatialIndex
  
 
# Set up the ADQL statement
 
adql_string = []

adql_string.append("select TOP 100 psc.ra, psc.dec distance(point(\"\", psc.ra, psc.dec), point('GALACTIC', 234.56, 34.567)) as dist from iraspsc as psc where contains(point('icrs', psc.ra, psc.dec), circle('GALACTIC', 234.56, 34.567, 0.006)) = 1 and psc.glat > 34.567 order    by dec desc group   \n  by dec")

for i in range(len(adql_string)):

    print('------------------------------------')
    print('')
    print('BEFORE: ', adql_string[i])
    print('')

    adql = ADQL(dbms='oracle', level=20, debugfile='test.debug',
                racol='RA', deccol='Dec',
                xcol = 'X', ycol='Y', zcol='Z', indxcol='HTM20',
                mode=SpatialIndex.HTM, encoding=SpatialIndex.BASE10)

    try:
        sql_string = adql.sql(adql_string[i])

        print('AFTER:  ', sql_string)
        print('')

    except Exception as adql_error:
        print('ERROR: ', adql_error)

for i in range(len(adql_string)):

    print('------------------------------------')
    print('')
    print('BEFORE: ', adql_string[i])
    print('')

    adql = ADQL(dbms='sqlite3', level=20, debugfile='test.debug',
                racol='RA', deccol='Dec',
                xcol = 'X', ycol='Y', zcol='Z', indxcol='HTM20',
                mode=SpatialIndex.HTM, encoding=SpatialIndex.BASE10)

    try:
        sql_string = adql.sql(adql_string[i])

        print('AFTER:  ', sql_string)
        print('')

    except Exception as adql_error:
        print('ERROR: ', adql_error)
