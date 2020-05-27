import pytest

from adql import ADQL
from spatial_index import SpatialIndex


infile  = open('tests/indata.txt',  'r')
outfile = open('tests/outdata.txt', 'r')
 
values = []

while True:

    adql_in = infile.readline().rstrip('\n')

    if not adql_in:
        break

    sql_out = outfile.readline().rstrip('\n')
    
    if not sql_out:
        break

    values.append((adql_in, sql_out))



@pytest.mark.parametrize("adqlstr,sqlstr", values)

def test_adql(adqlstr,sqlstr):

    adql = ADQL(level=20, debugfile=None,
                racol='ra', deccol='dec',
                xcol = 'x', ycol='y', zcol='z', indxcol='htm20',
                mode=SpatialIndex.HTM, encoding=SpatialIndex.BASE10)

    try:
        sqlnew = adql.sql(adqlstr)

    except Exception as adql_error:

        sqlnew = str(adql_error)

    assert sqlnew == sqlstr
