import math
import sqlite3

import pytest

from ADQL.adql import ADQL
from spatial_index import SpatialIndex


def translate(adql, dbms='oracle'):
    return ADQL(dbms=dbms, level=20, racol='ra', deccol='dec',
                xcol='x', ycol='y', zcol='z', indxcol='htm20',
                mode=SpatialIndex.HTM, encoding=SpatialIndex.BASE10).sql(adql)


def xyz(ra, dec):
    ra, dec = math.radians(ra), math.radians(dec)
    return math.cos(ra) * math.cos(dec), math.sin(ra) * math.cos(dec), math.sin(dec)


@pytest.mark.parametrize("adql,sql", [
    ("select top 3 pl_name from pscomppars where pl_orbper is not null order by pl_orbper desc",
     "select pl_name from pscomppars where pl_orbper is not null order by pl_orbper desc FETCH FIRST 3 ROWS ONLY"),
    ("select top 1 discoverymethod, count(*) from pscomppars group by discoverymethod order by count(*) desc",
     "select discoverymethod, count(*) from pscomppars group by discoverymethod order by count(*) desc FETCH FIRST 1 ROWS ONLY"),
    ("select top 1 count(*) from pscomppars",
     "select count(*) from pscomppars FETCH FIRST 1 ROWS ONLY"),
])
def test_top(adql, sql):
    assert translate(adql) == sql


@pytest.mark.parametrize("region,inside,outside", [
    ("circle('icrs', 10., 20., 5.)", [(10., 20.)], [(10., 26.), (190., -20.)]),
    ("polygon('icrs', 10., 10., 14., 10., 14., 14., 10., 14.)", [(12., 12.)], [(15., 12.), (12., 9.)]),
    ("box('icrs', 30., 40., 2., 2.)", [(30., 40.)], [(30., 42.), (32.5, 40.)]),
])
def test_contains_zero(region, inside, outside):
    db = sqlite3.connect(':memory:')
    db.execute('create table pts (ra, dec, x, y, z)')
    db.executemany('insert into pts values (?, ?, ?, ?, ?)', [p + xyz(*p) for p in inside + outside])
    sql = translate("select ra, dec from pts where contains(point('icrs', ra, dec), " + region + ") = 0", 'sqlite3')
    assert sorted(db.execute(sql)) == sorted(outside)


def test_distance_fixed_point_first():
    query = "select distance({}, {}) as d from pscomppars"
    column, fixed = "point('icrs', ra, dec)", "point('icrs', 290.67, 44.5)"
    assert translate(query.format(fixed, column)) == translate(query.format(column, fixed))
