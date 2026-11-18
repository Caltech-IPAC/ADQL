import pytest

from ADQL.adql import ADQL
from spatial_index import SpatialIndex


def translate(adql, dbms='oracle'):
    return ADQL(dbms=dbms, level=20, racol='ra', deccol='dec',
                xcol='x', ycol='y', zcol='z', indxcol='htm20',
                mode=SpatialIndex.HTM, encoding=SpatialIndex.BASE10).sql(adql)


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
