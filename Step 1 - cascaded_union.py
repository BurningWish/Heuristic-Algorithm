"""
We do the cascaded_union operation on the buildings,

i.e. buildings are merged if they are topologically connected with each other,

as a result, in the postgis database, we will create a table called terraces
"""

from shapely.ops import cascaded_union
from shapely.wkt import loads
import psycopg2

"""
=======================  Parameters for this script   =========================
"""
dbname = "sample_hackthon"
password = "19891202"
user = "postgres"
"""
===============================================================================
"""

building_polys = []
conn = psycopg2.connect("dbname = %s password = %s user = %s" % (dbname, password, user))  # NOQA
cur = conn.cursor()
cur.execute("SELECT bid, ST_AsText(geom) from buildings order by bid asc")
results = cur.fetchall()

for line in results:
    coordsData = line[1][13: -1]
    wkt = 'POLYGON ' + coordsData
    polygon = loads(wkt)
    building_polys.append(polygon)

c_union = cascaded_union(building_polys)

cur.execute("drop table if exists terraces")
cur.execute("create table terraces (tid integer, geom geometry(Polygon, 27700))")  # NOQA
tid = 0
for terrace in c_union:
    wkt = terrace.wkt
    cur.execute("insert into terraces (tid, geom) \
                values (%d, st_geomfromtext('%s', 27700))" % (tid, wkt))
    tid += 1

cur.execute("alter table terraces add column cp geometry(Point, 27700)")
cur.execute("update terraces set cp = st_centroid(geom)")

conn.commit()
cur.close()
conn.close()
