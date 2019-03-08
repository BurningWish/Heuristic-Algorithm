"""
In the last script, we generated terraces,

this script will further combine terraces into unions,

if terraces are very close with each other (for example, distance < 5 meters),

as a result, this script will generate a table in the postgis database called unions
"""

from shapely.ops import cascaded_union
from shapely.wkt import loads
import psycopg2
import networkx as nx
from shapely.geometry import MultiPolygon

"""
========================== Parameters of this script  =========================
"""
dbname = "sample_hackthon"
password = "19891202"
user = "postgres"
threshold_distance = 5  # unit is meter. For example, it means if two or more terraces are within 5 meters with each other, they will produce a union
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

# Be careful, that every element in the c_union is a Polygon object

# Let's say 5 meters is the threshold now
cur.execute("select t1.tid, t2.tid, st_distance(t1.geom, t2.geom) as dist \
            from terraces as t1, terraces as t2 \
            where t1.tid < t2.tid \
            and st_dwithin(t1.geom, t2.geom, %d) \
            order by t1.tid, dist" % threshold_distance)
results = cur.fetchall()
print("We use %d meters as the threshold distance" % threshold_distance)

# this solution of buffered union finding solution is GENIUS!!!!!!!!!!!!!!!!!!!
group_graph = nx.Graph()
for result in results:
    group_graph.add_edge(result[0], result[1])

graph_list = list(nx.connected_component_subgraphs(group_graph))

# now the remaining thing is much straight forward
# create use a table to store the buffered union

buffered_tids = []

for graph in graph_list:
    for tid in graph.nodes():
        buffered_tids.append(tid)

buffered_tids.sort()
tup_buf_tids = (tuple(buffered_tids),)

cur.execute("drop view if exists unions_to_triangles")
cur.execute("drop table if exists unions")
# we want accurate spatial geometry of unions, so let's still use multipolygon instead of convex hull  # NOQA
cur.execute("create table unions (uid integer, geom geometry(MultiPolygon, 27700))")  # NOQA
conn.commit()

uid = 0

# first let's copy the normal terraces to the unions table
cur.execute("select tid, st_astext(geom) from terraces \
            where tid not in %s order by tid" % tup_buf_tids)
terraces_results = cur.fetchall()
for terrace_result in terraces_results:
    wkt = terrace_result[1]
    fake_polys = list()
    fake_polys.append(loads(wkt))
    multi_poly = MultiPolygon(fake_polys)
    new_wkt = multi_poly.wkt
    cur.execute("insert into unions (uid, geom) \
                values (%d, st_geomfromtext('%s', 27700))" % (uid, new_wkt))
    uid += 1
    # print("copy terrace into unions table, %d" % uid)
conn.commit()


# now the hard part, generate the true unions and save it to the unions table
for graph in graph_list:
    current_poly_list = []

    for tid in graph.nodes():
        cur.execute("select st_astext(geom) from terraces \
                    where tid = %d" % tid)
        old_wkt = cur.fetchone()[0]
        # remember terraces table is polygon geometry
        temp_poly = loads(old_wkt)
        current_poly_list.append(temp_poly)

    current_multi_poly = MultiPolygon(current_poly_list)
    wkt = current_multi_poly.wkt

    cur.execute("insert into unions (uid, geom) \
                values (%d, st_geomfromtext('%s', 27700))" % (uid, wkt))

    uid += 1
    # print("generate union", uid)
conn.commit()

cur.execute("alter table unions add column cp geometry(Point, 27700)")
cur.execute("update unions set cp = st_centroid(geom)")
conn.commit()

cur.close()
conn.close()
