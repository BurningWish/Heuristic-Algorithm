"""
Do the delaunay triangulation work
"""

import psycopg2
import numpy as np
from shapely.wkt import loads
from shapely.geometry import Polygon
from scipy.spatial import Delaunay


def uidToPotentialSids(dbname, password, user):
    conn = psycopg2.connect("dbname = %s password = %s user = %s" % (dbname, password, user))  # NOQA
    cur = conn.cursor()

    # Get substation points
    cur.execute("SELECT sid, ST_AsText(geom) from substations order by sid asc")  # NOQA
    results = cur.fetchall()

    points = []

    for line in results:
        wkt = line[1]
        sub_point = loads(wkt)
        points.append((sub_point.x, sub_point.y))

    points = np.array(points)

    # Merge the points for triangulation, note "points" contain sub and convex hull  # NOQA
    tri = Delaunay(points)

    cur.execute("drop view if exists unions_to_triangles")
    # Write the triangles to the database
    cur.execute("drop table if exists triangles")
    cur.execute("create table triangles (trid integer, geom geometry(Polygon, 27700), vxid integer[])")  # NOQA
    conn.commit()

    trid = 0
    for indices in tri.simplices:
        triangle = Polygon(points[indices])
        wkt = triangle.wkt
        cur.execute("insert into triangles (trid, geom, vxid) \
                    values(%d, st_geomfromtext('%s', 27700), '{%d, %d, %d}')"
                    % (trid, wkt, indices[0], indices[1], indices[2]))
        trid = trid + 1
    conn.commit()

    # Search which union belongs to which triangle
    cur.execute("create view unions_to_triangles as select u.uid, tr.trid, tr.geom, tr.vxid \
                from unions as u, triangles as tr \
                where st_within(u.cp, tr.geom) \
                order by uid")  # NOQA

    conn.commit()

    cur.execute("select utt.uid, utt.trid, tr.trid, tr.vxid \
                from unions_to_triangles as utt, triangles as tr \
                where st_intersects(utt.geom, tr.geom) \
                order by utt.uid")
    results = cur.fetchall()

    # Note I have to make it as a dictionary instead of a list anymore
    uidToPotentialSids = {}
    for result in results:
        uid = result[0]
        if uid not in uidToPotentialSids:
            # if uid appears for the first time, make a new key for the dict
            uidToPotentialSids[uid] = []
        vertices = result[3]
        for vertice in vertices:
            if vertice not in uidToPotentialSids[uid]:
                uidToPotentialSids[uid].append(vertice)

    cur.close()
    conn.close()

    return uidToPotentialSids
