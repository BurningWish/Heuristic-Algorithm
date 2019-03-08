"""
After the data (buildings, substations, roads) have been loaded into postgis,
run this script. It will do some data preparations (such as giving each
building a unique id)
"""
import psycopg2

"""
===========================  Parameters for this script =======================
"""
dbname = "sample_hackthon"
password = "19891202"
user = "postgres"
"""
===============================================================================
"""

conn = psycopg2.connect("dbname = %s password = %s user = %s" % (dbname, password, user))  # NOQA
cur = conn.cursor()

# assign a unique id to each building
cur.execute("alter table buildings drop column if exists bid")
cur.execute("alter table buildings add column bid integer")
cur.execute("update buildings set bid = gid - 1")

# calculate the centroid of each building
cur.execute("alter table buildings drop column if exists cp")
cur.execute("alter table buildings add column cp geometry(Point, 27700)")
cur.execute("update buildings set cp = st_centroid(geom)")

# assign a unqiue id to each substation
cur.execute("alter table substations drop column if exists sid")
cur.execute("alter table substations add column sid integer")
cur.execute("update substations set sid = gid - 1")

# assign a unique id to each road
cur.execute("alter table roads drop column if exists rid")
cur.execute("alter table roads add column rid integer")
cur.execute("update roads set rid = gid - 1")

conn.commit()

cur.close()
conn.close()
