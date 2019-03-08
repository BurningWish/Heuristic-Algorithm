import networkx as nx
import fiona
import shapely
from shapely.geometry import LineString, Point

def read_shp(edge_path, node_path):
    network = nx.Graph()    
    
    with fiona.open(edge_path, 'r') as c:
        for record in c:
            old_coords = record['geometry']['coordinates']
            new_coords = []
            for coord in old_coords:
                x = round(coord[0], 3)
                y = round(coord[1], 3)
                new_coords.append((x, y))
            new_line = LineString(new_coords)
            startNode = new_coords[0]
            endNode = new_coords[-1]
            network.add_edge(startNode, endNode)
            network.edge[startNode][endNode]['Wkb'] = new_line.wkb
            network.edge[startNode][endNode]['Wkt']= new_line.wkt
            network.edge[startNode][endNode]['feeder'] = record['properties']['feeder']
            network.edge[startNode][endNode]['length'] = record['properties']['Length']
            network.edge[startNode][endNode]['netid'] = record['properties']['netID']
            network.edge[startNode][endNode]['type'] = record['properties']['type']

    with fiona.open(node_path, 'r') as c:
        for record in c:
            point = record['geometry']['coordinates']
            node = (round(point[0], 3), round(point[1], 3))
            network.node[node]['netid'] = record['properties']['netID']
            network.node[node]['toid'] = record['properties']['toid']
            network.node[node]['type']= record['properties']['type']
            
    return network
