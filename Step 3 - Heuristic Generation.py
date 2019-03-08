"""
In the last step, we generated a table called unions,

now we can generate the heuristic network, this is what this script does,


"""

import sys
import psycopg2
import pg_read
import copy
import networkx as nx
import fiona
from shapely.geometry import LineString, Point
from math import acos, degrees
import Triangulation
import pickle
import os


"""
=========================  Parameters for this script  ========================
"""
dbname = "sample_hackthon"
password = "19891202"
user = "postgres"
"""
===============================================================================
"""

def finalCheckValidity(distriG):
    nodeSubAc = 0    
    
    for node in distriG:
        if distriG.node[node]['type'] == 'substationAccess':
            nodeSubAc = node
    
    for node in distriG:
        if distriG.node[node]['type'] == 'distribution':
            if not nx.has_path(distriG, nodeSubAc, node):
                distriG.add_edge(nodeSubAc, node)
                distriG.edge[nodeSubAc][node]['Length'] = LineString([nodeSubAc, node]).length
                coords1 = distriG.node[nodeSubAc]['Coords']
                coords2 = distriG.node[node]['Coords']
                distriG.edge[nodeSubAc][node]['Coords'] = [coords1, coords2]
                distriG.edge[nodeSubAc][node]['type'] = 'feeder'

def identifyFeeder(distriG):
    """
    This function is used to identify feeders given a distribution network
    """
    
    """
    First remove the edges which are service lines,
    in other words, we keep the feeders only
    """
    servEdges = []
    rmServLines = []
    for edge in distriG.edges():
        startNode = edge[0]
        endNode = edge[1]
        if distriG.edge[startNode][endNode]['type'] == 'service':
            servEdges.append(edge)
    
    for edge in servEdges:
        startNode = edge[0]
        endNode = edge[1]
        edgeDict = distriG.edge[startNode][endNode]
        rmServLines.append([edge, edgeDict])
        distriG.remove_edge(startNode, endNode)
    
    """
    Then add a new field to the feeder edges
    """    
    for edge in distriG.edges():
        startNode = edge[0]
        endNode = edge[1]
        distriG.edge[startNode][endNode]['feeder'] = []
    
    """
    Now find all the non-substation nodes, whose degree are 1
    """
    destiNodes = []
    subNode = (0,0)
    
    # find the substation node on the distribution network
    for node in distriG.nodes():
        if distriG.node[node]['type'] == 'substation':
            subNode = node
            break
    
    for node in distriG.nodes():
        if node != subNode and distriG.degree(node) == 1 and nx.has_path(distriG, subNode, node):
            destiNodes.append(node)
    
    """
    For each destiNode, resolve the dijkPath, and assign feeder index to edge on the path
    """
    for i in range(len(destiNodes)):
        thisNode = destiNodes[i]
        thisPath = nx.dijkstra_path(distriG, subNode, thisNode, weight='Length')  # NOQA
        for index in range(len(thisPath) - 1):
            startNode = thisPath[index]
            endNode = thisPath[index + 1]
            distriG.edge[startNode][endNode]['feeder'].append(i)
            
    """
    Finally add the removed edges back to the distriG
    """
    for row in rmServLines:
        edge = row[0]
        edgeDict = row[1]
        startNode = edge[0]
        endNode = edge[1]
        distriG.add_edge(startNode, endNode)
        distriG.edge[startNode][endNode] = edgeDict
        distriG.edge[startNode][endNode]['feeder'] = "NULL"

def subAccessPoint(substation, road):
    """
    given a substation and the nearest road to the substation,

    the function will return the best access point on the road

    for the substation. this access point can possibily be

    the projection point or the middle point on the road.
    """
    reference = road.geom.project(substation.geom)
    originAccessPoint = road.geom.interpolate(reference)
    PointS = (substation.geom.x, substation.geom.y)
    PointOAP = (originAccessPoint.x, originAccessPoint.y)
    originAccessLine = LineString((PointS, PointOAP))
    cur.execute("select * from buildings \
                where st_intersects(geom, st_geomfromtext('%s', 27700))" % originAccessLine.wkt)  # NOQA
    originIntersects = len(cur.fetchall())
    if originIntersects == 0:
        # which means the projection point doesn't cause any intersection
        return originAccessPoint
    else:
        # projection point causes intersection, need further decision
        newReference = road.geom.length * 0.5
        newAccessPoint = road.geom.interpolate(newReference)
        PointNAP = (newAccessPoint.x, newAccessPoint.y)
        newAccessLine = LineString((PointS, PointNAP))
        cur.execute("select * from buildings \
                    where st_intersects(geom, st_geomfromtext('%s', 27700))" % newAccessLine.wkt)  # NOQA
        newIntersects = len(cur.fetchall())
        if newIntersects < originIntersects:
            # which means the newAccessPoint is better
            return newAccessPoint
        else:
            # we still use the originAccessPoint
            return originAccessPoint


def simpleObjPickRoad(obj, roads):
    # in here the obj (either union or terrace) consists of one building
    """
    the function will return the rid and accessPoint of road for this obj

    we don't have angle here, and we want the deriving edge don't intersect
    other objs
    """
    fittestRid = -1
    accessPoint = Point(0, 0)

    findRoad = False

    for road in roads:
        reference = road.geom.project(obj.centroid)
        tempAccessPoint = road.geom.interpolate(reference)
        PointC = (obj.centroid.x, obj.centroid.y)
        PointD = (tempAccessPoint.x, tempAccessPoint.y)
        lineCD = LineString((PointC, PointD))

        if type(obj) == pg_read.Union:
            # now we are using a union
            uid = obj.id
            cur.execute("select * from unions \
                        where st_intersects(geom, st_geomfromtext('%s', 27700)) \
                        and uid != %d" % (lineCD.wkt, uid))  # NOQA
            results = cur.fetchall()
            if not results:
                # which means no other unions intersects lineCD
                findRoad = True
                fittestRid = road.id
                accessPoint = tempAccessPoint
                break
        else:
            # which means we are using a terrace
            tid = obj.id
            cur.execute("select * from terraces \
                        where st_intersects(geom, st_geomfromtext('%s', 27700)) \
                        and tid != %d" % (lineCD.wkt, tid))  # NOQA
            results = cur.fetchall()
            if not results:
                # which means no other terraces intersects lineCD
                findRoad = True
                fittestRid = road.id
                accessPoint = tempAccessPoint
                break

    if findRoad:
        # if findRoad == True:
        return fittestRid, accessPoint
    else:
        # which means findRoad == False
        # I know it's sad, but we need to have default option here
        # we use the roads[0] as accessRoad anyway
        road = roads[0]
        reference = road.geom.project(obj.centroid)
        fittestRid = road.id
        accessPoint = road.geom.interpolate(reference)
        return fittestRid, accessPoint


def pickFittestRoad(obj, roads):
    """
    We say union or terrace can choose fittest road from 3 potential roads

    The function will return the rid and accesspoint for the fittest road
    """
    fittestRid = -1
    accessPoint = Point(0, 0)

    PointA, PointB = objDiameter(obj)

    findRoad = False

    for road in roads:
        reference = road.geom.project(obj.centroid)
        tempAccessPoint = road.geom.interpolate(reference)
        PointC = (obj.centroid.x, obj.centroid.y)
        PointD = (tempAccessPoint.x, tempAccessPoint.y)
        deltaX = PointC[0] - PointA[0]
        deltaY = PointC[1] - PointA[1]
        PointE = (PointD[0] - deltaX, PointD[1] - deltaY)
        sideA = LineString((PointA, PointB)).length
        sideB = LineString((PointA, PointE)).length
        sideC = LineString((PointB, PointE)).length

        angle = getAngle(sideA, sideB, sideC)
        if angle > 90:
            angle = 180 - angle

        if angle > 30:
            # we think this angle is large enough
            # one more check, this straightLine CD had better not intersect another obj  # NOQA
            lineCD = LineString((PointC, PointD))
            tid = obj.id
            cur.execute("select * from terraces \
                        where st_intersects(geom, st_geomfromtext('%s', 27700)) \
                        and tid != %d" % (lineCD.wkt, tid))  # NOQA
            results = cur.fetchall()
            if not results:
                # which means no other terraces intersects lineCD
                findRoad = True
                fittestRid = road.id
                accessPoint = tempAccessPoint
                break

    if findRoad:
        return fittestRid, accessPoint

    else:
        # which means findRoad == False
        # we use the middle point of the roads[0] as access point
        terraceList[obj.id].projectType = 'special'
        road = roads[0]
        reference = road.geom.length * 0.5
        accessPoint = road.geom.interpolate(reference)
        fittestRid = road.id
        return fittestRid, accessPoint


def getAngle(A, B, C):
    """
    This function will return the angle of side a and b,
    if you know the length of a, b, c
    """
    if A * B == 0:
        return 180
    else:
        return degrees(acos((A * A + B * B - C * C)/(2.0 * A * B)))


def objDiameter(obj):
    """
    This function return a coordinate (PointA, PointB)
    where PointA is (Ax, Ay) and PointB is (Bx, By)
    """
    bids = obj.buildingIdList
    if len(bids) == 2:
        # this obj has only 2 buildings
        buildingA = buildingList[bids[0]]
        buildingB = buildingList[bids[1]]
        PointA = (buildingA.centroid.x, buildingA.centroid.y)
        PointB = (buildingB.centroid.x, buildingB.centroid.y)
        return PointA, PointB
    else:
        # this obj has at least 3 buildings
        bidA = 0
        bidB = 1
        longestDist = buildingList[bids[0]].centroid.distance(buildingList[bids[1]].centroid)  # NOQA
        for i in range(0, len(bids)):
            for j in range(i+1, len(bids)):
                currentDist = buildingList[bids[i]].centroid.distance(buildingList[bids[j]].centroid)  # NOQA
                if currentDist > longestDist:
                    bidA = i
                    bidB = j
                    longestDist = currentDist

        buildingA = buildingList[bids[bidA]]
        buildingB = buildingList[bids[bidB]]
        PointA = (buildingA.centroid.x, buildingA.centroid.y)
        PointB = (buildingB.centroid.x, buildingB.centroid.y)
        return PointA, PointB


# for a buildings with bid, decide the nearest substation with sid using a given sids  # NOQA
def bdToNearestSub(bid, sids):
    building = buildingList[bid]
    initialSubstation = substationList[sids[0]]
    nearestDistance = building.geom.distance(initialSubstation.geom)
    nearestSid = sids[0]
    for sid in sids:
        testSubstation = substationList[sid]
        testDistance = building.geom.distance(testSubstation.geom)
        if testDistance < nearestDistance:
            nearestDistance = testDistance
            nearestSid = sid
    return nearestSid


# Cut a line using a break point, with the help of Shapely library
def cut(line, point):
    distance = line.project(point)
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i + 1]),
                LineString(coords[i:])]
        if pd > distance:
            return [
                LineString(coords[:i] + [(point.x, point.y)]),
                LineString([(point.x, point.y)] + coords[i:])]


# using a list of points to cut the line segment
def multiCut(line, points):

    cutLines = []
    referenceToPoint = {}  # dictionary to map the refererence to point on this line  # NOQA
    references = []

    for point in points:
        reference = line.project(point)
        referenceToPoint[reference] = point
        references.append(reference)

    tempReferences = []  # remove the duplicate references
    for reference in references:
        if reference not in tempReferences:
            tempReferences.append(reference)

    references = tempReferences
    references.sort()  # sort the references

    untilEnd = False  # it's used to indicate if the last break point just hits the end vertex of the line  # NOQA
    tempLine = line  # Let's start to do multi cut at this stage
    while len(references) != 0:
        reference = references.pop(0)
        breakPoint = referenceToPoint[reference]
        currentReference = tempLine.project(breakPoint)
        if currentReference != 0 and currentReference != tempLine.length:
            cutLine, tempLine = cut(tempLine, breakPoint)
            cutLines.append(cutLine)
        elif currentReference == 0:
            pass
        elif currentReference == tempLine.length:
            cutLines.append(tempLine)
            untilEnd = True

    if untilEnd == False:  # the last break point doesn't hit the end vertex of the original line  # NOQA
        cutLines.append(tempLine)

    return cutLines


def roundNode(coordinate):
    return round(coordinate[0], 3), round(coordinate[1], 3)


def roundPoint(point):
    return roundNode((point.x, point.y))


# for case A, use the cutLines (line lists) to modify the base network, AND return break node list  # NOQA
def modifyGraph(graph, cutLines, extraNodeList):
    startNode = roundNode(list(cutLines[0].coords)[0])
    endNode = roundNode(list(cutLines[-1].coords)[-1])
    if(startNode, endNode) in graph.edges():
        graph.remove_edge(startNode, endNode)
    elif(endNode, startNode) in graph.edges():
        graph.remove_edge(endNode, startNode)
    else:
        print("for substation %d, faliure to remove edge whose rid is %d" % (sid, rid))  # NOQA

    breakNodeList = []
    for line in cutLines:
        coordinate = list(line.coords)
        tempStartNode = roundNode(coordinate[0])
        tempStartGeom = coordinate[0]
        tempEndNode = roundNode(coordinate[-1])
        tempEndGeom = coordinate[-1]
        graph.add_edge(tempStartNode, tempEndNode)
        graph.edge[tempStartNode][tempEndNode]['Coords'] = coordinate
        graph.edge[tempStartNode][tempEndNode]['Length'] = line.length
        graph.node[tempStartNode]['Coords'] = tempStartGeom
        graph.node[tempEndNode]['Coords'] = tempEndGeom
        breakNodeList.append(tempEndNode)
    breakNodeList.pop(-1)  # we don't need the last node

    for extraNode in extraNodeList:  # what if we have at most two extra nodes here?  # NOQA
        breakNodeList.append(extraNode)

    return breakNodeList


def caseA():
    # this is the case where substation.roadId is not in ridToBids.keys()
    # which means the road accessed by subsation will NOT be accessed by other buildings  # NOQA
    edgeGetNewNodes = {}

    # modify the graph by the roads needed
    for rid in ridToBids.keys():
        breakPoints = []  # used to store the break points on this road of given rid  # NOQA
        extraNodeList = []
        for bid in ridToBids[rid]:
            building = buildingList[bid]
            breakPoints.append(building.accessPoint)
            testNode = roundPoint(building.accessPoint)
            if testNode in roadList[rid].edge:
                extraNodeList.append(testNode)
        cutLines = multiCut(roadList[rid].geom, breakPoints)  # we get a list of lines that is cut by breakpoints  # NOQA
        edgeGetNewNodes[roadList[rid].edge] = modifyGraph(roadNetCopy,
                                                          cutLines,
                                                          extraNodeList)

    # now modify the graph finally using accessPoint of substation
    rid = substation.roadId
    subAccessNode = roundPoint(substation.accessPoint)
    road = roadList[rid]
    startNode, endNode = road.edge
    if startNode != subAccessNode and endNode != subAccessNode:
        line1, line2 = cut(road.geom, substation.accessPoint)
        roadNetCopy.remove_edge(startNode, endNode)
        roadNetCopy.add_edge(startNode, subAccessNode)
        roadNetCopy.edge[startNode][subAccessNode]['Coords'] = list(line1.coords)  # NOQA
        roadNetCopy.edge[startNode][subAccessNode]['Length'] = line1.length
        roadNetCopy.add_edge(subAccessNode, endNode)
        roadNetCopy.edge[subAccessNode][endNode]['Coords'] = list(line2.coords)
        roadNetCopy.edge[subAccessNode][endNode]['Length'] = line2.length

        roadNetCopy.node[startNode]['Coords'] = list(line1.coords)[0]
        roadNetCopy.node[subAccessNode]['Coords'] = list(line1.coords)[-1]
        roadNetCopy.node[endNode]['Coords'] = list(line2.coords[-1])

    # Let's create a new graph instance and draw its edges!
    distriGraph = nx.Graph()
    allPaths = []

    # again, need to loop here via the node pair
    for nodePair in edgeGetNewNodes.keys():
        node0 = nodePair[0]
        node1 = nodePair[1]
        nodeX = None  # need to define nodeX, otherwise can't be referenced

        # if there is path ...
        if nx.has_path(roadNetCopy, subAccessNode, node0):
            distance0 = nx.shortest_path_length(roadNetCopy, subAccessNode,
                                                node0, weight='Length')
            distance1 = nx.shortest_path_length(roadNetCopy, subAccessNode,
                                                node1, weight='Length')
            if distance0 < distance1:
                nodeX = node0
            else:
                nodeX = node1
            # get the path from subAccessNode to nodeX
            if subAccessNode != nodeX:
                if nx.has_path(roadNetCopy, subAccessNode, nodeX):  # this check is very important  # NOQA
                    pathToX = nx.dijkstra_path(roadNetCopy, subAccessNode,
                                               nodeX, weight='Length')
                    allPaths.append(pathToX)
        # ok sadly, node0 or node1 cannot connect subAccessNode
        else:
            distance0 = Point(node0).distance(Point(subAccessNode))
            distance1 = Point(node1).distance(Point(subAccessNode))
            if distance0 < distance1:
                nodeX = node0
            else:
                nodeX = node1
            distriGraph.add_edge(nodeX, subAccessNode)
            coords1 = roadNetCopy.node[nodeX]['Coords']
            coords2 = roadNetCopy.node[subAccessNode]['Coords']
            distriGraph.edge[nodeX][subAccessNode]['Coords'] = [coords1, coords2]  # NOQA
            distriGraph.edge[nodeX][subAccessNode]['Length'] = Point(nodeX).distance(Point(subAccessNode))  # NOQA
            distriGraph.edge[nodeX][subAccessNode]['type'] = 'feeder'  # NOQA

        # get the path within the node pair
        for breakNode in edgeGetNewNodes[nodePair]:
                pathBreakToX = nx.dijkstra_path(roadNetCopy, breakNode,
                                                nodeX, weight='Length')
                allPaths.append(pathBreakToX)

    # Yep, we got allPaths, it's time to write them in DistriGraph!
    for path in allPaths:
        for index in range(len(path) - 1):
            startNode = path[index]
            endNode = path[index + 1]
            if (startNode, endNode) not in distriGraph.edges():
                distriGraph.add_edge(startNode, endNode)
                distriGraph.edge[startNode][endNode]['Coords'] = roadNetCopy.edge[startNode][endNode]['Coords']  # NOQA
                distriGraph.edge[startNode][endNode]['Length'] = roadNetCopy.edge[startNode][endNode]['Length']  # NOQA
                distriGraph.edge[startNode][endNode]['type'] = 'feeder'

    for node in distriGraph.nodes():
        distriGraph.node[node]['type'] = 'distribution'
        distriGraph.node[node]['Coords'] = roadNetCopy.node[node]['Coords']

    # not to forget add edges linking buildings to the access nodes
    for rid in ridToBids.keys():
        for bid in ridToBids[rid]:
            building = buildingList[bid]
            startNode = roundPoint(building.centroid)
            endNode = roundPoint(building.accessPoint)
            coords1 = (building.centroid.x, building.centroid.y)
            coords2 = (building.accessPoint.x, building.accessPoint.y)
            distriGraph.add_edge(startNode, endNode)
            distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
            distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
            distriGraph.edge[startNode][endNode]['type'] = 'service'
            distriGraph.node[startNode]['type'] = 'building'
            distriGraph.node[startNode]['toid'] = buildingList[bid].toid
            distriGraph.node[endNode]['type'] = 'buildingAccess'

            distriGraph.node[startNode]['Coords'] = coords1
            distriGraph.node[endNode]['Coords'] = coords2

    # finally add the last edge, which links substation with subAccessNode
    startNode = roundPoint(substation.geom)
    endNode = subAccessNode
    coords1 = (substation.geom.x, substation.geom.y)
    coords2 = (substation.accessPoint.x, substation.accessPoint.y)
    distriGraph.add_edge(startNode, endNode)
    distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
    distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
    distriGraph.edge[startNode][endNode]['type'] = 'feeder'
    distriGraph.node[startNode]['type'] = 'substation'
    distriGraph.node[endNode]['type'] = 'substationAccess'

    distriGraph.node[startNode]['Coords'] = coords1
    distriGraph.node[endNode]['Coords'] = coords2

    # well really finally, check if we need to draw edges directly connecting sub and building  # NOQA
    for bid in substation.directBids:
        print("direct edge: sid is %d, bid is %d" % (substation.id, bid))
        building = buildingList[bid]
        startNode = roundPoint(substation.geom)
        endNode = roundPoint(building.centroid)
        coords1 = (substation.geom.x, substation.geom.y)
        coords2 = (building.centroid.x, building.centroid.y)
        distriGraph.add_edge(startNode, endNode)
        distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
        distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
        distriGraph.edge[startNode][endNode]['type'] = 'service'
        distriGraph.node[startNode]['type'] = 'substation'
        distriGraph.node[endNode]['type'] = 'building'
        distriGraph.node[endNode]['toid'] = buildingList[bid].toid

        distriGraph.node[startNode]['Coords'] = coords1
        distriGraph.node[endNode]['Coords'] = coords2

    # Make sure that all the non-building nodes' toids are 0
    for node in distriGraph.nodes():
        if distriGraph.node[node]['type'] != 'building':
            distriGraph.node[node]['toid'] = '0'

    # Get the Distribution Network Finally!
    return distriGraph


def caseB():
    # this is the case where substation.roadId is in ridToBidss.keys()
    # which means the road accessed by subsation will also be accessed by other buildings  # NOQA

    edgeGetNewNodes = {}

    specBidList = ridToBids.pop(substation.roadId)

    # use the remaining NON-special roads to modify the graph
    for rid in ridToBids.keys():
        breakPoints = []  # used to store the break points on this road of given rid  # NOQA
        extraNodeList = []
        for bid in ridToBids[rid]:
            building = buildingList[bid]
            breakPoints.append(building.accessPoint)
            testNode = roundPoint(building.accessPoint)
            if testNode in roadList[rid].edge:
                extraNodeList.append(testNode)
        cutLines = multiCut(roadList[rid].geom, breakPoints)
        edgeGetNewNodes[roadList[rid].edge] = modifyGraph(roadNetCopy, cutLines, extraNodeList)  # NOQA

    # now use that SPECIAL road to modify the graph
    rid = substation.roadId
    subAccessNode = roundPoint(substation.accessPoint)
    road = roadList[rid]
    breakPoints = []
    extraNodeList = []
    for bid in specBidList:
        breakPoint = buildingList[bid].accessPoint
        breakPoints.append(breakPoint)
        testNode = roundPoint(breakPoint)
        if testNode in road.edge:
            extraNodeList.append(testNode)
    breakPoints.append(substation.accessPoint)  # not to forget the substation accessPoint  # NOQA

    cutLines = multiCut(road.geom, breakPoints)
    startNode = roundNode(list(cutLines[0].coords)[0])
    endNode = roundNode(list(cutLines[-1].coords)[-1])
    roadNetCopy.remove_edge(startNode, endNode)

    specNodeList = []  # different from case A, here we need a list to store the nodes at the SAME edge as subAccessNode  # NOQA

    for line in cutLines:
        coordinate = list(line.coords)
        tempStartNode = roundNode(coordinate[0])
        tempEndNode = roundNode(coordinate[-1])
        roadNetCopy.add_edge(tempStartNode, tempEndNode)
        roadNetCopy.edge[tempStartNode][tempEndNode]['Coords'] = coordinate
        roadNetCopy.edge[tempStartNode][tempEndNode]['Length'] = line.length

        roadNetCopy.node[tempStartNode]['Coords'] = coordinate[0]
        roadNetCopy.node[tempEndNode]['Coords'] = coordinate[-1]

        specNodeList.append(tempEndNode)
    specNodeList.pop(-1)  # we don't need the last node
    for extraNode in extraNodeList:
        specNodeList.append(extraNode)

    # Let's create a new graph instance and draw its edges!
    distriGraph = nx.Graph()
    allPaths = []

    # first let's deal with the SPECIAL node list
    for specNode in specNodeList:
        if subAccessNode != specNode:
            if nx.has_path(roadNetCopy, subAccessNode, specNode):
                pathToSpecNode = nx.dijkstra_path(roadNetCopy, subAccessNode,
                                                  specNode, weight='Length')
                allPaths.append(pathToSpecNode)

    # now deal with the normal nodes, here need to loop via the node pair
    for nodePair in edgeGetNewNodes.keys():
        node0 = nodePair[0]
        node1 = nodePair[1]
        nodeX = None

        if nx.has_path(roadNetCopy, node0, subAccessNode):
            distance0 = nx.shortest_path_length(roadNetCopy, subAccessNode,
                                                node0, weight='Length')
            distance1 = nx.shortest_path_length(roadNetCopy, subAccessNode,
                                                node1, weight='Length')
            if distance0 < distance1:
                nodeX = node0
            else:
                nodeX = node1
            # get the path from subAccessNode to nodeX
            if subAccessNode != nodeX:
                if nx.has_path(roadNetCopy, subAccessNode, nodeX):
                    pathToX = nx.dijkstra_path(roadNetCopy, subAccessNode,
                                               nodeX, weight='Length')
                    allPaths.append(pathToX)

        else:  # ok sadly, node0 or node1 cannot connect subAccessNode
            distance0 = Point(node0).distance(Point(subAccessNode))
            distance1 = Point(node1).distance(Point(subAccessNode))
            if distance0 < distance1:
                nodeX = node0
            else:
                nodeX = node1
                distriGraph.add_edge(nodeX, subAccessNode)
                coords1 = roadNetCopy.node[nodeX]['Coords']
                coords2 = roadNetCopy.node[subAccessNode]['Coords']
                distriGraph.edge[nodeX][subAccessNode]['Coords'] = [coords1, coords2]  # NOQA
                distriGraph.edge[nodeX][subAccessNode]['Length'] = Point(nodeX).distance(Point(subAccessNode))  # NOQA
                distriGraph.edge[nodeX][subAccessNode]['type'] = 'feeder'  # NOQA

        # get the path within the node Pair
        for breakNode in edgeGetNewNodes[nodePair]:
            if nx.has_path(roadNetCopy, breakNode, nodeX):
                pathBreakToX = nx.dijkstra_path(roadNetCopy, breakNode,
                                                nodeX, weight='Length')
                allPaths.append(pathBreakToX)

    # Yep, we got allPaths, it's time to write them in DistriGraph!
    for path in allPaths:
        for index in range(len(path) - 1):
            startNode = path[index]
            endNode = path[index + 1]
            if (startNode, endNode) not in distriGraph.edges():
                distriGraph.add_edge(startNode, endNode)
                distriGraph.edge[startNode][endNode]['Coords'] = roadNetCopy.edge[startNode][endNode]['Coords']  # NOQA
                distriGraph.edge[startNode][endNode]['Length'] = roadNetCopy.edge[startNode][endNode]['Length']  # NOQA
                distriGraph.edge[startNode][endNode]['type'] = 'feeder'

    for node in distriGraph.nodes():
        distriGraph.node[node]['type'] = 'distribution'
        distriGraph.node[node]['Coords'] = roadNetCopy.node[node]['Coords']

    # not to forget add edes linking buildings to the access node
    for rid in ridToBids.keys():  # normal buildings
        for bid in ridToBids[rid]:
            building = buildingList[bid]
            startNode = roundPoint(building.centroid)
            endNode = roundPoint(building.accessPoint)
            coords1 = (building.centroid.x, building.centroid.y)
            coords2 = (building.accessPoint.x, building.accessPoint.y)
            distriGraph.add_edge(startNode, endNode)
            distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
            distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
            distriGraph.edge[startNode][endNode]['type'] = 'service'
            distriGraph.node[startNode]['type'] = 'building'
            if startNode == (415709.827, 566713.464):
                print("bingo")
            distriGraph.node[startNode]['toid'] = buildingList[bid].toid
            distriGraph.node[endNode]['type'] = 'buildingAccess'

            distriGraph.node[startNode]['Coords'] = coords1
            distriGraph.node[endNode]['Coords'] = coords2

    for bid in specBidList:  # special buildings
        building = buildingList[bid]
        startNode = roundPoint(building.centroid)
        endNode = roundPoint(building.accessPoint)
        coords1 = (building.centroid.x, building.centroid.y)
        coords2 = (building.accessPoint.x, building.accessPoint.y)
        distriGraph.add_edge(startNode, endNode)
        distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
        distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
        distriGraph.edge[startNode][endNode]['type'] = 'service'
        distriGraph.node[startNode]['type'] = 'building'
        distriGraph.node[startNode]['toid'] = buildingList[bid].toid
        distriGraph.node[endNode]['type'] = 'buildingAccess'

        distriGraph.node[startNode]['Coords'] = coords1
        distriGraph.node[endNode]['Coords'] = coords2

    # finally add the last edge, which links substation with subAccessNode
    startNode = roundPoint(substation.geom)
    endNode = subAccessNode
    coords1 = (substation.geom.x, substation.geom.y)
    coords2 = (substation.accessPoint.x, substation.accessPoint.y)
    distriGraph.add_edge(startNode, endNode)
    distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
    distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
    distriGraph.edge[startNode][endNode]['type'] = 'feeder'
    distriGraph.node[startNode]['type'] = 'substation'
    distriGraph.node[endNode]['type'] = 'substationAccess'

    distriGraph.node[startNode]['Coords'] = coords1
    distriGraph.node[endNode]['Coords'] = coords2

    # well really finally, check if we need to draw edges directly connecting sub and building  # NOQA
    for bid in substation.directBids:
        print("direct edge: sid is %d, bid is %d" % (substation.id, bid))
        building = buildingList[bid]
        startNode = roundPoint(substation.geom)
        endNode = roundPoint(building.centroid)
        coords1 = (substation.geom.x, substation.geom.y)
        coords2 = (building.centroid.x, building.centroid.y)
        distriGraph.add_edge(startNode, endNode)
        distriGraph.edge[startNode][endNode]['Coords'] = [coords1, coords2]
        distriGraph.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA
        distriGraph.edge[startNode][endNode]['type'] = 'service'
        distriGraph.node[startNode]['type'] = 'substation'

        distriGraph.node[endNode]['type'] = 'building'
        distriGraph.node[endNode]['toid'] = buildingList[bid].toid

        distriGraph.node[startNode]['Coords'] = coords1
        distriGraph.node[endNode]['Coords'] = coords2

    # Make sure that all the non-building nodes' toids are 0
    for node in distriGraph.nodes():
        if distriGraph.node[node]['type'] != 'building':
            distriGraph.node[node]['toid'] = '0'

    # Get the Distribution Network Finally!
    return distriGraph

def writeGraph(graph, loopid):
    # Write the edges first
    sourceDriver = 'ESRI Shapefile'
    sourceCrs = {'y_0': -100000, 'units': 'm', 'lat_0': 49,
                 'lon_0': -2, 'proj': 'tmerc', 'k': 0.9996012717,
                 'no_defs': True, 'x_0': 400000, 'datum': 'OSGB36'}

    result_folder = "result//Edges//"

    # write the network edges
    sourceSchema = {'properties': {'Length': 'float:19.11', 'type': 'str', 'netID': 'int'},  # NOQA
                    'geometry': 'LineString'}
    fileName = result_folder + 'Edges' + str(loopid) + '.shp'
    with fiona.open(fileName,
                    'w',
                    driver=sourceDriver,
                    crs=sourceCrs,
                    schema=sourceSchema) as source:
        for edge in graph.edges():
            startNode = edge[0]
            endNode = edge[1]
            record = {}
            thisEdge = graph.edge[startNode][endNode]
            record['geometry'] = {'coordinates': thisEdge['Coords'], 'type': 'LineString'}  # NOQA
            record['properties'] = {'Length': thisEdge['Length'],
                                    'type': thisEdge['type'],
                                    'netID': loopid}  # NOQA
            source.write(record)

    # write the nodes then
    sourceDriver = 'ESRI Shapefile'
    sourceCrs = {'y_0': -100000, 'units': 'm', 'lat_0': 49,
                 'lon_0': -2, 'proj': 'tmerc', 'k': 0.9996012717,
                 'no_defs': True, 'x_0': 400000, 'datum': 'OSGB36'}

    result_folder = "result//Nodes//"

    # write the network edges
    sourceSchema = {'properties': {'type': 'str', 'toid': 'str', 'netID': 'int'},  # NOQA
                    'geometry': 'Point'}
    fileName = result_folder + 'Nodes' + str(loopid) + '.shp'
    with fiona.open(fileName,
                    'w',
                    driver=sourceDriver,
                    crs=sourceCrs,
                    schema=sourceSchema) as source:
        for node in graph.nodes():
            thisNode = graph.node[node]
            record = {}
            record['geometry'] = {'coordinates': thisNode['Coords'], 'type': 'Point'}  # NOQA
            record['properties'] = {'type': thisNode['type'], 'toid': thisNode['toid'], 'netID': loopid}  # NOQA
            source.write(record)

"""
===============================================================================
                     From here is the main programme
===============================================================================
"""

"""
================  First let's read all necessary data into memory =============
"""


"""
Read the road network data into networkx instance in python FROM POSTGIS,
this will be the base network in the whole algorithm
"""

roadNet = pg_read.readNet(dbname, password, user)

print("=====================================================")
print("Reading the input data for this algorithm............")
print("=====================================================")

"""
Read each road into RoadList
"""
roadList = pg_read.readRoad(dbname, password, user)
print("%d roads have been read\n" % len(roadList))

"""
Read each building into buildingList
"""
buildingList = pg_read.readBuilding(dbname, password, user)
print("%d buildings have been read\n" % len(buildingList))

"""
Read each substation into substationList
"""
substationList = pg_read.readSubstation(dbname, password, user)
print("%d substations have been read\n" % len(substationList))

"""
Read each terrace into terraceList
"""
terraceList = pg_read.readTerrace(dbname, password, user)
print("%d terraces have been read\n" % len(terraceList))

"""
Read each union into unionList
"""
unionList = pg_read.readUnion(dbname, password, user)
print("%d unions have been read\n" % len(unionList))

print("=====================================================")
print("Data Reading is completed............................")
print("=====================================================")
print("\n")


"""
========= Now we will relate different class objects with each other ==========
"""

conn = psycopg2.connect("dbname = %s password = %s user = %s" % (dbname, password, user))  # NOQA
cur = conn.cursor()

"""
Map each road to the network edges
"""
for road in roadList:
    startNodeX = round(list(road.geom.coords)[0][0], 3)
    startNodeY = round(list(road.geom.coords)[0][1], 3)
    endNodeX = round(list(road.geom.coords)[-1][0], 3)
    endNodeY = round(list(road.geom.coords)[-1][1], 3)

    # road.edge: a tuple, it's the network edge representation of this road segment   # NOQA
    road.edge = ((startNodeX, startNodeY), (endNodeX, endNodeY))
print("Mapping each road to the network edges done \n")

"""
Assign a road to each substation, ie,
calculate access point and get id of the nearest road
"""
cur.execute("select distinct on (sid) sid, rid \
            from substations as s, roads as r \
            where st_dwithin(s.geom, r.geom, 100000) \
            order by sid, st_distance(s.geom, r.geom)")
results = cur.fetchall()
if len(results) != len(substationList):
    print("critical error in assigning a road to each substation")
    sys.exit()
for sid in range(len(substationList)):
    rid = results[sid][1]
    road = roadList[rid]
    substationList[sid].roadId = rid
    """
    I add the modified code here!
    """
    reference = road.geom.project(substationList[sid].geom)
    substationList[sid].accessPoint = road.geom.interpolate(reference)

print("Assign a road to each substation done \n")

"""
Assignment among buildings, unions and roads
"""
# a union should contain buildings
cur.execute("select bid, uid \
            from buildings as b, unions as u \
            where st_intersects(b.geom, u.geom) \
            order by bid, st_distance(b.geom, u.geom)")  # NOQA
results = cur.fetchall()
if len(results) != len(buildingList):
    print("critical error in assigning a union to each building")
    sys.exit()

for line in results:
    bid = line[0]
    uid = line[1]
    # Note these two lines make sure buildings and unions are linked together
    buildingList[bid].uid = uid
    unionList[uid].buildingIdList.append(bid)

# DECIDE the access road for each union
# In here, we only consider the nearest road segment for each union
# We don't consider angle, because the union cannot necessarily be
# a line-structure! Therefore, assigning the accessRoad and accessPoint
# to each union is still straightforward
for union in unionList:
    uid = union.id
    # for this union, the nearest road segment to it in 4000 meters
    cur.execute("select uid, rid, st_distance(u.cp, r.geom) \
                from unions as u, roads as r \
                where st_dwithin(u.cp, r.geom, 4000) \
                and uid = %d \
                order by uid, st_distance(u.cp, r.geom) \
                limit 1" % uid)
    results = cur.fetchall()
    result = results[0]
    rid = result[1]
    road = roadList[rid]
    unionList[uid].roadId = rid
    reference = road.geom.project(unionList[uid].centroid)
    unionList[uid].accessPoint = road.geom.interpolate(reference)
print("Assign union to building, and road to union done \n")


"""
Assign a substation to each union using the base network!
"""
ridToBreakPoints = {}
for union in unionList:
    rid = union.roadId
    if rid not in ridToBreakPoints.keys():
        ridToBreakPoints[rid] = []
    if union.accessPoint not in ridToBreakPoints[rid]:
        ridToBreakPoints[rid].append(union.accessPoint)

for substation in substationList:
    rid = substation.roadId
    if rid not in ridToBreakPoints.keys():
        ridToBreakPoints[rid] = []
    if substation.accessPoint not in ridToBreakPoints[rid]:
        ridToBreakPoints[rid].append(substation.accessPoint)

# Now the really tricky part, use all the BreakPoints to modify necessary edges
edgeGetNewNodes = {}

baseNet = copy.deepcopy(roadNet)

print("before modifying, the baseNet has %d edges" % baseNet.number_of_edges())

for rid in ridToBreakPoints.keys():
    extraNodeList = []
    breakPoints = ridToBreakPoints[rid]
    for breakPoint in breakPoints:
        testNode = roundPoint(breakPoint)
        if testNode in roadList[rid].edge:
            extraNodeList.append(testNode)
    cutLines = multiCut(roadList[rid].geom, breakPoints)
    edgeGetNewNodes[roadList[rid].edge] = modifyGraph(baseNet, cutLines, extraNodeList)  # NOQA
print("after modyfing, the baseNet has %d edges" % baseNet.number_of_edges())

# Now derive edges from terraces/substations to the baseNetwork as well
for union in unionList:
    startNode = roundPoint(union.centroid)
    endNode = roundPoint(union.accessPoint)
    coords1 = (union.centroid.x, union.centroid.y)
    coords2 = (union.accessPoint.x, union.accessPoint.y)
    baseNet.add_edge(startNode, endNode)
    baseNet.edge[startNode][endNode]['Coords'] = [coords1, coords2]
    baseNet.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA

for substation in substationList:
    startNode = roundPoint(substation.geom)
    endNode = roundPoint(substation.accessPoint)
    coords1 = (substation.geom.x, substation.geom.y)
    coords2 = (substation.accessPoint.x, substation.accessPoint.y)
    baseNet.add_edge(startNode, endNode)
    baseNet.edge[startNode][endNode]['Coords'] = [coords1, coords2]
    baseNet.edge[startNode][endNode]['Length'] = LineString([coords1, coords2]).length  # NOQA

print("after deriving addtional edges, the baseNet has %d edges" % baseNet.number_of_edges())  # NOQA

"""
Visualize the basenetwork!
"""

"""
sourceDriver = 'ESRI Shapefile'
sourceCrs = {'y_0': -100000, 'units': 'm', 'lat_0': 49,
             'lon_0': -2, 'proj': 'tmerc', 'k': 0.9996012717,
             'no_defs': True, 'x_0': 400000, 'datum': 'OSGB36'}

result_folder = "result//"

# write the network edges
sourceSchema = {'properties': {'Length': 'float:19.11'},
                'geometry': 'LineString'}
fileName = result_folder + 'Base Network' + '.shp'
with fiona.open(fileName,
                'w',
                driver=sourceDriver,
                crs=sourceCrs,
                schema=sourceSchema) as source:
    for edge in baseNet.edges():
        startNode = edge[0]
        endNode = edge[1]
        record = {}
        record['geometry'] = {'coordinates': baseNet.edge[startNode][endNode]['Coords'], 'type': 'LineString'}  # NOQA
        record['properties'] = {'Length': baseNet.edge[startNode][endNode]['Length']}  # NOQA
        source.write(record)
sys.exit()
"""

# Let's now calculate the nearest substation (path distance) to each union  # NOQA
uidSids = Triangulation.uidToPotentialSids(dbname, password, user)

if len(uidSids.keys()) != len(unionList):
    print("it seems not all the unions are covered by triangles, error!")
    sys.exit()

for uid in range(len(unionList)):
    sids = uidSids[uid]
    print("union-sub", uid)

    """
    Basically, this approach is very fine.

    However, if we can't find a path from the union to the subsation,
    what should we do?

    Case A: at least one substation has direct path to the union, choose that
    Case B: Unfortunately, all the substations is not connected to the union
            via the path, so let's ... choose the nearest sub of Euc Distance
            This might not be an optimal approach, but let's try it first
    """
    union = unionList[uid]
    startNode = roundPoint(union.centroid)
    chosenSid = -1
    minDist = 100000000000
    for sid in sids:
        substation = substationList[sid]
        endNode = roundPoint(substation.geom)
        if nx.has_path(baseNet, startNode, endNode):
            currentDist = nx.shortest_path_length(baseNet, startNode, endNode, weight='Length')  # NOQA
        else:  # in some cases, no network path avaible
            currentDist = 10000000000 + Point(startNode).distance(Point(endNode))  # NOQA
        if currentDist < minDist:
            chosenSid = sid
            minDist = currentDist
    unionList[uid].substationId = chosenSid

print("Assign a substation to each union done \n")


"""
Double direction assignment between buildings and substations!
"""
for union in unionList:
    bids = union.buildingIdList
    sid = union.substationId
    for bid in bids:
        buildingList[bid].substationId = sid
        substationList[sid].buildingIdList.append(bid)
print("Double direction assignment beweetn buildings to substations done \n")

"""
Check if a union(has no more than 1 building) is too close to its servicing substation (< 50m)  # NOQA
"""
cur.execute("select distinct on (u.uid) u.uid, s.sid \
            from unions as u, substations as s \
            where st_dwithin(u.geom, s.geom, 50) \
            order by u.uid, st_distance(u.geom, s.geom)")
results = cur.fetchall()
for result in results:
    uid = result[0]
    sid = result[1]
    union = unionList[uid]
    if union.substationId == sid and len(union.buildingIdList) == 1:
        bid = union.buildingIdList[0]
        pt1 = union.centroid
        pt2 = substationList[sid].geom
        testLine = LineString([(pt1.x, pt1.y), (pt2.x, pt2.y)])
        cur.execute("select * from buildings \
                    where st_intersects(geom, st_geomfromtext('%s', 27700)) \
                    and bid != '%d'" % (testLine.wkt, bid))
        results1 = cur.fetchall()
        cur.execute("select * from roads \
                    where st_intersects(geom, st_geomfromtext('%s', 27700))" % testLine.wkt)  # NOQA
        results2 = cur.fetchall()
        if not results1 and not results2:
            buildingList[bid].connectType = "direct"
            substationList[sid].directBids.append(bid)
print("checking if a union too close to substation done!")
print("\n")

"""
Assign a road to each building but using terrace instead
"""
# first assign each building a terrace, easy
cur.execute("select distinct on(bid) bid, tid  \
            from buildings as b, terraces as t \
            where st_intersects(b.geom, t.geom) \
            and st_area(st_intersection(b.geom, t.geom)) > 0.9 * st_area(b.geom) \
            order by bid")
results = cur.fetchall()
if len(results) != len(buildingList):
    print("error happens when assiging each building a terrace")
    sys.exit()

for result in results:
    bid = result[0]
    tid = result[1]
    buildingList[bid].terraceId = tid
    terraceList[tid].buildingIdList.append(bid)

print("Potential Error Check - step1 done")

# then assign each terrace a road
# this is the hardest part in the entire algorithm
# the ITN road network is bad somewhere in the city
# that's why assigning an approriate road to a terrace is nightmare

# for each terrace, find nearest 3 road segments
for terrace in terraceList:
    tid = terrace.id
    cur.execute("select tid, rid, st_distance(t.cp, r.geom) \
                from terraces as t, roads as r \
                where st_dwithin(t.cp, r.geom, 2000) \
                and tid = %d \
                order by tid, st_distance(t.cp, r.geom) \
                limit 3" % tid)
    results = cur.fetchall()
    if len(results) < 3:
        # which means less than 3 nearest roads to the terraces
        # or the terrace itself consists of 1 building
        result = results[0]
        # only use the nearest road
        rid = result[1]
        road = roadList[rid]
        terraceList[tid].roadId = rid
        reference = road.geom.project(terraceList[tid].centroid)
        terraceList[tid].accessPoint = road.geom.interpolate(reference)
    else:
        # now things are really complicated, the terrace has 3 nearest roads        
        roads = []
        for result in results:
            rid = result[1]
            roads.append(roadList[rid])

        if len(terrace.buildingIdList) == 1:
            # only one building in the terrace
            rid, accessPoint = simpleObjPickRoad(terrace, roads)
        else:
            # several buildings in the terrace
            rid, accessPoint = pickFittestRoad(terrace, roads)

        terraceList[tid].roadId = rid
        terraceList[tid].accessPoint = accessPoint

print("Potential Error Check, step 2 - done")

# last step is to assign each building a road!!!!
for terrace in terraceList:
    rid = terrace.roadId
    road = roadList[rid]
    bids = terrace.buildingIdList
    for bid in bids:
        if terrace.projectType == 'normal':
            buildingList[bid].roadId = rid
            reference = road.geom.project(buildingList[bid].centroid)
            buildingList[bid].accessPoint = road.geom.interpolate(reference)
        else:
            buildingList[bid].roadId = rid
            buildingList[bid].accessPoint = terrace.accessPoint


print("Potential Error Check, step 3 - done")

print("congratualtions, you have assigned a road to each building!")

cur.close()
conn.close()


"""
===============================================================================
     The main loop. For each substation calculate its distribution network
===============================================================================
"""

# Create the folders
if not os.path.exists("result"):
    os.mkdir("result")

if not os.path.exists("result//Edges"):
    os.mkdir("result//Edges")

if not os.path.exists("result//Nodes"):
    os.mkdir("result//Nodes")

networks = []

print("=====================================================")
print("Generating network for each substation...............")
print("=====================================================")
print("\n")

for substation in substationList:

    print("generating network No.%d substation......" % substation.id)

    numOfEndConsumers = len(substation.buildingIdList)

    # this substation unluckily serves no buildings
    # so it has no own distribution network
    if numOfEndConsumers == 0:
        print("This substation serves no buildings, and thus has no distribution network!\n")  # NOQA
    else:
        ridToBids = {}  # use a dictionary to group buildings by their accessing road  # NOQA
        # ridToBids will be something like this: {0:[0,1,...], 1:[2,4...]}

        for bid in substation.buildingIdList:
            # this means the building will access the road, breakpoint is useful  # NOQA
            if buildingList[bid].connectType == "indirect":
                rid = buildingList[bid].roadId
                if rid not in ridToBids.keys():
                    ridToBids[rid] = []
                ridToBids[rid].append(bid)

        roadNetCopy = copy.deepcopy(roadNet)
        substation.ridToBids = ridToBids

        # Now we need to divide the situation into 2 cases: see below

        # 1st case: the road accessed by substation will NOT be accessed by buildings  # NOQA
        if substation.roadId not in ridToBids.keys():
            print("substation no.%d goes to caseA" % substation.id)
            distriGraph = caseA()
            finalCheckValidity(distriGraph)
            # before write the graph let's generate feeders
            # identifyFeeder(distriGraph)
            networks.append(distriGraph)
            writeGraph(distriGraph, substation.id)

        # 2nd case: The road accessed by substation will ALSO be accessed by buildings  # NOQA
        else:
            print("substation no.%d goes to caseB" % substation.id)
            distriGraph = caseB()
            finalCheckValidity(distriGraph)
            # before write the graph let's generate feeders
            # identifyFeeder(distriGraph)
            networks.append(distriGraph)
            writeGraph(distriGraph, substation.id)

        print("distribution network servicing %d buildings is generated" % numOfEndConsumers)  # NOQA
        print("\n")

print("=====================================================")
print("Network generation algorithm is finished.............")
print("=====================================================")
