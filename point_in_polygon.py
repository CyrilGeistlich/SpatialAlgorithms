# Spatial Algorithm here... 

## Libraries ##

from classes import *
import pandas as pd # Used to import csv
import geopandas as gpd # Used to import gpkg
import rtree # Used for spatial index
import random
import matplotlib.pyplot as plt

## Import Data ## 

#mun = gpd.read_file('data/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg')

# The SwissNames3D Gazetteer comes with 3 seperate csv files, all containing points.
# The names_lin and names_pol are already reduced to points in the downloaded data set.
# But we devide to ignore since it does not make sense when comparing to other points. 

#names_pt = pd.read_csv('data/swissnames3d_2023_2056/swissNAMES3D_PKT.csv')
#names_lin = pd.read('data/swissnames3d_2023_2056/swissNAMES3D_LIN.csv')
#names_pol = pd.read('data/swissnames3d_2023_2056/swissNAMES3D_PLY.csv')

    
# create an R-tree index for polygons
def create_spatial_index(polygons):
    index = rtree.index.Index() # initialize an R-tree index
    for idx, polygon in enumerate(polygons):    # the index of the polygon in the list, and the actual polygon object       
        index.insert(idx, (polygon.bbox.ll.x,polygon.bbox.ll.y,polygon.bbox.ur.x,polygon.bbox.ur.y)) # get the bounding box and add it to the spatial index
    return index

# which polygon does each point belong to? some points may not belong to any polygon
# to set polygon_id in the class Point
# if activate those lines, can return a list of polygon_id (some of them could be None), keep the same index with the input list "points"
def assign_polygons_to_points(points, polygons, index):
    #results = []
    for point in points:
        candidate_ids = list(index.intersection((point.x, point.y, point.x, point.y))) # get the index of potential polygons that might contain the point
        found = False
        for idx in candidate_ids:
            if polygons[idx].containsPoint(point): # should check our input polygons have this method, otherwise use shapely and def contains(self, point): return self.geometry.contains(point.geometry)
                point.set_polygon_id(polygons[idx])    # if using this line, this method should be add to the class Point as def set_polygon_id(self, polygon_id):self.polygon_id = polygon_id
                #results.append(idx)
                #found = True
                break
        #if not found: 
            #results.append(None)
    #return results 


## POINT IN POLYGON

# TEST DATA
# =========
# a small group of points, spatially ordered

sample1 = [[0,10], [5,0], [10,10], [15,0], [20,10], [25, 0],
             [30, 20], [35, 15], [45, 0], [50, 50], [45, 40], 
             [40, 50], [30, 45], [25, 40], [20, 30], [15, 50],
             [10,35], [5, 50],[5,50], [0, 10]]

samplePolygon = Polygon(sample1, xcol=0, ycol=1)

index=create_spatial_index(samplePolygon)
assign_polygons_to_points(points, samplePolygon,index)


