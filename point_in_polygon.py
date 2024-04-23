# Spatial Algorithm here... 

## Libraries ##

import pandas as pd # Used to import csv
import geopandas as gpd # Used to import gpkg



## Import Data ## 

mun = gpd.read_file('data/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg')

# The SwissNames3D Gazetteer comes with 3 seperate csv files, all containing points.
# The names_lin and names_pol are already reduced to points in the downloaded data set.
# But we devide to ignore since it does not make sense when comparing to other points. 

names_pt = pd.read('data/swissnames3d_2023_2056/swissNAMES3D_PKT.csv')
#names_lin = pd.read('data/swissnames3d_2023_2056/swissNAMES3D_LIN.csv')
#names_pol = pd.read('data/swissnames3d_2023_2056/swissNAMES3D_PLY.csv')