import classes
from classes import *
from import_json import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import geopandas as gpd
from shapely.geometry import Polygon as ShapelyPolygon
import json
import folium

## EXECTUE FUNCTIONS ##

# 1. READ DATA (MAREK)

process_json_file()

# (1.5 SPATIAL INDEX???) (TAO)

# 2. POLYGON CLIPPING (SEBASTIAN)


# 3. POINT IN POLYGON (TAO)

# 4. COUNT, ANALYSIS PREPROCESSING (CYRIL)

# 5. ANALYSIS AND VISUALIZATION (JON)

interactive_map = classes.Polygon.viz_interactive()