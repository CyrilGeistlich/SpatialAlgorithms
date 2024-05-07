from classes import Point, PointGroup, Segment, Polygon, Bbox
import json

def process_json_file(json_files_data):

    all_data = {}
    for json_file, return_name in json_files_data:
        if return_name not in all_data:
            all_data[return_name] = []

        with open(json_file) as f:
            data = json.load(f)

            for feature in data['features']:
                if feature['geometry']['type'] == 'Point':
                    name = feature['properties']['NAME']
                    coordinates = feature['geometry']['coordinates']
                    x, y = coordinates[:2] 
                    point = Point(name, x, y)
                    all_data[return_name].append(point)
                elif feature['geometry']['type'] == 'Polygon':
                    name = feature['properties']['name']
                    points = feature['geometry']['coordinates'][0]
                    x_coords, y_coords = zip(*points)
                    polygon = Polygon(name, list(zip(x_coords, y_coords)), 0, 1)
                    all_data[return_name].append(polygon)

    return all_data