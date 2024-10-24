#!/usr/bin/env python3

# Python script to draw the NHC cone of uncertainty using the forecast error data from the NHC
# Bret Jones & Matthew V. Bilskie, PhD
# Help with MS Copilot
# Oct. 24, 2024

import sys
import shapely
from shapely.geometry import Point, Polygon, LineString, shape
import numpy as np
import json
import geopandas as gpd
import argparse

def parse_coordinates(coord, direction):
    value = float(coord[:-1])
    if direction in ['S', 'W']:
        value = -value
    return round(value * 0.1, 3)

# Create a function to interpolate error for a given forecast period
def interpolate_error(forecast_periods, errors, known_forecast_period):
    # Use numpy's interpolation function
    interpolated_error = np.interp(known_forecast_period, forecast_periods, errors)
    return interpolated_error

def calculate_radius(time, error_list):
    if time > len(error_list):
        return error_list[-1] + 30 * (time - len(error_list))
    return error_list[time]

def create_geojson(union, output_cone_geojson):
    geojson_data = shapely.to_geojson(union)
    geojson_dict = json.loads(geojson_data)
    with open(output_cone_geojson, "w") as f:
        json.dump(geojson_dict, f, indent=4)

def create_shapefile(union, points, output_cone_shapefile, output_points_shapefile):
    geojson_data = shapely.to_geojson(union)
    geojson_dict = json.loads(geojson_data)
    if geojson_dict['type'] == 'Polygon':
        gdf = gpd.GeoDataFrame([{'geometry': shape(geojson_dict)}], crs="EPSG:4326")
        gdf.to_file(output_cone_shapefile)
    else:
        raise ValueError("Unsupported GeoJSON type")

    # Create a GeoDataFrame for the points
    point_geometries = [Point(y, x) for x, y, _ in points]  # Swapped x and y
    points_gdf = gpd.GeoDataFrame(geometry=point_geometries, crs="EPSG:4326")
    #print("Points GeoDataFrame:\n", points_gdf)  # Debug print
    points_gdf.to_file(output_points_shapefile)
    #print(f"{output_points_shapefile} created successfully")  # Debug print

def create_line_shapefile(points, output_line_shapefile):
    # Create a LineString from the points
    line = LineString([(y, x) for x, y, _ in points])  # Swapped x and y
    line_gdf = gpd.GeoDataFrame(geometry=[line], crs="EPSG:4326")
    line_gdf.to_file(output_line_shapefile)
    #print(f"{output_line_shapefile} created successfully")  # Debug print

def main(input_file, output_cone_geojson, output_cone_shapefile, output_points_shapefile, output_line_shapefile):
    with open(input_file, "r") as f:
        lines = f.read().strip().split("\n")

    points = []
    forecast_period = [0, 12, 24, 36, 48, 60, 72, 96, 120] # 2024 Atlantic Basin Forecast Error https://www.nhc.noaa.gov/aboutcone.shtml
    error = [0, 26, 41, 55, 70, 88, 102, 151, 220]
    error_2D = [[forecast_period[i], value] for i, value in enumerate(error)]

    for line in lines:
        parts = line.split(", ")
        if len(parts) >= 7:
            forecast_hour = int(parts[5].strip())
            interpolated_error = interpolate_error(forecast_period, error, forecast_hour)
            x = parse_coordinates(parts[6].strip(), parts[6][-1])
            y = parse_coordinates(parts[7].strip(), parts[7][-1])
            r = interpolated_error
            new_point = [x, y, r]
            if not points or points[-1] != new_point:
                print(f"The interpolated error for a forecast period of {forecast_hour} is {interpolated_error}.")
                points.append(new_point)

    circles = []
    rings = []
    RESOLUTION = 36

    for i in range(1, len(points) - 1):
        item = points[i]
        ring = [
            Point(
                item[1] + np.sin(2 * np.pi * j / RESOLUTION) * item[2] / 60,
                item[0] + np.cos(2 * np.pi * j / RESOLUTION) * item[2] / 60
            )
            for j in range(RESOLUTION)
        ]
        if i > 1:
            combined_ring = rings[-1] + ring
            hull = Polygon(combined_ring).convex_hull
            circles.append(hull)
        rings.append(ring)

    union = shapely.union_all(circles)
    
    create_geojson(union, output_cone_geojson)
    create_shapefile(union, points, output_cone_shapefile, output_points_shapefile)
    create_line_shapefile(points, output_line_shapefile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some coordinates.")
    parser.add_argument("--input", required=True, help="Path to the input file")
    parser.add_argument("--outputPoints", required=True, help="Path to the output points file (Shapefile)")
    parser.add_argument("--outputConeGeoJSON", required=True, help="Path to the output cone file (GeoJSON)")
    parser.add_argument("--outputConeShapefile", required=True, help="Path to the output cone file (Shapefile)")
    parser.add_argument("--outputLine", required=True, help="Path to the output line file (Shapefile)")

    args = parser.parse_args()
    main(args.input, args.outputConeGeoJSON, args.outputConeShapefile, args.outputPoints, args.outputLine)
