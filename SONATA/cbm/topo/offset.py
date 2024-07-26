# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:37:44 2016

@author: TPflumm
"""
# Third party modules
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as shp

# Local modules
from .utils import (P2Pdistance, Polygon_orientation,
                    calc_DCT_angles, isclose, unique_rows,)


# Function to check if two line segments (p1, q1) and (p2, q2) intersect
def do_intersect(p1, q1, p2, q2):
    def orientation(p, q, r):
        val = (float(q[1] - p[1]) * (r[0] - q[0])) - (float(q[0] - p[0]) * (r[1] - q[1]))
        if val > 0:
            return 1  # Clockwise
        elif val < 0:
            return 2  # Counterclockwise
        else:
            return 0  # Collinear

    def on_segment(p, q, r):
        if min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1]):
            return True
        return False

    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if o1 != o2 and o3 != o4:
        return True

    # Special cases
    if o1 == 0 and on_segment(p1, q1, p2):
        return True
    if o2 == 0 and on_segment(p1, q1, q2):
        return True
    if o3 == 0 and on_segment(p2, q2, p1):
        return True
    if o4 == 0 and on_segment(p2, q2, q1):
        return True

    return False

# Function to check if a shape intersects itself
def shape_intersects_itself(coords):
    num_points = len(coords)
    for i in range(num_points - 1):  # Loop to (num_points - 1) to avoid out-of-bounds access
        for j in range(i + 2, num_points - (1 if i == 0 else 0)):  # Adjusting the range to avoid out-of-bounds
            if do_intersect(coords[i], coords[i + 1], coords[j], coords[(j + 1) % num_points]):
                return [True, coords[i:j]]
    return [False, coords]

def shp_parallel_offset(arrPts, dist, join_style=1, side="right", res=16):
    # OFFSET ALGORITHM
    # join_style = 1#( 1:round,2:mitre,3:bevels)
    closed = None

    # ==============SHAPELY-OFFSET ALGORITHM====================================
    if P2Pdistance(arrPts[0], arrPts[-1]) <= 1e-6:
        closed = True
        try:
            afpoly = shp.Polygon(arrPts)
            noffafpoly = afpoly.buffer(-dist)  # Inward offset
            data = np.array(noffafpoly.exterior.xy).T
        except:
            intersects = True
            while intersects:
                [intersects, arrPts] = shape_intersects_itself(arrPts)
            afpoly = shp.Polygon(arrPts)
            noffafpoly = afpoly.buffer(-dist)  # Inward offset
            data = np.array(noffafpoly.exterior.xy).T

    else:
        closed = False
        line = shp.LineString(arrPts)
        offset = line.parallel_offset(dist, side, res, join_style)

        if isinstance(offset, shp.MultiLineString):
            parts = hasattr(offset, "geoms") and offset or [offset]
            for part in parts:
                if part.is_closed:
                    x, y = part.xy
                    data = np.vstack((x, y))
                    data = data.T
                else:
                    print("ERROR: \t A multilinestring has been created by shp_parallel_offset that is not closed!")

        elif isinstance(offset, shp.LineString):
            data = np.array(offset.coords)

        else:
            data = np.array(offset.coords)

    # ==============CHECK ORIENTATION if closed=================================
    # Check Orientation and reverse if neccessary
    # TODO: Be careful not to reverse Linestring!
    
    if closed == True:
        Orientation = Polygon_orientation(data)
        if Orientation == True:
            data = np.flipud(data)

    # ==============Interpolate large linear spaces=============================
    seg_P2Plength = []
    cumm_length = 0
    Resolution = 100

    for j in range(0, len(data) - 1):
        seg_P2Plength.append(P2Pdistance(data[j], data[j + 1]))
        cumm_length += P2Pdistance(data[j], data[j + 1])

    # Check if Refinement is necessary:
    if len(seg_P2Plength) > 0 and max(seg_P2Plength) > cumm_length / Resolution:
        Refinement = True
    else:
        Refinement = False

    while Refinement == True:
        temp_data = []
        for i in range(0, len(data) - 1):
            if P2Pdistance(data[i], data[i + 1]) > (cumm_length / Resolution):
                p0 = data[i]
                p1 = data[i + 1]
                v1 = p1 - p0
                p05 = p0 + v1 / 2
                temp_data.append(p0)
                temp_data.append(p05)
            else:
                temp_data.append(data[i])

        temp_data.append(data[-1])
        data = np.vstack(temp_data)

        # Check if further Refinement is necessary
        seg_P2Plength = []
        cumm_length = 0
        for j in range(0, len(data) - 1):
            seg_P2Plength.append(P2Pdistance(data[j], data[j + 1]))
            cumm_length += P2Pdistance(data[j], data[j + 1])

        if max(seg_P2Plength) > cumm_length / Resolution:
            Refinement = True
        else:
            Refinement = False
        # if closed:
        #     np.vstack((data, data[0]))

    return data


# ==============================================================================
if __name__ == "__main__":
    exec(compile(open("SONATA.py").read(), "SONATA.py", "exec"))
