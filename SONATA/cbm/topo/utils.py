## 2016 Tobias Pflumm (tobias.pflumm@tum.de)
## Note that this is adapted from pythonocc-Utils.
""" This Module describes all helpful utility functionalities of the code"""

# Core Library modules
# Basic Libraries:
import math

# Third party modules
import numpy as np
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeWire,
                                     BRepBuilderAPI_Transform,)
# Python OCC Libraries
from OCC.Core.gp import gp_Pnt, gp_Pnt2d, gp_Trsf, gp_Vec2d
from OCC.Core.TColgp import (TColgp_Array1OfPnt, TColgp_Array1OfPnt2d,
                             TColgp_HArray1OfPnt, TColgp_HArray1OfPnt2d,)
from OCC.Core.TopoDS import topods
from scipy.spatial import cKDTree

# Own Modules:


#######################UTILITY FUNCTIONS######################################
def P2Pdistance(P1, P2):
    """Calculates the distance between two 2D Points"""
    return math.sqrt((P1[0] - P2[0]) ** 2 + (P1[1] - P2[1]) ** 2)

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a, idx = np.unique(a.view([("", a.dtype)] * a.shape[1]), return_index=True)
    unique_a = unique_a[np.argsort(idx)]
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


def isclose(a, b, rel_tol=2e-09, abs_tol=2e-09):
    """
    rel_tol: is the relative tolerance -- it is the amount of error allowed, relative to the larger absolute value of a or b. For
    example, to set a tolerance of 5%, pass tol=0.05. The default tolerance is 1e-9, which assures that the two values are the same
    within about 9 decimal digits. rel_tol must be greater than 0.0
    
    abs_tol: is a minimum absolute tolerance level -- useful for comparisons near zero.
    The name, isclose, is selected for consistency with the existing isnan and isinf .
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def fuse_rows(a, tol=1e-5, keep_closed=True):
    """
    fuses rows if the coordiante points are closer than tol but keeps the set 
    of coordinates closed if keep_closed is selected 
    
    Parameters
    --------
    a : ndarray
        input array of coordinate points
    tol : float
         tolerance distance
    keep_closed : bool
        to keep the set of coordinates closed after the fuse
        
    Returns
    -------
    d : ndarray 
        returns the fused array
    
    """
    tree = cKDTree(a)
    # print(a)
    rows_to_fuse = tree.query_pairs(tol, output_type="ndarray")
    # print(rows_to_fuse)
    # to delete but respect first and last rows:
    b = np.unique(rows_to_fuse[:, 1])

    if keep_closed and len(b) > 0 and b[-1] == len(a) - 1:
        c = b[:-1]
    else:
        c = b

    d = np.asarray([e for i, e in enumerate(a) if i not in c])
    return d


def Polygon_orientation(npArray):
    # Calculate of Polygon.
    area = PolygonArea(npArray)
    # print area
    if area > 0:  # counter-clockwise=True
        return True
    elif area < 0:  # clockwise=False
        return False
    else:
        return None


def PolygonArea(corners):
    """ Area of Polygon using Shoelace formula
    http://en.wikipedia.org/wiki/Shoelace_formula
     corners must be ordered in clockwise or counter-clockwise direction
    """

    n = len(np.atleast_1d(corners))  # of corners
    if n >> 1:
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += corners[i][0] * corners[j][1]
            area -= corners[j][0] * corners[i][1]
        area = area / 2.0
    else:
        area = 0
    return area


def TColgp_HArray1OfPnt_from_nparray(data):
    # Create TColgp_HArray1OfPnt from np.array
    # data[0] = x coord,  data[1] = y coord,  data[2] = z coord
    harray = TColgp_HArray1OfPnt(1, np.ma.size(data, 1))
    for index in range(0, np.ma.size(data, 1)):
        i = index + 1
        harray.SetValue(i, gp_Pnt(float(data[0, index]), float(data[1, index]), float(data[2, index])))
    return harray


def TColgp_HArray1OfPnt2d_from_nparray(data):
    # Create TColgp_HArray1OfPnt from np.array
    # data[0] = x coord,  data[1] = y coord,
    harray = TColgp_HArray1OfPnt2d(1, np.ma.size(data, 1))
    for index in range(0, np.ma.size(data, 1)):
        i = index + 1
        harray.SetValue(i, gp_Pnt2d(float(data[0, index]), float(data[1, index])))
    return harray


def point2d_list_to_TColgp_Array1OfPnt2d(li):
    return _Tcol_dim_1(li, TColgp_Array1OfPnt2d)


def _Tcol_dim_1(li, _type):
    pts = _type(0, len(li) - 1)
    for n, i in enumerate(li):
        pts.SetValue(n, i)
    return pts


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def calc_angle_between(v1, v2):
    """Returns the angle in degrees between vectors 'v1' and 'v2'::"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    if angle < 0:
        angle = 360 + angle
    return angle


def calc_DCT_angles(DCT_data):
    temp = []
    if np.array_equal(DCT_data[0], DCT_data[-1]):  # closed
        for i in range(0, DCT_data.shape[0]):
            if i == 0:  # first point
                v1 = DCT_data[i - 2] - DCT_data[i]
                v2 = DCT_data[i + 1] - DCT_data[i]
            elif i == DCT_data.shape[0] - 1:  # last point
                v1 = DCT_data[i - 1] - DCT_data[i]
                v2 = DCT_data[1] - DCT_data[i]
            else:
                v1 = DCT_data[i - 1] - DCT_data[i]
                v2 = DCT_data[i + 1] - DCT_data[i]
            # print(v1,v2)
            temp.append(calc_angle_between(v1, v2))

    else:  # open
        for i in range(0, DCT_data.shape[0]):
            if i == 0:  # first point
                v2 = DCT_data[i + 1] - DCT_data[i]
                v1 = -v2
            elif i == DCT_data.shape[0] - 1:  # last point
                v1 = DCT_data[i - 1] - DCT_data[i]
                v2 = -v1
            else:
                v1 = DCT_data[i - 1] - DCT_data[i]
                v2 = DCT_data[i + 1] - DCT_data[i]
            # print(v1,v2)
            temp.append(calc_angle_between(v1, v2))

    return np.array(temp)


def Pnt2dLst_to_npArray(Pnt2dLst):
    lst_tmp = []
    for i, item in enumerate(Pnt2dLst):
        lst_tmp.append([item.X(), item.Y()])
    return np.asarray(lst_tmp)


def PntLst_to_npArray(PntLst):
    lst_tmp = [p.Coord() for p in PntLst]
    return np.asarray(lst_tmp)


def Array_to_PntLst(array):
    PntLst = [gp_Pnt(a[0], a[1], a[2]) for a in array]
    return PntLst


def getID(custom):
    return custom.ID


def lin_pln_intersect(n0, p, p1, p2):
    """
    calculates the intersection point of the a line and plane. 
    The plane is defined by its normalized normal vector n0 and a plane point p.
    the line is defined from p1 to p2. The function not only returns the 
    intersection point but also the line coordinate lambda
    
    Parameters
    ---------
    n0 : array
    p : array
    p1 : array
    p2 : array
    
    Returns
    --------
    Pnt : array
    lamb : float 
    
    """
    n0 = np.asarray(n0)
    p = np.asarray(p)
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)

    v = p2 - p1
    d = np.dot(n0, p)
    lamb = (d - np.dot(n0, p1)) / (np.dot(n0, p2 - p1))
    pnt = p1 + lamb * v
    return (pnt, lamb)
