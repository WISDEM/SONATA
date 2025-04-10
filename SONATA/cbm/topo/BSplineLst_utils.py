# Basic Libraries:
# Core Library modules

# Third party modules
import matplotlib.pyplot as plt
import numpy as np
from OCC.Core.GCPnts import (GCPnts_AbscissaPoint, 
                             GCPnts_QuasiUniformDeflection,)
from OCC.Core.Geom import Geom_BSplineCurve, Geom_Curve
from OCC.Core.Geom2d import (Geom2d_BSplineCurve, Geom2d_Curve, Geom2d_Line,
                             Handle_Geom2d_BSplineCurve_DownCast,)
from OCC.Core.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Core.Geom2dAPI import (Geom2dAPI_Interpolate,
                                Geom2dAPI_ProjectPointOnCurve,Geom2dAPI_InterCurveCurve,)
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.GeomAPI import (GeomAPI_IntCS, GeomAPI_Interpolate,
                              GeomAPI_ProjectPointOnCurve,)
# PythonOCC Libraries
from OCC.Core.gp import (gp_Pnt,
                         gp_Pnt2d, gp_Vec, gp_Vec2d,)

# First party modules
# Own Libraries:
from SONATA.cbm.topo.utils import (Pnt2dLst_to_npArray,
                                   TColgp_HArray1OfPnt2d_from_nparray,
                                   TColgp_HArray1OfPnt_from_nparray,
                                   calc_DCT_angles, fuse_rows, isclose,
                                   unique_rows,)

###############################################################################
# BSpline and BSplineLst Utilities
###############################################################################

def ProjectPointOnBSplineLst(BSplineLst, Pnt2d, tolerance_distance=100):
    """
    Projets a Pnt2d onto a list of BSplines with given tolerance_distance and 
    returns the a tuple of information of the closest Point
    
    Parameters
    ----------
    BSplineLst : list
        consecutive list of Geom2d_BSplineCurves 
    
    Pnt2d : gp_Pnt2d
        2d Point object of OCC.gp.gp_Pnt2d
    
    Returns
    ----------
    p2 : list 
     list of projection information [gp_Pnt2d, index of BSplineLst, curve parameter u, distance]
         
    
    """
    p2 = []
    for idx, item in enumerate(BSplineLst):
        projection2 = Geom2dAPI_ProjectPointOnCurve(Pnt2d, item)
        for j in range(1, projection2.NbPoints() + 1):
            if projection2.Distance(j) <= tolerance_distance:
                p2.append([projection2.Point(j), idx, projection2.Parameter(j), projection2.Distance(j)])
            else:
                None

    if p2:
        p2 = np.asarray(p2)
        min_index = p2[:, 3].argmin()
        return p2[min_index, :]
    else:
        return []

def findPnt_on_curve(Pnt, Curve):
    if isinstance(Pnt, gp_Pnt2d) and issubclass(Curve.__class__, Geom2d_Curve):
        twoD = True
    elif isinstance(Pnt, gp_Pnt) and issubclass(Curve.__class__, Geom_Curve):
        twoD = False

    if twoD:
        projection = Geom2dAPI_ProjectPointOnCurve(Pnt, Curve)
    else:
        projection = GeomAPI_ProjectPointOnCurve(Pnt, Curve)
    u = projection.LowerDistanceParameter()
    return u

def isPnt_on_curve(Pnt, Curve, tolerance=1e-6):
    """
    checks if a point (either gp_Pnt or gp_Pnt2d) is on a curve (subclass of 
    Geom2d_Curve or Geom_Curve) by projecting the Point orthogonal on to the 
    curve and checing ifthe distance is below a tolerance tol.
    
    Parameter
    ---------
    Pnt : OCC.gp_Pnt or OCC.gp_Pnt2d
    Curve : OCC.Geom_Curve or OCC.Geom_Curve Handle
    tolerance : float, optional
        default tolerance 1e-6
        
    Return
    --------
    Trigger : bool    
    """
    if isinstance(Pnt, gp_Pnt2d) and issubclass(Curve.__class__, Geom2d_Curve):
        twoD = True
    elif isinstance(Pnt, gp_Pnt) and issubclass(Curve.__class__, Geom_Curve):
        twoD = False

    if twoD:
        projection = Geom2dAPI_ProjectPointOnCurve(Pnt, Curve)
    else:
        projection = GeomAPI_ProjectPointOnCurve(Pnt, Curve)
    Trigger = False
    for i in range(1, projection.NbPoints() + 1):
        if projection.Distance(i) <= tolerance:
            Trigger = True
        else:
            None
    return Trigger

def findPnt_on_BSplineLst(Pnt, BSplineLst):
    for i, item in enumerate(BSplineLst):
        if isPnt_on_curve(Pnt, item):
            u = findPnt_on_curve(Pnt, item)
            coordinates = [i, u]
            break
        else:
            coordinates = None
    return coordinates


def get_BSpline_length(BSpline):
    """
    Returns the Length of the BSpline
    """
    if isinstance(BSpline, Geom_BSplineCurve):
        twoD = False
    elif isinstance(BSpline, Geom2d_BSplineCurve):
        twoD = True

    tolerance = 1e-10
    first = BSpline.FirstParameter()
    last = BSpline.LastParameter()

    if twoD:
        Adaptor = Geom2dAdaptor_Curve(BSpline)
    else:
        Adaptor = GeomAdaptor_Curve(BSpline)

    Length = GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance)
    return Length


def copy_BSpline(BSpline):
    BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(BSpline.Copy())
    return BSplineCopy


def copy_BSplineLst(BSplineLst):
    copied_BSplineLst = []
    for i, item in enumerate(BSplineLst):
        copied_BSplineLst.append(copy_BSpline(item))
    return copied_BSplineLst


def get_BSplineLst_length(BSplineLst):

    """
    Returns the cummulated length of the BSplineLst
    """

    CummLength = 0
    for i, item in enumerate(BSplineLst):
        CummLength += get_BSpline_length(item)
    return CummLength


def find_BSplineLst_coordinate(BSplineLst, S, start, end):
    """The function find_BSplineLst_coordinate returns the list index of 
    the bspline and its parameter where the coordinate is located.
    
    Parameters
    ----------
    S : float
        Coordinate to be found 
    start : float
        start of the interval where the BSplineLst is defined
    end : float
        end of the interval where the BSplineLst is defined
    
    Returns
    ----------
    [i,U]: list
        Return [index of the bspline, parameter U on the bspline]
        
    """

    if all(isinstance(s, Geom_BSplineCurve) for s in BSplineLst):
        twoD = False
    elif all(isinstance(s, Geom2d_BSplineCurve) for s in BSplineLst):
        twoD = True

    BSplineLstLength = get_BSplineLst_length(BSplineLst)
    x = BSplineLstLength * (S - start) / (end - start)
    CummLength = 0
    tolerance = 1e-8
    for i, item in enumerate(BSplineLst):
        first = item.FirstParameter()
        last = item.LastParameter()

        if twoD:
            Adaptor = Geom2dAdaptor_Curve(item)
        else:
            Adaptor = GeomAdaptor_Curve(item)

        CummLength += get_BSpline_length(item)
        if x < CummLength or isclose(x, CummLength):
            dist = x - (CummLength - GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance))
            tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, dist, first)
            U = tmp.Parameter()
            break
    return [i, U]


def reverse_BSplineLst(BSplineLst):
    BSplineLst = list(reversed(BSplineLst))
    for i, item in enumerate(BSplineLst):
        item.Reverse()
    return BSplineLst


def get_BSplineLst_Pnt2d(BSplineLst, S, start, end):
    if start > end:  # transfrom the interval if start>end to a basic interval [0..]
        D = 1 - start
        S = (S + D) % 1
        end = (end + D) % 1
        start = (start + D) % 1

    P = gp_Pnt2d()
    [idx, U] = find_BSplineLst_coordinate(BSplineLst, S, start, end)
    BSplineLst[idx].D0(U, P)
    return P


def get_BSplineLst_D2(BSplineLst, S, start, end, **kwargs):
    if start > end:  # transfrom the interval if start>end to a basic interval [0..]
        D = 1 - start
        S = (S + D) % 1
        end = (end + D) % 1
        start = (start + D) % 1

    if all(isinstance(s, Geom_BSplineCurve) for s in BSplineLst):
        twoD = False
    elif all(isinstance(s, Geom2d_BSplineCurve) for s in BSplineLst):
        twoD = True

    if twoD:
        P = gp_Pnt2d()
        V1 = gp_Vec2d()
        V2 = gp_Vec2d()
    else:
        P = gp_Pnt()
        V1 = gp_Vec()
        V2 = gp_Vec()

    [idx, U] = find_BSplineLst_coordinate(BSplineLst, S, start, end, **kwargs)
    BSplineLst[idx].D2(U, P, V1, V2)
    return P, V1, V2


def intersect_BSplineLst_with_plane(BSplineLst, plane):
    """
    intersects a 3D-BSplineLst (OCC.Geom_BSplineCurve) with a plane (Handle<Geom_Surface>) and returns the
    a list of intersection coordinates and the corresponding list of 
    intersection points
    
    Paramters
    ---------
    BSplineLst : list
         list of OCC.Geom_BSplineCurve
    plane : <Geom_Surface>
         Plane or surface - Geom_Surface, Geom_Plane is a subclass of Geom_Surface
    
    Return
    --------
    IntCoords : list
        list of intersection coordinates
        [index of the bspline, parameter U on the bspline]
    IntPnts : list 
        list of intersection points (OCC.gp_Pnt)
    """
    # plane = Geom_Plane(gp_Pnt(float(x),0,0), gp_Dir(1,0,0))
    IntCoords = []
    IntPnts = []
    for i, item in enumerate(BSplineLst):
        intersection = GeomAPI_IntCS(item, plane)
        for j in range(1, intersection.NbPoints() + 1):
            IntPnt = intersection.Point(j)
            IntPnts.append(IntPnt)
            u = findPnt_on_curve(IntPnt, item)
            IntCoords.append([i, u])

    return IntCoords, IntPnts


def trim_BSplineLst(BSplineLst, S1, S2, start, end):
    """ the trim_BSplineLst function trims the BSplineLst that is defined on 
    the interval between start and end to the interval defined by S1 and S2.
    
    Parameters
    ----------
    BSplineLst: list of geom2d_bsplines
        is the object to be modified.
    S1: float
        start of the interval to which the object is trimmed.
    S2: float
        end of the interval to which the object is trimmed.
    start: float
        start of the interval where the BSplineLst is defined
    end: float
        end of the interval where the BSplineLst is defined
            
    Returns
    ----------
    trimmed_BSplineLst: list of geom2d_bsplines
        is the modified object.
            
    Example
    ----------
    a BSplineLst is defined between 0.0 and 1.0 and shall be trimmed 
    to a interval S1=0.3 and S2=0.7
    """

    if start > end:  # transfrom the interval if start>end to a basic interval [0..]
        D = 1 - start
        S1 = (S1 + D) % 1
        S2 = (S2 + D) % 1
        end = (end + D) % 1
        start = (start + D) % 1

    if S1 > S2:
        trimmed_BSplineLst = []
        front_BSplineLst = []
        rear_BSplineLst = []
        para1 = find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 = find_BSplineLst_coordinate(BSplineLst, S2, start, end)

        # check if para1 and para2 is out of bounds!
        if isclose(para1[1], BSplineLst[para1[0]].FirstParameter()):
            para1[1] = BSplineLst[para1[0]].FirstParameter()
        if isclose(para2[1], BSplineLst[para2[0]].LastParameter()):
            para2[1] = BSplineLst[para2[0]].LastParameter()

        for i, item in enumerate(BSplineLst):
            First = item.FirstParameter()
            Last = item.LastParameter()
            if para1[0] == i and para2[0] != i:  # Okay
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                if isclose(para1[1], Last):
                    pass
                elif isclose(para1[1], First):
                    BSplineCopy.Segment(First, Last)
                    trimmed_BSplineLst.append(BSplineCopy)
                else:
                    BSplineCopy.Segment(para1[1], Last)
                    trimmed_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                front_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                rear_BSplineLst.append(BSplineCopy)

            elif para1[0] == i and para2[0] == i:  # Okay
                BSplineCopy1 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy1.Segment(para1[1], Last)
                trimmed_BSplineLst.append(BSplineCopy1)
                BSplineCopy2 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy2.Segment(First, para2[1])
                trimmed_BSplineLst.append(BSplineCopy2)
                break

            elif para1[0] != i and para2[0] == i:  # Okay
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, para2[1])
                front_BSplineLst.append(BSplineCopy)
                # break

        if front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst + front_BSplineLst

        elif front_BSplineLst and not rear_BSplineLst:
            trimmed_BSplineLst = front_BSplineLst

        elif not front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst

    elif S1 == start and S2 == end:
        trimmed_BSplineLst = copy_BSplineLst(BSplineLst)

    elif S2 == start or S1 == S2:
        trimmed_BSplineLst = []

    elif S2 > S1:
        trimmed_BSplineLst = []

        para1 = find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 = find_BSplineLst_coordinate(BSplineLst, S2, start, end)

        # check if para1 and para2 is out of bounds!
        if isclose(para1[1], BSplineLst[para1[0]].FirstParameter()):
            para1[1] = BSplineLst[para1[0]].FirstParameter()
        if isclose(para2[1], BSplineLst[para2[0]].LastParameter()):
            para2[1] = BSplineLst[para2[0]].LastParameter()

        for i, item in enumerate(BSplineLst):
            First = item.FirstParameter()
            Last = item.LastParameter()
            if para1[0] == i and para2[0] != i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                if isclose(para1[1], Last):
                    pass
                elif isclose(para1[1], First):
                    BSplineCopy.Segment(First, Last)
                    trimmed_BSplineLst.append(BSplineCopy)
                else:
                    BSplineCopy.Segment(para1[1], Last)
                    trimmed_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                trimmed_BSplineLst.append(BSplineCopy)

            elif para1[0] == i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(para1[1], para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
                break

            elif para1[0] != i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                if isclose(para2[1], Last):
                    trimmed_BSplineLst.append(BSplineCopy)
                else:
                    BSplineCopy.Segment(First, para2[1])
                    trimmed_BSplineLst.append(BSplineCopy)
                break

    return trimmed_BSplineLst


def discretize_BSplineLst(BSplineLst, Deflection=2e-4, AngularDeflection=0.02, CurvatureDeflection=0.06, MinimumOfPoints=50):
    Pnt2dLst = []
    for i, item in enumerate(BSplineLst):
        Adaptor = Geom2dAdaptor_Curve(item)
        discretization = GCPnts_QuasiUniformDeflection(Adaptor, Deflection, 1)  # GeomAbs_Shape Continuity: 1=C0, 2=G1, 3=C1, 3=G2,...
        # discretization = GCPnts_TangentialDeflection(Adaptor,AngularDeflection,CurvatureDeflection,MinimumOfPoints)
        NbPoints = discretization.NbPoints()
        for j in range(1, NbPoints + 1):
            Pnt = discretization.Value(j)
            Pnt2dLst.append(Pnt)

    a = Pnt2dLst_to_npArray(Pnt2dLst)
    b = np.around(a, 10)  # Evenly round to the given number of decimals.

    # check if it is closed:
    if np.array_equal(b[0], b[-1]):
        closed = True
    else:
        closed = False

    npArray = unique_rows(b)  # Remove possible doubles!

    if closed:
        npArray = np.vstack((npArray, b[0]))
    else:
        None
    return npArray


def BSplineLst_from_dct(DCT_data, angular_deflection=15,  cutoff_style = 2, closed=False, tol_interp=1e-6, twoD=True):
    if closed and not np.allclose(DCT_data[0], DCT_data[-1]):
        print("INFO:\t Closing open discrete definition")
        DCT_data = np.concatenate((DCT_data, DCT_data[0:1, :]), axis=0)

    if cutoff_style == 0:
        DCT_data = fuse_rows(DCT_data, 1e-3)  # check all datapoints and merge if allclose is true
    else:
        DCT_data = fuse_rows(DCT_data,1e-6)

    # Find corners and edges of data
    DCT_angles = calc_DCT_angles(DCT_data)
    corners = []
    for i in range(0, DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection):
            corners.append(i)
    NbCorners = np.size(corners)

    # Segmenting data according to corners
    # ====================================
    DCT_Segments = []

    # Closed:
    if np.allclose(DCT_data[0], DCT_data[-1]):
        closed = True
        if NbCorners == 0:
            DCT_Segments.append(DCT_data)

        elif NbCorners > 0:

            if corners[0] == 0 and corners[-1] == len(DCT_data) - 1:
                for i in range(0, NbCorners - 2):
                    DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
                DCT_Segments.append(DCT_data[corners[-2] : corners[-1] + 1])  # last segment

            # regular case:
            else:
                DCT_Segments.append(DCT_data[: corners[0] + 1])  # first
                if NbCorners > 1:  # intermediate segments:
                    for i in range(0, NbCorners - 1):
                        DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
                DCT_Segments.append(DCT_data[corners[-1] :])  # last
                DCT_Segments[0] = np.concatenate((DCT_Segments[-1][:-1], DCT_Segments[0]))
                del DCT_Segments[-1]

    else:  # Open:
        closed = False
        if NbCorners == 0:
            DCT_Segments.append(DCT_data)

        if NbCorners > 0:
            DCT_Segments.append(DCT_data[: corners[0] + 1])  # first
            if NbCorners > 1:  # intermediate segments:
                for i in range(0, NbCorners - 1):
                    DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
            DCT_Segments.append(DCT_data[corners[-1] :])  # last segment

    list_of_bsplines = []
    for i, item in enumerate(DCT_Segments):
        data = item.T

        if twoD:
            tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
        else:
            tmp_harray = TColgp_HArray1OfPnt_from_nparray(data)

        try:
            if twoD:
                tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray, False, tol_interp)  # Interpolate datapoints to bspline
            else:
                tmp_interpolation = GeomAPI_Interpolate(tmp_harray, False, tol_interp)
        except:
            print("ERROR: BSplineLst_from_dct did not perform the Interpolation with tol_interp", tol_interp)  # RAISED ERROR
            plt.figure()
            plt.clf()
            plt.axis("equal")
            plt.plot(*item.T, marker=".")
            for i, it in enumerate(item):
                plt.annotate(i, (it[0], it[1]), color="black")
            plt.show()

        tmp_interpolation.Perform()
        tmp_bspline = tmp_interpolation.Curve()
        list_of_bsplines.append(tmp_bspline)

    return list_of_bsplines


def set_BSplineLst_to_Origin2(BSplineLst, gp_Pnt2d, tol=1e-1):
    """
    this procedure reorders the self.BSplineLst to an origin if the layer 
    is closed. The Origin is detected by searching for an orthogonal 
    projection of the StartPoint of the self.Boundary_BSplineLst. If no 
    projection is found it takes the closest neighbor of the discrete 
    offlinepts (Offset Line Points).
            
    """
    # Determine Origin as point
    if BSplineLst[0].StartPoint().IsEqual(BSplineLst[-1].EndPoint(), 1e-5):
        proj = ProjectPointOnBSplineLst(BSplineLst, gp_Pnt2d, tol)
        if len(proj) > 0:
            OriPnt = proj[0]
        elif len(proj) == 0:
            print("WARNING : No Projection found to produce origin!")

    else:
        print("ERRROR: Only apply set_BSplineLst_to_Origin2 function when BSplineLst is Closed!")

    # Reorder Sequence of BSplines of BSplinesLst
    OriPara = findPnt_on_BSplineLst(OriPnt, BSplineLst)
    OBSplineLst = []
    CorrectOrigin = False

    for i, item in enumerate(BSplineLst):
        if i == OriPara[0]:
            First = item.FirstParameter()
            Last = item.LastParameter()

            if isclose(OriPara[1], First) == True:
                OBSplineLst.append(item)
                CorrectOrigin = True

            elif isclose(OriPara[1], Last) == True:
                CorrectOrigin = False
                BSplineCurve2 = item

            else:
                CorrectOrigin = False
                BSplineCurve1 = copy_BSpline(item)
                BSplineCurve1.Segment(OriPara[1], Last)
                BSplineCurve2 = copy_BSpline(item)
                BSplineCurve2.Segment(First, OriPara[1])
                OBSplineLst.append(BSplineCurve1)

        elif i > OriPara[0]:
            OBSplineLst.append(item)
        else:
            None

    for i, item in enumerate(BSplineLst):
        if i < OriPara[0]:
            OBSplineLst.append(item)
        else:
            None

    if CorrectOrigin == False:
        OBSplineLst.append(BSplineCurve2)
    else:
        None

    BSplineLst = OBSplineLst

    return BSplineLst


def intersect_BSplineLst_with_BSpline(BSplineLst, BSpline, tolerance=1e-10):
    IntPnts = []
    IntPnts_Pnt2d = []
    for i, item in enumerate(BSplineLst):
        intersection = Geom2dAPI_InterCurveCurve(BSpline, item, tolerance)
        for j in range(1, intersection.NbPoints() + 1):
            IntPnt = intersection.Point(j)
            IntPnts_Pnt2d.append(IntPnt)
            u = findPnt_on_curve(IntPnt, item)
            IntPnts.append([i, u])
    return IntPnts, IntPnts_Pnt2d

if __name__ == "__main__":
    pass
