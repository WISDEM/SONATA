# PythonOCC Libraries
# Third party modules
from OCC.Core.Geom2d import Geom2d_BSplineCurve
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.TColgp import TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.Core.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.Display.SimpleGui import init_display

# First party modules
# Own Libraries:
from SONATA.cbm.topo.utils import (_Tcol_dim_1,)

# ==============================================================================
# STORE BSPLINE INFORMATION:
# ==============================================================================


class Para_Geom2d_BSplineCurve(object):
    """This Object stores the parameters as np.array,int,and boolean value, necessary to construct a Geom2d_BsplineCurve Opencascade Object"""

    def __init__(self, Geom2d_BSplineCurve):
        NbPoles = Geom2d_BSplineCurve.NbPoles()
        NbKnots = Geom2d_BSplineCurve.NbKnots()
        Poles = TColgp_Array1OfPnt2d(1, NbPoles)
        Weights = TColStd_Array1OfReal(1, NbPoles)
        Knots = TColStd_Array1OfReal(1, NbKnots)
        Multiplicities = TColStd_Array1OfInteger(1, NbKnots)

        Geom2d_BSplineCurve.Knots(Knots)
        Geom2d_BSplineCurve.Poles(Poles)
        Geom2d_BSplineCurve.Weights(Weights)
        Geom2d_BSplineCurve.Multiplicities(Multiplicities)

        self.Poles_bin = TColgp_Array1OfPnt2d_to_array(Poles)  # np.array
        self.Weights_bin = TColStd_to_array(Weights)  # np.array
        self.Knots_bin = TColStd_to_array(Knots)  # np.array
        self.Multiplicities_bin = TColStd_to_array(Multiplicities)  # np.array
        self.Degree_bin = Geom2d_BSplineCurve.Degree()  # np.array
        self.Periodic_bin = Geom2d_BSplineCurve.IsPeriodic()  # np.array

    def BSplineCurve2d(self):
        Poles = TColgp_Array1OfPnt2d_from_nparray(self.Poles_bin.T)
        Weights = _Tcol_dim_1(self.Weights_bin.tolist(), TColStd_Array1OfReal)
        Knots = _Tcol_dim_1(self.Knots_bin.tolist(), TColStd_Array1OfReal)

        Multiplicities = _Tcol_dim_1(self.Multiplicities_bin.tolist(), TColStd_Array1OfInteger)
        Degree = self.Degree_bin
        Periodic = self.Periodic_bin

        return Geom2d_BSplineCurve(Poles, Weights, Knots, Multiplicities, Degree, Periodic)


def ParaLst_from_BSplineLst(BSplineLst):
    Para_BSplineLst = []
    for item in BSplineLst:
        Para_BSpline = Para_Geom2d_BSplineCurve(item)
        Para_BSplineLst.append(Para_BSpline)

    return Para_BSplineLst


def BSplineLst_from_ParaLst(ParaLst):
    BSplineLst = []
    for item in ParaLst:
        BSplineLst.append(item.BSplineCurve2d())

    return BSplineLst


# =============================================================================
#            IF __MAIN__
# =============================================================================
if __name__ == "__main__":
    display, start_display, add_menu, add_function_to_menu = init_display("wx")
    display.Context.SetDeviationAngle(0.000001)  # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.000001)  # 0.001 default. Be careful to scale it to the problem.
    display.set_bg_gradient_color(20, 6, 111, 200, 200, 200)

    array = TColgp_Array1OfPnt2d(1, 5)
    array.SetValue(1, gp_Pnt2d(0, 0))
    array.SetValue(2, gp_Pnt2d(1, 2))
    array.SetValue(3, gp_Pnt2d(2, 3))
    array.SetValue(4, gp_Pnt2d(4, 3))
    array.SetValue(5, gp_Pnt2d(5, 5))
    bspline_1 = Geom2dAPI_PointsToBSpline(array).Curve()

    ParType = 2
    DegMin = 3
    DegMax = 8
    Continuity = 4
    Tol2D = 1.0e-3
    bspline_2 = Geom2dAPI_PointsToBSpline(array, ParType, DegMin, DegMax, Continuity, Tol2D).Curve()

    # the second one
    harray = TColgp_HArray1OfPnt2d(1, 5)
    harray.SetValue(1, gp_Pnt2d(0, 0))
    harray.SetValue(2, gp_Pnt2d(1, 2))
    harray.SetValue(3, gp_Pnt2d(2, 3))
    harray.SetValue(4, gp_Pnt2d(4, 3))
    harray.SetValue(5, gp_Pnt2d(5, 5))

    anInterpolation = Geom2dAPI_Interpolate(harray, False, 0.01)
    anInterpolation.Perform()
    bspline_3 = anInterpolation.Curve()

    for j in range(array.Lower(), array.Upper() + 1):
        p = array.Value(j)
    #    display.DisplayShape(p, update=False)
    for j in range(harray.Lower(), harray.Upper() + 1):
        p = harray.Value(j)
    #    display.DisplayShape(p, update=False)

    # ============Para_Geom2d_BSplineCurve======================================
    test = Para_Geom2d_BSplineCurve(bspline_3)
    bspline_4 = test.BSplineCurve2d()

    display.DisplayShape(bspline_1, update=False)
    display.DisplayShape(bspline_2, update=False, color="BLUE")
    display.DisplayShape(bspline_3, update=False, color="GREEN")
    display.DisplayShape(bspline_4, update=False, color="BLACK")

    display.View_Top()
    display.FitAll()
    start_display()
