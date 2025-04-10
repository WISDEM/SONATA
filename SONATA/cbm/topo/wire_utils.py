# Basic Libraries:
# Third party modules
from OCC.Core.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeWire,
                                     BRepBuilderAPI_Transform,)
from OCC.Core.GCPnts import (GCPnts_QuasiUniformAbscissa,)
from OCC.Core.Geom import Geom_Plane
# PythonOCC Libraries
from OCC.Core.gp import (gp_Pnt,gp_Dir, gp_Vec,)
from OCC.Core.TopoDS import topods
from OCC.Core.TopTools import (TopTools_ListOfShape,)

# First party modules
# Own Libraries:

###############################################################################
# Wire Utilities
###############################################################################


def equidistant_Points_on_wire(wire, NbPoints):
    AdaptorComp = BRepAdaptor_CompCurve(wire, True)
    discretization = GCPnts_QuasiUniformAbscissa(AdaptorComp, NbPoints)

    PntLst = []
    for j in range(1, NbPoints + 1):
        para = discretization.Parameter(j)
        Pnt = gp_Pnt()
        AdaptorComp.D0(para, Pnt)
        PntLst.append(Pnt)
    return PntLst


def trsf_wire(wire, gp_Trsf, copy=False):
    aBRespTrsf = BRepBuilderAPI_Transform(wire, gp_Trsf, copy)
    aBRespTrsf.Build()
    atrsfedShape = aBRespTrsf.Shape()
    trsfedWire = topods.Wire(atrsfedShape)
    return trsfedWire


def build_wire_from_BSplineLst2(BSplineLst):  # Builds TopoDS_Wire from connecting BSplineSegments and returns it
    tmp_wire = BRepBuilderAPI_MakeWire()
    for i, item in enumerate(BSplineLst):
        P = gp_Pnt(0, 0, 0)
        V = gp_Dir(gp_Vec(0, 0, 1))
        Plane = Geom_Plane(P, V)
        tmp_edge = BRepBuilderAPI_MakeEdge(item, Plane)
        # print(tmp_edge.IsDone())
        tmp_wire.Add(tmp_edge.Edge())
    return tmp_wire.Wire()


def build_wire_from_BSplineLst(BSplineLst, twoD=True):  # Builds TopoDS_Wire from connecting BSplineSegments and returns it
    P = gp_Pnt(0, 0, 0)
    V = gp_Dir(gp_Vec(0, 0, 1))
    Plane = Geom_Plane(P, V)
    TopTools_EdgeLst = TopTools_ListOfShape()
    for i, item in enumerate(BSplineLst):
        if twoD:
            tmp_edge = BRepBuilderAPI_MakeEdge(item, Plane)
        else:
            tmp_edge = BRepBuilderAPI_MakeEdge(item)
        # print(tmp_edge.IsDone())
        TopTools_EdgeLst.Append(tmp_edge.Edge())

    wire = BRepBuilderAPI_MakeWire()
    wire.Add(TopTools_EdgeLst)
    wire.Build()
    return wire.Wire()
