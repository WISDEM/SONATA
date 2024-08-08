import itertools
import matplotlib.pyplot as plt

# Third party modules
import numpy as np
from OCC.Core.Geom2dAPI import (Geom2dAPI_PointsToBSpline,
                                Geom2dAPI_ProjectPointOnCurve,)
from OCC.Core.gp import gp_Pnt2d, gp_Vec2d, gp_Dir2d
from OCC.Core.Geom2d import Geom2d_Line
from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve

# First party modules
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.node import Node
from SONATA.cbm.topo.BSplineLst_utils import (
    ProjectPointOnBSplineLst, find_BSplineLst_coordinate, intersect_BSplineLst_with_BSpline,)
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d

def display_bsplinelst(bsplinelst, color = 'red'):
    for bspline in bsplinelst:
        u_min, u_max = bspline.FirstParameter(), bspline.LastParameter()
        # Extract points for plotting
        num_points = 100  # Number of points to plot
        u_values = [u_min + (u_max - u_min) * i / (num_points - 1) for i in range(num_points)]
        x_values = [bspline.Value(u).X() for u in u_values]
        y_values = [bspline.Value(u).Y() for u in u_values]
        plt.plot(x_values, y_values, color = color)

def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'"""
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Clip to handle numerical errors
    return np.arccos(cos_theta)


def three_pnt_projection_method(points, bspline, tolerance):
    """
    Project node onto target Bspline along the bisecting vector based on the two adjacent nodes.
    ----------
    Input:
    points: array of three points in order. The middle point is the point being projected and the
    first and last points are the two points adjacent to the point being projected.
    bspline: The targeted bspline.
    tolerance: A distance tolerance. If the point is projected a distance farther than the tolerance,
    it is not projected.

    Output:
    The point at which the bisecting vector intersects the target bspline.
    """
    prev = points[0]
    curr = points[1]
    start_point_2d = gp_Pnt2d(curr[0],curr[1])
    next = points[2]
    vec1 = curr - prev
    vec2 = next - curr
    angle = angle_between(vec1,vec2)
    
    # Find the two normal vectors pointing from the middle point to the adjacent points
    vec1_perp = np.array([-vec1[1], vec1[0]]) / np.linalg.norm([-vec1[1], vec1[0]])
    vec2_perp = np.array([-vec2[1], vec2[0]]) / np.linalg.norm([-vec2[1], vec2[0]])
    
    # Finding the vector that bisects the two normal vectors. The project point will be in the direction of the bisecting vector.
    bisecting_vector = - np.sqrt((1 + np.dot(vec1_perp, vec2_perp)) / 2) * (vec1_perp + vec2_perp) / np.linalg.norm(vec1_perp + vec2_perp)
    direction_2d = gp_Dir2d(bisecting_vector[0],bisecting_vector[1])
    line_2d = Geom2d_Line(start_point_2d, direction_2d)

    # Finding where the bisecting vector intersects the target Bspline
    intersector = Geom2dAPI_InterCurveCurve(line_2d, bspline)
    if intersector.NbPoints() > 0:
        intersection_point = intersector.Point(1)
        distance = np.sqrt((curr[0]-intersection_point.X())**2+(curr[1]-intersection_point.Y())**2)
        if distance < tolerance*(2+np.cos(angle)):
            return intersection_point
    return None

def mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst, a_nodes, b_BSplineLst, layer_thickness, tol=1e-2, crit_angle=95, LayerID=0, refL=1.0, **kw):
    """
    *function to mesh the SONATA topologies by projecting nodes onto the 
    generated BSplineLists.   
    
    :rtype: [list,list,list]
    
    Variables and arguments:
    ----------
    a_BSplineLst & a_nodes: the projection source
    b_BSplineLst & b_nodes: the projection destination curve and its resulting nodes
    layer_thickness & tol: is exactly that and will be used with the tolerance(tol) to determine
                     a distance, in which the resulting projection point has to be.
    crit_angle: is the critical angle to determine a corner if 2 projection points are found.
    display: the kwargs display object can be passed to plot/display within the main OCC3DViewer
                    
    
    Workflow:
    ----------        
      
    
    Notes & Comments:
    ----------  
    TODO: * scale distance not only to layerthickenss but also to min_len.
            or adapt the distance individually for each node. 
    """
    # display_bsplinelst(a_BSplineLst, 'blue')
    # display_bsplinelst(b_BSplineLst, 'green')

    # KWARGS:
    if kw.get("display") != None:
        display = kw.get("display")

    # LayerID = 'T_' + a_nodes[0].parameters[0]
    # TODO:
    b_nodes = []
    cellLst = []
    distance = (1 + tol) * layer_thickness
    flag_integrate_leftover_interior_nodes = False

    # Is a_BSplineLst closed?
    
    closed_a = False
    if a_BSplineLst[0].StartPoint().IsEqual(a_BSplineLst[-1].EndPoint(), 1e-5):
        closed_a = True

    # ==================PROJECT POINTS ON LOWER BOUNDARY =======================
    if closed_a == True:
        prj_nodes = a_nodes
    else:
        prj_nodes = a_nodes[1:-1]

    for i, node in enumerate(prj_nodes, start=1):
        Pnt2d = node.Pnt2d
        pPnts = []
        pPara = []
        pIdx = []
        projected = False

        # Projects current node onto the target Bspline where the vector between the current node
        # and the projected node is perpendicular.
        for idx, item in enumerate(b_BSplineLst):
            first = item.FirstParameter()
            last = item.LastParameter()
            tol = 2e-6 * refL
            Umin = first - (last - first) * tol
            Umax = last + (last - first) * tol
            projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d, item, Umin, Umax)

            for j in range(1, projection.NbPoints() + 1):
                if projection.Distance(j) <= distance:
                    pPnts.append(projection.Point(j))
                    pPara.append(projection.Parameter(j))
                    pIdx.append(idx)
                    projected = True
                else:
                    None
        # If no node is found that creates a vector perpendicular to the target Bspline, a different
        # projection method is tried using the current node and the neighboring nodes to create a
        # bisecting vector. The point at which the bisecting vector intersects the target Bspline
        # is where the projected point is placed. This helps for projecting nodes on corners.
        if not projected and i-1 != 0 and not i-1 > len(prj_nodes)-2:
            for idx, item in enumerate(b_BSplineLst):
                prev_point = np.array([prj_nodes[i-2].coordinates[0], prj_nodes[i-2].coordinates[1]])
                curr_point = np.array([prj_nodes[i-1].coordinates[0], prj_nodes[i-1].coordinates[1]])
                next_point = np.array([prj_nodes[i].coordinates[0], prj_nodes[i].coordinates[1]])
                projected_point = three_pnt_projection_method([prev_point,curr_point,next_point], item, distance*10)
                if projected_point:
                    pPnts.append(projected_point)
                    pPara.append(node.parameters[2])
                    pIdx.append(idx)
                    node.corner = False
                    break


        # ==================making sure the pPnts are unique:
        """It happend that somehow the same points were found multiple times"""
        unique_tol = 5e-5 * refL
        rm_idx = []
        for a, b in itertools.combinations(enumerate(pPnts), 2):
            if a[1].IsEqual(b[1], unique_tol):
                rm_idx.append(a[0])

        pPnts = [j for k, j in enumerate(pPnts) if k not in rm_idx]
        pPara = [j for k, j in enumerate(pPara) if k not in rm_idx]
        pIdx = [j for k, j in enumerate(pIdx) if k not in rm_idx]

        # =========if 3 Points are found. Select the 2 points that create the larger angle.
        if len(pPnts) == 3:
            v0 = gp_Vec2d(Pnt2d, pPnts[0])
            v1 = gp_Vec2d(Pnt2d, pPnts[1])
            v2 = gp_Vec2d(Pnt2d, pPnts[2])
            angle01 = abs(v0.Angle(v1) * 180 / np.pi)
            angle02 = abs(v0.Angle(v2) * 180 / np.pi)
            angle12 = abs(v1.Angle(v2) * 180 / np.pi)
            aLst = [angle01, angle02, angle12]
            pPnts.pop(aLst.index(max(aLst)))
            pPara.pop(aLst.index(max(aLst)))
            pIdx.pop(aLst.index(max(aLst)))

        # ==================DETECT CORNERS======================================
        if len(pPnts) == 0:
            print("ERROR:\t No Projection found for node,", node.id)

        elif len(pPnts) == 1:
            b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))

        elif len(pPnts) == 2:
            # =======================determine the angle of the potential corner
            v1 = gp_Vec2d(Pnt2d, pPnts[0])
            v2 = gp_Vec2d(Pnt2d, pPnts[1])
            angle = 180 - abs(v1.Angle(v2) * 180 / np.pi)
            vres = v1.Added(v2)

            if angle < crit_angle and not node.corner:
                node.corner = True

            if angle >= crit_angle:
                node.corner = False
            # print 'Node:', node.id, 'corner angle: ', angle

            # =======================determine the corner character - regular_corner type:(Bool)
            node.regular_corner = True
            vp0 = gp_Vec2d()
            vp1 = gp_Vec2d()
            p_tmp = gp_Pnt2d()
            b_BSplineLst[pIdx[0]].D1(pPara[0], p_tmp, vp0)
            b_BSplineLst[pIdx[1]].D1(pPara[1], p_tmp, vp1)
            dot0 = vres.Dot(vp0)
            # dot1 = vres.Dot(vp1)

            if dot0 < 0:
                node.regular_corner = False
                z_BSplineLst = b_BSplineLst[pIdx[1] :] + b_BSplineLst[: pIdx[0] + 1]

            elif dot0 > 0:
                node.regular_corner = True

            else:
                print("ERROR: cannot determine regular_corner because vres and v0 are orthogonal")

            # =======================determine the exterior corners on z_BsplineLst
            # TODO: DETECT ALL EXTERIOR CORNERS WITHIN THAT INTERVAL pIdx[0],pPara[0],'|',pIdx[1], pPara[1]
            # TODO: IF NO EXTERIOR CORNER IS FOUND USE BISECTOR!
            # TODO: USE node.cornertype = 1,2,3,4 to detemine Shape of Elements
            # TODO: THIS ALGORITHM DOESNT FIND ALL
            exterior_corners = []
            exterior_corners_para = []
            aglTol = 5.0
            if node.regular_corner == True:
                for j, item in enumerate(b_BSplineLst[pIdx[0] : pIdx[1]], start=pIdx[0]):
                    spline1 = item
                    # print node, pIdx[0], pIdx[1]
                    spline2 = b_BSplineLst[j + 1]
                    u1, p1, v1 = spline1.LastParameter(), gp_Pnt2d(), gp_Vec2d()
                    u2, p2, v2 = spline2.FirstParameter(), gp_Pnt2d(), gp_Vec2d()
                    spline1.D1(u1, p1, v1)
                    spline2.D1(u2, p2, v2)

                    Angle = abs(v1.Angle(v2)) * 180 / np.pi
                    if Angle > aglTol:
                        exterior_corners.append(item.EndPoint())
                        exterior_corners_para.append([LayerID, j, u1])
                        # display.DisplayShape(item.EndPoint(),color='WHITE')

            else:
                for j, item in enumerate(z_BSplineLst[:-1]):
                    spline1 = item
                    spline2 = z_BSplineLst[j + 1]
                    u1, p1, v1 = spline1.LastParameter(), gp_Pnt2d(), gp_Vec2d()
                    u2, p2, v2 = spline2.FirstParameter(), gp_Pnt2d(), gp_Vec2d()
                    spline1.D1(u1, p1, v1)
                    spline2.D1(u2, p2, v2)

                    Angle = abs(v1.Angle(v2)) * 180 / np.pi
                    if Angle > aglTol:
                        exterior_corners.append(item.EndPoint())
                        if len(b_BSplineLst) > j + pIdx[1]:
                            idx = j + pIdx[1]
                        else:
                            idx = j + pIdx[1] - len(b_BSplineLst)
                        exterior_corners_para.append([LayerID, idx, u1])
                        # display.DisplayShape(item.EndPoint(),color='WHITE')            

            # =======================generate b_nodes===========================
            # print node,'corner: ', node.corner, ', regular_corner = ',node.regular_corner, ',  Len:exterior_corners =',len(exterior_corners)

            # ===CORNERSTYLE 0======
            if len(exterior_corners) == 0 and node.corner == False:
                node.cornerstyle = 0
                # TODO: a more robust possibilit is to use a bisector vres

                if pIdx[0] == pIdx[1]:  #
                    v = gp_Vec2d(pPnts[0], pPnts[1])
                    cP = pPnts[0].Translated(v.Multiplied(0.5))
                    p2 = ProjectPointOnBSplineLst(b_BSplineLst, cP, 1)
                    newPnt = p2[0]
                    newPara = [LayerID, p2[1], p2[2]]
                    b_nodes.append(Node(newPnt, newPara))

                elif node.regular_corner == False:
                    print("WARNING: cornerstyle 0: this possibility has not been implemented yet.")
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))

                else:
                    if pPnts[0].IsEqual(pPnts[1], 1e-7 * refL):
                        b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))
                    else:
                        print("ERROR: cornerstyle 0: this possibility has not been implemented yet. pIdx[0] != pIdx[0]. Create Bisector and intersect with bsplinelist!")

            # ===CORNERSTYLE 1======
            elif len(exterior_corners) == 1 and node.corner == False:
                node.cornerstyle = 1
                # print('Node ID: ',node.id,', Len(exterior_corners):', len(exterior_corners))
                b_nodes.append(Node(exterior_corners[0], exterior_corners_para[0]))

            # ===CORNERSTYLE 2======
            elif len(exterior_corners) == 0 and node.corner == True:
                node.cornerstyle = 2
                v = gp_Vec2d(pPnts[0], pPnts[1])
                mid_Pnt = pPnts[0].Translated(v.Multiplied(0.5))
                v2 = gp_Vec2d(node.Pnt2d, mid_Pnt)
                v2.Normalize()
                far_Pnt = node.Pnt2d.Translated(v2.Multiplied(10))

                straigh_bspline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([node.Pnt2d, far_Pnt])).Curve()
                (iPnts, iPnts_Pnt2d) = intersect_BSplineLst_with_BSpline(b_BSplineLst, straigh_bspline)

                newPara = iPnts[0]
                newPnt = iPnts_Pnt2d[0]

                if node.regular_corner:
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))
                    b_nodes.append(Node(newPnt, [LayerID, newPara[0], newPara[1]]))
                    b_nodes.append(Node(pPnts[1], [LayerID, pIdx[1], pPara[1]]))

                else:
                    b_nodes.append(Node(pPnts[1], [LayerID, pIdx[1], pPara[1]]))
                    b_nodes.append(Node(newPnt, [LayerID, newPara[0], newPara[1]]))
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))

            # ===CORNERSTYLE 3======
            elif len(exterior_corners) == 1 and node.corner == True:
                node.cornerstyle = 3
                # print 'node.cornerstyle = 3 @', node
                if node.regular_corner == True:
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))
                    b_nodes.append(Node(exterior_corners[0], [exterior_corners_para[0][0], exterior_corners_para[0][1], exterior_corners_para[0][2]]))
                    b_nodes[-1].corner = True
                    b_nodes.append(Node(pPnts[1], [LayerID, pIdx[1], pPara[1]]))

                else:
                    b_nodes.append(Node(pPnts[1], [LayerID, pIdx[1], pPara[1]]))
                    b_nodes.append(Node(exterior_corners[0], [exterior_corners_para[0][0], exterior_corners_para[0][1], exterior_corners_para[0][2]]))
                    b_nodes[-1].corner = True
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))

            # ===CORNERSTYLE 4======
            elif len(exterior_corners) == 2 and node.corner == True:
                node.cornerstyle = 4

                if node.regular_corner == True:
                    # print 'R',[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]],[exterior_corners_para[1][0],exterior_corners_para[1][1],exterior_corners_para[1][2]]
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))
                    b_nodes.append(Node(exterior_corners[0], [exterior_corners_para[0][0], exterior_corners_para[0][1], exterior_corners_para[0][2]]))
                    b_nodes[-1].corner = True
                    # Find Middle between the two exterior corners on b_BsplineLst
                    newPnt = gp_Pnt2d()
                    c_BSplineLst = b_BSplineLst[exterior_corners_para[0][1] + 1 : exterior_corners_para[1][1] + 1]
                    [tmp_idx, tmp_u] = find_BSplineLst_coordinate(c_BSplineLst, 0.5, 0, 1)
                    newIdx = exterior_corners_para[0][1] + 1 + tmp_idx
                    newPara = tmp_u
                    b_BSplineLst[newIdx].D0(newPara, newPnt)

                    b_nodes.append(Node(newPnt, [LayerID, newIdx, newPara]))
                    b_nodes.append(Node(exterior_corners[1], [exterior_corners_para[1][0], exterior_corners_para[1][1], exterior_corners_para[1][2]]))
                    b_nodes[-1].corner = True
                    b_nodes.append(Node(pPnts[1], [LayerID, pIdx[1], pPara[1]]))

                else:
                    b_nodes.append(Node(pPnts[1], [LayerID, pIdx[1], pPara[1]]))
                    b_nodes.append(Node(exterior_corners[0], [exterior_corners_para[0][0], exterior_corners_para[0][1], exterior_corners_para[0][2]]))
                    b_nodes[-1].corner = True

                    v = gp_Vec2d(exterior_corners[0], exterior_corners[1])
                    cP = exterior_corners[0].Translated(v.Multiplied(0.5))
                    p2 = ProjectPointOnBSplineLst(b_BSplineLst, cP, distance)
                    newPnt = p2[0]
                    newPara = [LayerID, p2[1], p2[2]]

                    b_nodes.append(Node(newPnt, newPara))
                    b_nodes.append(Node(exterior_corners[1], [exterior_corners_para[1][0], exterior_corners_para[1][1], exterior_corners_para[1][2]]))
                    b_nodes[-1].corner = True
                    b_nodes.append(Node(pPnts[0], [LayerID, pIdx[0], pPara[0]]))

            # ===CORNERSTYLE 5======
            elif len(exterior_corners) > 2 and node.corner == True:
                node.cornerstyle = 5
                print("WARNING: cornerstyle 5 has not been implemented yet.")

            #                display.DisplayShape(pPnts[0], color='BLUE')
            #                display.DisplayShape(pPnts[1], color='RED')
            # TODO:
            # for p in exterior_corners:
            # display.DisplayShape(p,color='RED')
            # b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
            # b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
            # b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))

            # ===CORNERSTYLE 6======
            elif len(exterior_corners) == 4 and node.corner == True:
                node.cornerstyle = 6
                print("WARNING: cornerstyle 6 has not been implemented yet.")
                # TODO:
                # for p in exterior_corners:
                # display.DisplayShape(p,color='RED')
                # b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                # b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                # b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))

        else:
            print("Projection Error, number of projection points: ", len(pPnts))

    # ==============integrate_leftover_interior_nodes=========================================
    # determin_leftover_pnts
    if flag_integrate_leftover_interior_nodes:
        new_b_node = None
        new_a_node = None
        aglTol = 5.0
        linTol = 1e-9 * refL
        prjTol = 1e-5 * refL
        new_b_node = None
        insert_idx = None
        leftover_exterior_corners = []
        for i, item in enumerate(b_BSplineLst[:-1]):
            spline1 = item
            spline2 = b_BSplineLst[i + 1]
            u1, p1, v1 = spline1.LastParameter(), gp_Pnt2d(), gp_Vec2d()
            u2, p2, v2 = spline2.FirstParameter(), gp_Pnt2d(), gp_Vec2d()
            spline1.D1(u1, p1, v1)
            spline2.D1(u2, p2, v2)

            Angle = abs(v1.Angle(v2)) * 180 / np.pi
            if Angle > aglTol and not any(n.Pnt2d.IsEqual(item.EndPoint(), linTol) for n in b_nodes):
                leftover_exterior_corners.append((item.EndPoint(), [LayerID, i, u1]))

        # reversed projection and insert them into a_nodes and b_nodes
        for p1 in leftover_exterior_corners:
            p2 = ProjectPointOnBSplineLst(a_BSplineLst, p1[0], (1 + prjTol) * layer_thickness)
            if len(p2) > 0:
                # print 'WARNING: integrating leftover interior node.'
                new_b_node = Node(p1[0], p1[1])
                # display.DisplayShape(new_b_node.Pnt2d, color='RED')
                new_a_node = Node(p2[0], [LayerID, p2[1], p2[2]])
                # display.DisplayShape(new_a_node.Pnt2d, color='GREEN')

                # print 'new_a_node.parameters:',new_a_node.parameters,'a_nodes[0].parameters', a_nodes[0].parameters
                for i, n in enumerate(a_nodes):
                    if n.parameters[0] == new_a_node.parameters[0] and n.parameters[1] == new_a_node.parameters[1] and new_a_node.parameters[2] >= n.parameters[2]:
                        insert_idx = i
                        break

                if insert_idx != None:
                    b_nodes.append(new_b_node)
                    print("hello", new_b_node)
                    a_nodes.insert(insert_idx + 1, new_a_node)

        if new_b_node:
            print(new_b_node)
            # display.DisplayShape(new_b_node.Pnt2d, color='RED')
            # print(b_nodes[0].parameters[1], b_nodes[0].parameters[2])
            # b_nodes =  sorted(b_nodes, key=lambda n: (n.parameters[1], n.parameters[2]) )
            b_nodes = sorted(b_nodes)

    # ==============CREATE CELLS PROJECTION=========================================

    b = 0  # b_nodes idx
    if closed_a == True:
        start = 0
        end = len(a_nodes)
    else:
        start = 1
        end = len(a_nodes) - 1

    # for a,node in enumerate(a_nodes[1:-1], start=beginning):
    for a in range(start, end):
        try:
            # print 'Closed_a: ', closed_a, ', a: ', a, ', len(a_nodes): ', len(a_nodes),', b: ', b, ', len(b_nodes):', len(b_nodes), '\n',
            if closed_a == False and a == 1:  # Start Triangle
                cellLst.append(Cell([a_nodes[a], a_nodes[a - 1], b_nodes[b]]))

            elif closed_a == False and a == len(a_nodes) - 2:  # End Triangle
                cellLst.append(Cell([a_nodes[a - 1], b_nodes[b - 1], b_nodes[b], a_nodes[a]]))
                cellLst.append(Cell([a_nodes[a], b_nodes[b], a_nodes[a + 1]]))
                # print cellLst[-1]

            else:  # Regular Cell Creation
                if a_nodes[a].cornerstyle == 2 or a_nodes[a].cornerstyle == 3:
                    # print a, a_nodes[a], a_nodes[a].cornerstyle
                    cellLst.append(Cell([a_nodes[a - 1], b_nodes[b - 1], b_nodes[b], a_nodes[a]]))
                    b += 2
                    cellLst.append(Cell([a_nodes[a], b_nodes[b - 2], b_nodes[b - 1], b_nodes[b]]))

                elif a_nodes[a].cornerstyle == 4:
                    cellLst.append(Cell([a_nodes[a - 1], b_nodes[b - 1], b_nodes[b], a_nodes[a]]))
                    b += 2
                    cellLst.append(Cell([a_nodes[a], b_nodes[b - 2], b_nodes[b - 1], b_nodes[b]]))
                    b += 2
                    cellLst.append(Cell([a_nodes[a], b_nodes[b - 2], b_nodes[b - 1], b_nodes[b]]))

                else:
                    cellLst.append(Cell([a_nodes[a - 1], b_nodes[b - 1], b_nodes[b], a_nodes[a]]))

        except (IndexError):
            print("ERROR:\t IndexError: list index out of range", a_nodes[a])
            pass

        b += 1

    # ==============OCC3DVIEWER========================================
    if kw.get("display") != None:

        flag_display_a_nodes = True
        flag_display_b_nodes = True
        flag_display_a_BSplineLst = True
        flag_display_b_BSplineLst = True
        flag_display_cells = False

        if flag_display_a_nodes:
            for i, a in enumerate(a_nodes):
                if a.corner == True:
                    display.DisplayShape(a.Pnt, color="WHITE")
                    string = str(a.id) + " (cs=" + str(a.cornerstyle) + ", rg=" + str(a.regular_corner) + ")"
                    # display.DisplayMessage(a.Pnt,string,message_color=(1.0,0.0,0.0))

                elif a.cornerstyle == 1 or a.cornerstyle == 0:
                    display.DisplayShape(a.Pnt, color="WHITE")
                    string = str(a.id) + " (cs=" + str(a.cornerstyle) + ", rg=" + str(a.regular_corner) + ")"
                    # display.DisplayMessage(a.Pnt,string,message_color=(1.0,0.5,0.0))

                else:
                    display.DisplayShape(a.Pnt, color="WHITE")
                    # display.DisplayMessage(a.Pnt,str(a.id))

        if flag_display_b_nodes:
            for i, b in enumerate(b_nodes):
                display.DisplayShape(b.Pnt, color="GREEN")
                # display.DisplayMessage(b.Pnt,str(b.id),message_color=(1.0,0.5,0.0))

        if flag_display_a_BSplineLst:
            for i, a_spline in enumerate(a_BSplineLst):
                # display_custome_shape(display,a_spline,1.0,0.0,[0.2,0.9,0.8])
                display.DisplayShape(a_spline, color="CYAN")
                p = gp_Pnt2d()
                v = gp_Vec2d()
                u = (a_spline.LastParameter() - a_spline.FirstParameter()) / 2 + a_spline.FirstParameter()
                a_spline.D1(u, p, v)
                # display.DisplayMessage(p,str(i),height=30,message_color=(0,1,1))
                # display.DisplayVector(gp_Vec(v.X(),v.Y(),0), gp_Pnt(p.X(),p.Y(),0))

        if flag_display_b_BSplineLst:
            for i, b_spline in enumerate(b_BSplineLst):
                # display_custome_shape(display,b_spline,1.0,0.0,[0.1,0.5,1.0 ])
                display.DisplayShape(b_spline, color="BLUE")
                p = gp_Pnt2d()
                v = gp_Vec2d()
                u = (b_spline.LastParameter() - b_spline.FirstParameter()) / 2 + b_spline.FirstParameter()
                b_spline.D1(u, p, v)
                # display.DisplayMessage(p,str(i),height=30,message_color=(0,0,1))
                # display.DisplayVector(gp_Vec(v.X(),v.Y(),0), gp_Pnt(p.X(),p.Y(),0))

        if flag_display_cells:
            for c in cellLst:
                display.DisplayShape(c.wire, color="BLACK")

        display.View_Top()
        display.FitAll()

    return a_nodes, b_nodes, cellLst
