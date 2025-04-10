# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 10:30:01 2017

@author: TPflumm
"""

# Third party modules
import numpy as np
from OCC.Core.Geom2dAPI import (Geom2dAPI_ProjectPointOnCurve,)
from OCC.Core.gp import gp_Vec2d

# First party modules
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.mesh_utils import (theta_1_from_2nodes,)
from SONATA.cbm.mesh.node import Node
from SONATA.cbm.topo.BSplineLst_utils import ProjectPointOnBSplineLst


def modify_sharp_corners(cells, b_BSplineLst, global_minLen, layer_thickness, LayerID=0, tol=1e-2, alpha_crit=50, **kwargs):
    """
    
    """

    # KWARGS:
    if kwargs.get("display") != None:
        display = kwargs.get("display")

    enhanced_cells = []
    for i, c in enumerate(cells):
        if len(c.nodes) == 4:

            if c.nodes[0].cornerstyle == 2 or c.nodes[0].cornerstyle == 3:

                v21 = gp_Vec2d(c.nodes[2].Pnt2d, c.nodes[1].Pnt2d)
                v23 = gp_Vec2d(c.nodes[2].Pnt2d, c.nodes[3].Pnt2d)
                angle = abs(v21.Angle(v23) * 180 / np.pi)

                # print c,angle
                if angle < alpha_crit:
                    L = c.nodes[0].Pnt2d.Distance(c.nodes[2].Pnt2d) * 1.5
                    BS_Vec2d = gp_Vec2d(c.nodes[0].Pnt2d, c.nodes[2].Pnt2d)
                    MiddleNodes = []
                    if int(L // global_minLen)-1 >= 2:
                        for i in range(0, int(L // global_minLen) - 1):
                            P = c.nodes[0].Pnt2d.Translated(BS_Vec2d.Multiplied((1 + i) / float(int(L // global_minLen))))
                            MiddleNodes.append(Node(P))

                    FrontNodes = []
                    BackNodes = []
                    distance = (1 + tol) * layer_thickness

                    # Node indices of Middle nodes that fail to project and
                    # identify front/back.
                    rm_MiddleNodes = []

                    for mid_idx, n in enumerate(MiddleNodes):
                        pPnts = []
                        pPara = []
                        pIdx = []
                        for idx, item in enumerate(b_BSplineLst):
                            projection = Geom2dAPI_ProjectPointOnCurve(n.Pnt2d, item)
                            for j in range(1, projection.NbPoints() + 1):
                                if projection.Distance(j) <= distance:
                                    if not any(item.IsEqual(projection.Point(j), 1e-6) for item in pPnts):
                                        pPnts.append(projection.Point(j))
                                        pPara.append(projection.Parameter(j))
                                        pIdx.append(idx)
                                else:
                                    None

                        # print 'Nuber of Projected Middle Nodes pPnts:', len(pPnts)
                        trigger_f = True
                        trigger_b = True
                        for i, P in enumerate(pPnts):
                            v01 = gp_Vec2d(c.nodes[0].Pnt2d, c.nodes[1].Pnt2d)
                            vnP = gp_Vec2d(n.Pnt2d, P)

                            if len(pPnts) > 2:
                                print("vnP.Dot(v01): ", vnP.Dot(v01))

                            if vnP.Dot(v01) > 0 and trigger_f:
                                trigger_f = False
                                FrontNodes.append(Node(P, [LayerID, pIdx[i], pPara[i]]))

                            elif vnP.Dot(v01) < 0 and trigger_b:
                                trigger_b = False
                                BackNodes.append(Node(P, [LayerID, pIdx[i], pPara[i]]))

                            else:
                                # This need not be a complete error if other
                                # Middle nodes have worked.
                                print("WARNING: cannot determine FRONT and ",
                                      "BACK nodes @ ",
                                      c.nodes[0],
                                      ", so not using this point for sharp ",
                                      "corner improvement.")

                                rm_MiddleNodes.append(mid_idx)

                        # Make sure 'FrontNodes' and 'BackNodes' are the same
                        # length and fix.
                        # This fix is needed in case the front or back projection
                        # worked, but the other did not. In that case,
                        # rm_MiddleNodes should have this index.
                        if mid_idx in rm_MiddleNodes:

                            num_proj = min(len(FrontNodes), len(BackNodes))

                            if len(FrontNodes) - num_proj > 2\
                                or len(BackNodes) - num_proj > 2:
                                # Checking every for loop iteration, so these
                                # lists should not get off by a length of more
                                # than 1.
                                print("ERROR: Fixing shaper elements failed @ ", c.nodes[0])
                                MiddleNodes = []
                                FrontNodes = []
                                BackNodes = []
                            else:
                                FrontNodes = FrontNodes[:num_proj]
                                BackNodes = BackNodes[:num_proj]

                    # Remove MiddleNodes that have indices in 'rm_MiddleNodes'
                    MiddleNodes = [node for k, node in enumerate(MiddleNodes)
                                   if k not in rm_MiddleNodes]

                    if len(MiddleNodes) == 0 and len(rm_MiddleNodes) > 0:
                        print("ERROR: All refinements failed for sharp corner @ ",
                              c.nodes[0], "because vnp and v01 are orthogonal")

                    # =====================CREATE FRONT CELLS
                    FrontCellLst = []
                    # print '@', c.nodes[0],'  len(Middle):',len(MiddleNodes),'len(Front):',len(FrontNodes),'len(Back):',len(BackNodes)
                    if len(MiddleNodes) == len(FrontNodes) and len(BackNodes) == len(MiddleNodes) and MiddleNodes:
                            
                        for i in range(0, len(MiddleNodes)):

                            if i == 0:  # FIRST
                                nodeLst = [c.nodes[0], c.nodes[1], FrontNodes[i], MiddleNodes[i]]
                            else:
                                nodeLst = [MiddleNodes[i - 1], FrontNodes[i - 1], FrontNodes[i], MiddleNodes[i]]
                            FrontCellLst.append(Cell(nodeLst))

                        if len(MiddleNodes) > 0:  # LAST
                            nodeLst = [MiddleNodes[-1], FrontNodes[-1], c.nodes[2]]
                            FrontCellLst.append(Cell(nodeLst))

                        # =====================CREATE BACK CELLS
                        BackCellLst = []
                        for i in range(0, len(MiddleNodes)):

                            if i == 0:  # FIRST
                                nodeLst = [MiddleNodes[i], BackNodes[i], c.nodes[3], c.nodes[0]]
                            else:
                                nodeLst = [MiddleNodes[i], BackNodes[i], BackNodes[i - 1], MiddleNodes[i - 1]]
                            BackCellLst.append(Cell(nodeLst))

                        if len(MiddleNodes) > 0:  # LAST
                            nodeLst = [MiddleNodes[-1], c.nodes[2], BackNodes[-1]]
                            BackCellLst.append(Cell(nodeLst))

                        enhanced_cells.extend(FrontCellLst)
                        enhanced_cells.extend(reversed(BackCellLst))

                    if len(MiddleNodes) == 0:
                        # No refinement found, keep old cell.
                        enhanced_cells.append(c)
                else:
                    enhanced_cells.append(c)

            else:
                enhanced_cells.append(c)
                # print('modify_sharp_corners has not been implemented for cornerstyle ', c.nodes[0].cornerstyle)
        else:
            enhanced_cells.append(c)
    try:
        new_b_nodes = FrontNodes + BackNodes
    except:
        new_b_nodes = []

    return enhanced_cells, new_b_nodes

def second_stage_improvements(cells, b_BSplineLst, global_minLen, LayerID=0, factor1=1.8, factor2=0.15, **kw):

    if kw.get("display") != None:
        display = kw.get("display")

    enhanced_cells2 = []
    new_b_nodes = []
    for i, c in enumerate(cells):
        if len(c.nodes) == 4:
            v = gp_Vec2d(c.nodes[1].Pnt2d, c.nodes[2].Pnt2d)
            magnitude = v.Magnitude()

            # SPLIT CELLS INTO TRIANGLES AND ADD NODE!
            if magnitude >= factor1 * global_minLen:
                cP = c.nodes[1].Pnt2d.Translated(v.Multiplied(0.5))
                # display.DisplayColoredShape(cP, 'GREEN')
                p2 = ProjectPointOnBSplineLst(b_BSplineLst, cP, 1)
                # display.DisplayColoredShape(p2[0], 'YELLOW')
                nodeLst = c.nodes
                newNode = Node(p2[0], [LayerID, p2[1], p2[2]])
                # MODIFY EXISTING CELL
                c.nodes = [nodeLst[0], nodeLst[1], newNode]
                enhanced_cells2.append(c)
                enhanced_cells2[-1].calc_theta_1()
                # ADD NEW CELLS
                enhanced_cells2.append(Cell([nodeLst[0], newNode, nodeLst[3]]))
                enhanced_cells2[-1].theta_1 = theta_1_from_2nodes(nodeLst[0], nodeLst[3])
                # Append last triangle
                enhanced_cells2.append(Cell([nodeLst[3], newNode, nodeLst[2]]))
                enhanced_cells2[-1].calc_theta_1()

            # MERGE NODES when to small
            elif magnitude <= factor2 * global_minLen and c.nodes[1].corner == False and c.nodes[2].corner == False:
                cP = c.nodes[1].Pnt2d.Translated(v.Multiplied(0.5))
                # display.DisplayColoredShape(cP, 'GREEN')
                p2 = ProjectPointOnBSplineLst(b_BSplineLst, cP, 1)
                # display.DisplayColoredShape(p2[0], 'RED')
                nodeLst = c.nodes
                # Modify Node 2
                nodeLst[2].Pnt2d = p2[0]
                nodeLst[2].parameters = [LayerID, p2[1], p2[2]]
                # MODIFY EXISTING CELL
                c.nodes = [nodeLst[0], nodeLst[2], nodeLst[3]]
                c.theta_1 = theta_1_from_2nodes(nodeLst[0], nodeLst[3])
                enhanced_cells2.append(c)

                # MODIFY Last CELL
                cells[i - 1].nodes[2] = nodeLst[2]

            else:
                enhanced_cells2.append(c)

        else:
            enhanced_cells2.append(c)

        try:
            new_b_nodes.append(newNode)
        except:
            None

    return enhanced_cells2, new_b_nodes
