# -*- coding: utf-8 -*-
"""
Functionalities concerning the consolidation between the individual meshes in 
the different segments. Segments>0 are seperated by a web and meshed individually. 
The web is defining the interface between the segemnts. The functions in this 
module make sure that the nodes of both segments are matched and hanging nodes 
are eliminated by splitting cells.

Created on Thu Nov 02 10:46:29 2017
@author: TPflumm
"""

# Third party modules
import numpy as np
from OCC.Core.gp import gp_Vec2d

# First party modules
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.mesh_utils import find_cells_that_contain_node


def consolidate_mesh_on_web(web, w_tol, display=None):
    """ Consolidates mesh on the web interface.
    
    After the mesh has been generated for every segment. This function makes
    sure that no hanging nodes remain in the mesh.
    
    Args:
        mesh: The overall mesh (list of cells) to consolidate and also to return
        w_BSplineLst: is the definition of the web as a BSplineLst, usually
                    it is a List with one spline that is a straight line.
        w1_nodes: are the nodes of the left segment that lie on the web
        w2_nodes: are the nodes of the right segment that lie on the web 
        w_tol (float) : tolerance to decide whether to match a left node to a right node 
               or if the cell should rather be devided into triangles.
        display: 
            
    Returns: 
        mesh: the updated mesh formulation (list of cells)
                
    """
    w1_nodes = web.wl_nodes
    w2_nodes = web.wr_nodes
    mesh = web.wl_cells + web.wr_cells

    # Generate NodeMatching Matrix [NM]
    tmp = [[w1.Pnt2d.Distance(w2.Pnt2d), w1.id, w2.id] for w1 in w1_nodes for w2 in w2_nodes]
    tmp = np.asarray(tmp)

    NM = tmp[tmp[:, 0] < w_tol]
    NM = NM[NM[:, 0].argsort()]

    # remove possible double nodes in NodeMatching Matrix
    tmp, tmp_idx = np.unique(NM[:, 1], True, axis=0)
    tmp_idx.sort(axis=0)
    NM = NM[tmp_idx]
    tmp, tmp_idx = np.unique(NM[:, 2], True, axis=0)
    tmp_idx.sort(axis=0)
    NM = NM[tmp_idx]

    # MERGE Nodes according to NodeMatching(NM) Matrix:
    for match in NM:
        # print int(match[1]),int(match[2])
        n1 = [x for x in w1_nodes if x.id == match[1]][0]
        n2 = [x for x in w2_nodes if x.id == match[2]][0]
        # move the first node to the middle!
        n1.Pnt2d
        v = gp_Vec2d(n1.Pnt2d, n2.Pnt2d)
        v.Multiply(0.5)
        n1.Pnt2d.Translate(v)

        # find all cells that contain n2 and replace it with n1
        disco_cells = find_cells_that_contain_node(mesh, n2)
        for c in disco_cells:
            c.nodes = [n1 if x == n2 else x for x in c.nodes]

    # determine ramaining w1 and w2_nodes. If w1_node.id is not in NM[:,1]
    rem_w1_nodes = [w1 for w1 in w1_nodes if w1.id not in NM[:, 1]]
    rem_w2_nodes = [w2 for w2 in w2_nodes if w2.id not in NM[:, 2]]

    #    if display !=  None:
    #        for n in rem_w1_nodes:
    #            display.DisplayShape(n.Pnt2d)
    #
    #        for n in rem_w2_nodes:
    #            display.DisplayShape(n.Pnt2d, color = 'RED')

    newcells = []
    newcells.extend(split_cells_to_consolidate(web.wl_cells, rem_w2_nodes, display))
    newcells.extend(split_cells_to_consolidate(web.wr_cells, rem_w1_nodes, display))
    return newcells


def split_cells_to_consolidate(cells, rem_nodes, display):
    """Subfunction to split cells with hanging nodes.
    
    Subfunction of consolidate_mesh_on_web to split cells with hanging nodes.
    It first identifies cells that are beeing intersected by the hanging nodes
    and then splits quad cells into 1 triangle and one quad and a triangle cell
    into two triangles.
    
    Args:
        mesh: The overall mesh (list of cells) to consolidate and also to return
        rem_nodes: remaining hanging nodes, where no match was found
        display: passes OCC display environment
            
    Returns: 
        mesh: the updated mesh formulation (list of cells)
                        
    """
    newcells = []
    for c in cells:
        tmp_nodes = []
        tmp_newcells = []
        for rn in rem_nodes:
            if c.cell_node_distance(rn) < 1e-4 and (rn not in c.nodes):
                tmp_nodes.append(rn)
                edg_idx = c.closest_cell_edge(rn)

        if len(c.nodes) == 4 and len(tmp_nodes) > 0:
            edg_idx = c.closest_cell_edge(tmp_nodes[0])
            nodeLst = c.nodes

            if display:
                display.DisplayShape(c.wire, color="RED")

                for it, n in enumerate(nodeLst):
                    display.DisplayShape(n.Pnt2d)
                    display.DisplayMessage(n.Pnt, str(it), message_color=(1.0, 0.0, 0.0))
                for n in tmp_nodes:
                    display.DisplayShape(n.Pnt2d, color="RED")

            # TODO: sort tmp nodes based on web.BSplineLst postion
            if len(tmp_nodes) == 1:
                tmp_newcells.append(Cell([nodeLst[edg_idx - 1], tmp_nodes[0], nodeLst[edg_idx - 2]]))
                c.nodes = [nodeLst[edg_idx - 1], nodeLst[edg_idx], tmp_nodes[0]]
                tmp_newcells.append(Cell([nodeLst[edg_idx - 2], tmp_nodes[0], nodeLst[edg_idx - 3]]))

            elif len(tmp_nodes) == 2:
                tmp_newcells.append(Cell([nodeLst[edg_idx - 1], nodeLst[edg_idx], tmp_nodes[1]]))  # ADD NEW CELL
                c.nodes = [nodeLst[edg_idx - 1], tmp_nodes[1], tmp_nodes[0], nodeLst[edg_idx - 2]]
                tmp_newcells.append(Cell([nodeLst[edg_idx - 2], tmp_nodes[0], nodeLst[edg_idx - 3]]))

            elif len(tmp_nodes) == 3:
                tmp_newcells.append(Cell([nodeLst[edg_idx - 1], nodeLst[edg_idx], tmp_nodes[2]]))  # ADD NEW CELL
                tmp_newcells.append(Cell([nodeLst[edg_idx - 1], tmp_nodes[2], tmp_nodes[1]]))  # ADD NEW CELL
                c.nodes = [nodeLst[edg_idx - 1], tmp_nodes[1], nodeLst[edg_idx - 2]]  # MODIFY EXISTING CELL#ADD NEW CELL
                tmp_newcells.append(Cell([nodeLst[edg_idx - 2], tmp_nodes[1], tmp_nodes[0]]))
                tmp_newcells.append(Cell([nodeLst[edg_idx - 2], tmp_nodes[0], nodeLst[edg_idx - 3]]))  # ADD NEW CELL
                # tmp_newcells.append(Cell([nodeLst[edg_idx-2],tmp_nodes[2],nodeLst[edg_idx-3]]))       #ADD NEW CELL

            elif len(tmp_nodes) == 4:
                tmp_newcells.append(Cell([nodeLst[edg_idx - 1], nodeLst[edg_idx], tmp_nodes[3]]))  # ADD NEW CELL
                tmp_newcells.append(Cell([nodeLst[edg_idx - 1], tmp_nodes[3], tmp_nodes[2]]))  # ADD NEW CELL
                c.nodes = [nodeLst[edg_idx - 1], tmp_nodes[2], tmp_nodes[1], nodeLst[edg_idx - 2]]  # MODIFY EXISTING CELL#ADD NEW CELL
                tmp_newcells.append(Cell([nodeLst[edg_idx - 2], tmp_nodes[1], tmp_nodes[0]]))
                tmp_newcells.append(Cell([nodeLst[edg_idx - 2], tmp_nodes[0], nodeLst[edg_idx - 3]]))  # ADD NEW CELL
                # tmp_newcells.append(Cell([nodeLst[edg_idx-2],tmp_nodes[2],nodeLst[edg_idx-3]]))       #ADD NEW CELL

        for tc in tmp_newcells:
            tc.theta_1 = c.theta_1
            tc.theta_3 = c.theta_3
            tc.MatID = c.MatID
        newcells.extend(tmp_newcells)

    return newcells
