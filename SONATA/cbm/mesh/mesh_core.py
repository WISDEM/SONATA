# -*- coding: utf-8 -*-
"""
Functionalities concerning the unstructured discretization (triangulation) 
of the core or balance weight. 

Created on Thu Nov 02 10:46:29 2017
@author: TPflumm
"""

# Core Library modules
import pickle

# Third party modules
import matplotlib.pyplot as plt
import numpy as np
from OCC.Core.gp import gp_Pnt2d
from triangle import plot, triangulate

# First party modules
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.mesh_utils import equidistant_nodes_on_BSplineLst
from SONATA.cbm.mesh.node import Node
from SONATA.cbm.topo.BSplineLst_utils import get_BSplineLst_length

# ===================FUNCTIONALITIES==========================================


def triangle_mesh(array, options):
    """The triangle_mesh function generates the triagular mesh within the 
    a_nodes polygon. It uses the Python Triangle module, which is a python 
    wrapper around Jonathan Richard Shewchuks two-dimensional quality mesh
    generator and delaunay triangulator library, available here.  
    
    - http://dzhelil.info/triangle/
    - http://www.cs.cmu.edu/~quake/triangle.html
        
    Args:
        array: (array of nodes), the array contains all nodes that are on the
            innermost boundary of the generated topology. And are the boundary 
            for the triangulation
            
        options: the options can be passed to the triagle_mesh function,
            default is the optionstring 'pa%s' % (area)
            
   Returns: 
        mesh: a dictionary with 'vertices', 'segments', 'holes' and 'regions'
            as keys
    """

    tmp = []
    for i, v in enumerate(array):
        if i == len(array) - 1:
            tmp.append([i, 0])
        else:
            tmp.append([i, i + 1])
    segments_core = np.asarray(tmp)

    poly = {"vertices": None, "segments": None}
    poly["vertices"] = array
    poly["segments"] = segments_core
    # plt.plot(old_vertices[:,0],old_vertices[:,1],'.-')
    mesh = triangulate(poly, options)
    # plot.compare(plt, poly, mesh)
    return mesh

def gen_core_cells(a_nodes, area=1.0, **kwargs):
    """The gen_core_cells function generates the triagular mesh within the 
    a_nodes polygon.        
    
    Args:
        a_nodes: (list of nodes), the list contains all nodes that are on the
            innermost boundary of the generated topology.
        area: (float), is a resolution parameter / area constraint for the 
            triangle_mesh function.
    
    Kwargs: 
        options: the options can be passed to the triagle_mesh function,
            default is the optionstring 'pa%s' % (area)
            
        
    Returns: 
        [c_cells,c_nodes]: c_cells is a list of the newly generated cell 
            objects, while c_nodes is a list of the newly generated node 
            objects.              
    """

    # KWARGS:
    if kwargs.get("options") != None:
        options = kwargs.get("options")
    else:
        if area < 1.0:
            # options = 'pq'
            scalefactor = np.sqrt(1 / area) + 0.1
            area = area * scalefactor ** 2
            # print('area:', area, 'scalefactor:', scalefactor)
            options = "pa%f" % (area)  # Somehow crashing!?

        else:
            scalefactor = 1
            options = "pa%f" % (area)  # Somehow crashing!?

    old_vertices = np.asarray([[n.Pnt2d.X(), n.Pnt2d.Y(), n.id] for n in a_nodes])
    mesh = triangle_mesh(old_vertices[:, :2] * scalefactor, options)

    c_nodes = []
    connector = []  # connector stores the information [tri_vertex_id, old_node_id]
    for tri_id, v in enumerate(mesh["vertices"]):
        # print(type(v))
        v = v / scalefactor
        vec = np.linalg.norm(old_vertices[:, :2] - v, axis=1)
        idx = np.argmin(vec)
        value = vec[idx]
        if value <= 1e-5:
            # checks if the vertex exists allready in the a_nodes list
            connector.append([tri_id, old_vertices[idx, 2]])
        else:
            # create new node for every new vertice that is not in a_nodes:
            tmp_node = Node(gp_Pnt2d(v[0], v[1]))
            connector.append([tri_id, tmp_node.id])
            c_nodes.append(tmp_node)

    connector = np.asarray(connector)
    nodes = a_nodes + c_nodes

    # create cells from connector array
    c_cells = []
    for ele in mesh["triangles"]:
        nodeids = [connector[tri_id, 1] for tri_id in ele]
        nodedict = {n.id: n for n in nodes}
        nodeLst = [nodedict[i] for i in nodeids]
        c_cells.append(Cell(nodeLst))

    return [c_cells, nodes]


# ===================MAIN==========================================
if __name__ == "__main__":
    # LOAD .pkl data with SegmentLst
    filename = "C:\TPflumm_local\work\SONATA\sec_config.pkl"
    with open(filename, "rb") as handle:
        SegmentLst = pickle.load(handle)

    # Build wires for each layer and segment
    for seg in SegmentLst:
        seg.build_wire()
        for layer in seg.LayerLst:
            layer.build_wire()

    layer = SegmentLst[-1].LayerLst[-1]
    Resolution = 100  # Nb of Points on Segment0
    length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
    global_minLen = round(length / Resolution, 5)

    a_nodes = equidistant_nodes_on_BSplineLst(layer.Boundary_BSplineLst, True, True, True, minLen=global_minLen, LayerID=layer.ID[0])
    [c_cells, c_nodes] = gen_core_cells(a_nodes, 2)

    for c in c_cells:
        c.structured = False
        c.theta_3 = 0
        c.MatID = int(seg.CoreMaterial)
        c.calc_theta_1()

