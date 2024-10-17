# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:01:06 2017

@author: TPflumm
"""
# Core Library modules
import os
import math
# Third party modules
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
import matplotlib.cm as cm


# First party modules
# from SONATA.utl_openmdao.doe_utl import filename_generator


# mpl_logger = logging.getLogger("matplotlib")
# mpl_logger.setLevel(logging.WARNING)

def centroid(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    centroid = (sum(x) / len(points), sum(y) / len(points))
    return centroid

def plot_mesh(nodes, elements, theta_11, data, data_name, materials, title=None, VABSProperties=None, 
              show_element_number=False, show_node_number=False, invert_xaxis = True, lfactor=0.5e-2, **kw):
    
    """
    To be continued...
    
    Parameters
    ---------- 
    nodes : list
    elements : list
    data : 
    data_name : string
    
    """
    alpha = 1.
    if "cmap" in kw:
        cmap = plt.cm.get_cmap(kw["cmap"])

    elif data_name == "sf":
        colors = [(0.6, 0, 0), (1, 1, 0), (0, 0.5, 0)]  # R -> G -> B
        cmap_name = 'my_list'
        cmap = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=6)
        
    elif data_name == 'MatID':
        cmap = plt.cm.get_cmap()
        # extract all colors from the .jet map
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, max(data))
        
        cmap = cm.get_cmap('tab20',max(data))
    else:
        cmap = plt.cm.get_cmap()

    if "vmin" in kw:
        vmin = kw["vmin"]
    else:
        vmin = None

    if "vmax" in kw:
        vmax = kw["vmax"]
    else:
        vmax = None

    fig, ax = plt.subplots(1,1,figsize=(7.5, 6))
    patches = []
    centroids = []
    for i, ele in enumerate(elements):
        # print ele
        if int(0) in ele:
            array = np.vstack((nodes[ele[0] - 1], nodes[ele[1] - 1], nodes[ele[2] - 1]))
            centroids.append(centroid(array))
        else:
            array = np.vstack((nodes[ele[0] - 1], nodes[ele[1] - 1], nodes[ele[2] - 1], nodes[ele[3] - 1]))
            centroids.append(centroid(array))
        polygon = Polygon(array, closed=True, edgecolor="k")
        patches.append(polygon)
    
    p = PatchCollection(patches, alpha=alpha, cmap=cmap, edgecolors = 'k', linewidths=0.2)
    p.set_array(data)
    p.set_clim(vmin, vmax)
    _ = ax.add_collection(p)

    cbar = fig.colorbar(p, ax=ax, drawedges=(data_name == 'MatID'))
    cbar.ax.set_ylabel(data_name)

    if data_name == "MatID":
        cbar.set_ticks(np.linspace(1, max(data), max(data)))
        cbar_label = []
        for cbar_label_index in range(max(data)):
            cbar_label.append(materials[cbar_label_index+1].name)
        p.set_clim(0.5, max(data)+0.5)
        cbar.set_ticklabels(cbar_label)

    elif data_name[0:6] == 'stress':
        cbar.ax.set_ylabel('Stress, $\mathrm{Nm^2}$')

    elif data_name[0:6] == 'strain':
        cbar.ax.set_ylabel('Strain, m/m')

    if len(theta_11)==len(elements):
        for i,cent in enumerate(centroids):
            
            dx = lfactor*math.cos(math.radians(theta_11[i]))
            dy = lfactor*math.sin(math.radians(theta_11[i]))
            ax.arrow(cent[0], cent[1], dx, dy, width = 0.01e-2, head_width=0.1e-2, head_length=0.1e-2, fc='k', ec='k')
        
    plt.axis('equal')


    if title!=None:
        ax.set_title(title, fontweight = 'bold')
    ax.set_xlabel('x, m', fontweight = 'bold')
    ax.set_ylabel('y, m', fontweight = 'bold')
    # plt.grid(color=[0.8,0.8,0.8], linestyle='--')

    ##display element number
    if show_element_number == True:
        for i, item in enumerate(centroids):
            ax.annotate(i + 1, (item[0], item[1]))

    ##display node number
    if show_node_number == True:
        for i, item in enumerate(nodes):
            ax.annotate(i + 1, (item[0], item[1]), color="red")

    # if VABSProperties != None:
    #     pass
    #     (CG,) = plt.plot(VABSProperties.Xm[0], VABSProperties.Xm[1], "ro", label="CG: Mass Center")
    #     # ax.annotate('CG', (VABSProperties.Xm2,VABSProperties.Xm3),fontsize=20)
    #     (NA,) = plt.plot(VABSProperties.Xt[0], VABSProperties.Xt[1], "gs", label="NA: Neutral Axes")
    #     # ax.annotate('NA', (VABSProperties.Xt2,VABSProperties.Xt3),fontsize=20)
    #     try:
    #         (GC,) = plt.plot(VABSProperties.Xg[0], VABSProperties.Xg[1], "b^", label="GC: Geometric Center")
    #         # ax.annotate('GC', (VABSProperties.Xg2,VABSProperties.Xg3),fontsize=20)
    #         plt.legend(handles=[CG, GC, NA])
    #     except:
    #         plt.legend(handles=[CG, NA])

    #     if isinstance(VABSProperties.Xs, np.ndarray):
    #         (SC,) = plt.plot(VABSProperties.Xs[0], VABSProperties.Xs[1], "kD", label="SC: Generalized Shear Center")
    #         # ax.annotate('SC', (VABSProperties.Xs2,VABSProperties.Xs3),fontsize=20)
    #         try:
    #             plt.legend(handles=[CG, GC, NA, SC])
    #         except:
    #             plt.legend(handles=[CG, NA, SC])

    if invert_xaxis:
        ax.invert_xaxis()

    plt.show()
    return (fig, ax)


    return (fig,ax)
    
def plot_cells(cells,nodes, attr1, materials, VABSProperties=None, title='None', plotTheta11=False, plotDisplacement=False, **kw):
    """
    

    Parameters
    ----------
    cells : TYPE
        DESCRIPTION.
    nodes : TYPE
        DESCRIPTION.
    attr1 : TYPE
        DESCRIPTION.
    materials : TYPE
        DESCRIPTION.
    VABSProperties : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is 'None'.
    plotTheta11 : TYPE, optional
        DESCRIPTION. The default is False.
    plotDisplacement : TYPE, optional
        DESCRIPTION. The default is False.
    **kw : TYPE
        DESCRIPTION.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.

    """
    nodes_array = []
    for n in nodes:
        if plotDisplacement and n.displacement[0] is not None:
            nodes_array.append([[n.coordinates[0]+n.displacement[1]-1,n.coordinates[1]+n.displacement[2]-1]])
        else:
            nodes_array.append([n.coordinates[0], n.coordinates[1]])
    nodes_array = np.asarray(nodes_array)

    element_array = []
    for c in cells:
        tmp = []
        for i in range(0, 4):
            if i < len(c.nodes):
                tmp.append(c.nodes[i].id)
            else:
                tmp.append(0)
        element_array.append(tmp)
    element_array = np.asarray(element_array)

    data = []
    for c in cells:
        if "." in attr1:
            attrLst = attr1.split(".")
            tmp = getattr(c, attrLst[0])
            data.append(getattr(tmp, attrLst[1]))
        else:
            data.append(getattr(c, attr1))
    data = np.asarray(data)
    data_name = attr1

    theta_11 = []
    if plotTheta11 == True:
        for c in cells:
            theta_11.append(getattr(c, "theta_11"))
        theta_11 = np.asarray(theta_11)
    
    
    fig,ax = plot_mesh(nodes_array, element_array, theta_11, data, data_name, materials, title, VABSProperties, **kw)    
   
    if 'savepath' in kw:

        if not os.path.exists(os.path.join(kw['savepath'],'figures')):  # create 'figures' Folder if not already existing
            os.makedirs(os.path.join(kw['savepath'],'figures'))

        # datestr = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        # fname = kw['savepath'].split('.')[0]+'_'+datestr+'.'+kw['savepath'].split('.')[1]

        if 'opt_var' in kw:
            fname = kw['savepath'] + '/figures/blade_section_' + kw['section'] + '_optvar_' + kw['opt_var'] + '.pdf'
        else:
            fname = kw['savepath']+'/figures/blade_section_'+kw['section']+'.pdf'
        # print(fname)
        # tmp_fig = plt.gcf()
        tmp_fig = fig
        # tmp_fig.set_size_inches(11.69, 8.27)    #a4 landscape
        # tmp_fig.set_size_inches(10, 5)    #a4 landscape
        tmp_fig.savefig(fname, dpi=300, orientation='landscape')

    return (fig, ax)


# =============================================================
#       IF MAIN:
# ============================================================
if __name__ == "__main__":

    nodes = scipy.io.loadmat("nodes.mat")
    nodes = nodes["nodes"]
    elements = scipy.io.loadmat("elements.mat")
    elements = elements["elements"]
    layup_angle = scipy.io.loadmat("layupAngle.mat")
    layup_angle = layup_angle["layup_angle"]
    layup_angle = layup_angle[:, 0]

    plot_mesh(nodes, elements, layup_angle, False, False)
