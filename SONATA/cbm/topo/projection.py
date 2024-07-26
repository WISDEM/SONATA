# -*- coding: utf-8 -*-
"""
The SONATA.topo.projection module uses the python invervaltree package 
(Editable interval tree data structure for Python) to determine the cummulated 
layup boundaries of the SONATA Layup definition.

 consits of the following functions:
    
    - insert_interval_in_layup(layup,begin,end): insert_interval_in_layup 
            generates a intervaltree structure t1 from layup 
            and inserts a new intervaltree set t2 [begin, end] or 
            [0,end][begin,1] into the existing the set t1
    - chop_interval_from_layup(layup,begin,end): chops an intervaltree 
            structure t1 from layup to a given interval [begin, end] or 
            [0,end][begin,1].

Created on Mon Dec 19 16:32:15 2016
@author: TPflumm
"""
# Third party modules
import intervaltree
import matplotlib.pyplot as plt
import numpy as np

# =====================================================

def insert_interval_in_layup(layup, begin, end, **kw):
    """ insert_interval_in_layup generates a intervaltree structure t1 from layup 
    and inserts a new intervaltree set t2 [begin, end] or [0,end][begin,1] into
    the existing the set t1. First t1 is chopped with t2, then t2 is added to t1.
    
    for documentation of the intervaltree module refer to: 
        https://pypi.python.org/pypi/intervaltree
    
    Args:
        begin: float[0:1]: start of the interval 
        end: float[0:1]  end of the interval to be inserted. Note that begin
            can be larger than end (The layer is defined across the origin)
            
    Returns: 
        Projection: numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.7  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])                
    """

    # KWARGS:
    if kw.get("value") != None:
        idx = kw.get("value")
    else:
        idx = max(layup[:, 2]) + 1

    t1 = intervaltree.IntervalTree()
    for i, item in enumerate(layup):
        t1.addi(float(item[0]), float(item[1]), int(item[2]))

    # Initialize to be projected interval between begin and end
    t2 = intervaltree.IntervalTree()
    if begin > end:
        t2.addi(0, end)
        t2.addi(begin, 1)
    else:
        t2.addi(begin, end)

    # chop intervaltree t1
    for iv in t2:
        t1.chop(iv.begin, iv.end)

    # insert t2
    for iv in t2:
        t1.addi(iv.begin, iv.end, idx)

    # transform back to python np.array
    tmp_pj = []
    for iv in t1:
        tmp_pj.append([iv.begin, iv.end, iv.data])
    Projection = np.asarray(tmp_pj)
    return Projection


def chop_interval_from_layup(layup, begin, end):
    """ chop_interval_from_layup chops an intervaltree structure t1 from layup 
    to a given interval [begin, end] or [0,end][begin,1].
    
    for documentation of the intervaltree module refer to: 
        https://pypi.python.org/pypi/intervaltree
    
    Args:
        begin: float[0:1]: start of the interval 
        end: float[0:1]  end of the interval to which t1 is beeing chopped.
        Note that begin can be larger than end (The layer is defined across
        the origin)
            
    Returns: 
        Projection: numpy.array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.7  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])                
    """
    # get layup to intervaltree definition
    t1 = intervaltree.IntervalTree()
    for i, item in enumerate(layup):
        t1.addi(float(item[0]), float(item[1]), int(item[2]))

    # result_tree chop
    if begin > end:
        t1.chop(end, begin)
    else:
        t1.chop(0, begin)
        t1.chop(end, 1)

    # transform back to python np.array
    tmp_pj = []
    for iv in t1:
        tmp_pj.append([iv.begin, iv.end, iv.data])
    Projection = np.asarray(tmp_pj)
    return Projection


def sort_layup_projection(Projection):
    """The sort_layup_projection fuctions sorts the interval structure.
        
    Args:         
        Projection: list of numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.5  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])    
            
    Returns: 
        Projection: list of sorted numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]
                  [ 0.3  0.5  3. ]])               
    """
    # print(Projection)
    sorted_Projection = []
    for i, proj in enumerate(Projection):
        a = proj[:, :2].flatten()
        unique, counts = np.unique(a, return_counts=True)
        iv_boundaries = np.sort(unique[np.where(counts == 1)])
        # print(a, unique, counts, iv_boundaries)
        b = proj[proj[:, 0].argsort()]
        # regular interval:
        if np.all(a >= iv_boundaries[0]) and np.all(a <= iv_boundaries[1]):
            tmp = b
        # regular interval:
        else:
            splitter = (np.where(b[:, 0] == iv_boundaries[2])[0]).item()
            c = [b[0:splitter, :], b[splitter::, :]]
            # print i,"Split:",b, '@:', splitter, '\n C:', c
            tmp = np.vstack((c[1], c[0]))

        sorted_Projection.append(tmp)

    return sorted_Projection


def relevant_cummulated_layup_boundaries(Layup):
    """The relevant_cummulated_layup_boundaries function chops this valuable 
    information to the interval of the Layup itself to determin the offset
    origin in the given interval.        
    
    Args:        
        Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
                  #Start[-]	End[-] thickness [mm] Orientation [deg] MatID
            
        Returns: Projection: list of numpy arrays of the intervaltree structure 
            of the form:
                          start end layer#
                  np.array([[ 0.3  0.7  3. ]
                          [ 0.7  1.   2. ]
                          [ 0.   0.3  2. ]])                
    """
    if Layup.size == 0:
        return []

    else:
        nlayers = int(np.size(Layup, 0))
        # clean up Layup array
        a = np.delete(Layup, [2, 3, 4], 1)
        b = np.transpose(np.linspace(1, nlayers, nlayers))
        c = np.c_[a, b]
        clean_layup = np.insert(c, 0, np.array((0, 1, 0)), 0)
        tmp = clean_layup[[0]]
        # print clean_layup

        # Iterate over Layup
        relevant_projectionlist = []
        for i in range(1, nlayers + 1, 1):
            begin = float(clean_layup[i][0])
            end = float(clean_layup[i][1])
            relevant_projectionlist.append(chop_interval_from_layup(tmp, begin, end))
            tmp = insert_interval_in_layup(tmp, begin, end)
            tmp = tmp[np.lexsort(np.fliplr(tmp).T)]

        return sort_layup_projection(relevant_projectionlist)


def inverse_relevant_cummulated_layup_boundaries(Layup):
    """ The function inverse_relevant_cummulated_layup_boundaries is helpful to
    determine the a_nodes during the meshing process.
    It basically uses the relevant_cummulated_layup_boundaries with a flipped
    Layup definition and assigns the correct layer number.
    
    Args:         
        Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
                  #Start[-]	End[-] thickness [mm] Orientation [deg] MatID
            
    Returns: Projection: list of numpy arrays of the intervaltree structure 
            of the form:
                      start end layer#
            np.array([[ 0.3  0.7  3. ]
                      [ 0.7  1.   2. ]
                      [ 0.   0.3  2. ]])                
    """
    flipped_Projection = relevant_cummulated_layup_boundaries(np.flipud(Layup))
    b = np.flipud(np.arange(1, len(Layup) + 1))
    d = dict(enumerate(b))

    inverse_projection = []
    for iv_tree in flipped_Projection:
        part1 = iv_tree[:, :2]
        part2 = np.asarray([[d[l] for l in iv_tree[:, 2]]])
        inverse_projection.append(np.concatenate((part1, part2.T), axis=1))

    inverse_projection = reversed(inverse_projection)
    return sort_layup_projection(inverse_projection)

# ==============================================================================
if __name__ == "__main__":
    """Executes the following code if the file is exceuted as the main file 
    directly. This is used in this context to test the fuctions from abobe. """

    plt.close("all")
    Layup = np.array(
        [
            [0.0, 0.5, 0.23, 45.0, 1.0],
            [0.0, 1.0, 0.23, 45.0, 1.0],
            [0.3, 0.7, 0.6, 45.0, 2.0],
            [0.05, 0.6, 0.6, 45.0, 2.0],
            [0.85, 0.2, 0.6, 45.0, 2.0],
            [0.9, 0.4, 0.6, 45.0, 2.0],
            [0.4, 0.6, 0.6, 45.0, 2.0],
            [0.3, 0.7, 0.6, 45.0, 2.0],
            [0.20, 0.532, 0.6, 45.0, 2.0],
            [0.25, 0.3, 0.6, 45.0, 2.0],
        ],
    )

    plot_layup_projection(Layup)
    relevant_projectionlist = relevant_cummulated_layup_boundaries(Layup)
    inverse_projection = inverse_relevant_cummulated_layup_boundaries(Layup)
