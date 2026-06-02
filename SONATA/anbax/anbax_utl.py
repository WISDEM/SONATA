import numpy as np
from mpi4py import MPI

import b3_secfem
import basix.ufl
import ufl
from dolfinx import mesh as dmesh
from dolfinx.mesh import meshtags
from dolfinx import io

import pdb

def build_mat_library(cbm_materials):

    matLibrary = []

    maxE = 0.
    matid = 0
    matdict = {}
    for m in cbm_materials.values():
        if m.orth == 0:
            mat = b3_secfem.IsotropicMaterial(E=m.E, nu=m.nu, rho=m.rho, name=str(m.id))
        elif m.orth == 1:
            mat = b3_secfem.OrthotropicMaterial(E1=m.E[0], # Young's modulus, fibre [Pa] (Ezz along beam axis)
                                                E2=m.E[1], # Young's modulus, transverse-2 [Pa]
                                                E3=m.E[2], # Young's modulus, transverse-3 [Pa]
                                                G12=m.G[0], G13=m.G[1], G23=m.G[2],
                                                nu12=m.nu[0], nu13=m.nu[1], nu23=m.nu[2],
                                                rho=m.rho, name=str(m.id))
        elif m.orth == 2:
            raise ValueError('material type 2 (anysotropic) not supported by Anba')

            #mat = material.
        maxE = np.maximum(maxE, np.array(m.E).max())
        matLibrary.append(mat)
        matdict[m.id] = matid
        matid += 1
    return matLibrary, matdict, maxE

def build_dolfin_mesh(cbm_mesh, cbm_nodes, cbm_materials):
    """function to generate the dolfin.Mesh from a SONATA-CBM definition to run
    with anbax

    Parameters
    ----------
    cbm_mesh : list of cell instances
        from the SONATA-CBM preprocessor

    cbm_nodes : list of nodes
        from the SONATA-CBM preprocessor


    Returns
    ----------
    mesh : dolfin.Mesh
    matLibrary : vector of anbax materials
    materials : dolfin.MeshFunction definign cell materials
    plane_orientations : dolfin.MeshFunction defining cell plane orientations
    fiber_orientations : dolfin.MeshFunction defining cell material fiber orientation
    maxE : reference elastic modulus for scaling rigid mode constraints


    Notes
    ----------
    the cells of cbm_mesh already contain the nodes. So the information is
    currently passed twice. But consistent with the export_cells_for_vabs.

    """
    path = "temp.xdmf"
    (matLibrary, matdict, maxE) = build_mat_library(cbm_materials)

    elem = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
    domain = ufl.Mesh(elem)

    
    coords = np.array([n.coordinates for n in cbm_nodes])
    tris = np.zeros((len(cbm_mesh), 3), dtype=np.int32)
    for ic, c in enumerate(cbm_mesh):
        tris[ic,:] = [n.id-1 for n in c.nodes]
    mesh = dmesh.create_mesh(MPI.COMM_WORLD, tris, domain, coords)

    fiber_orientations_vec = np.array([c.theta_3 for c in cbm_mesh])
    plane_orientations_vec = np.array([c.theta_1[0] for c in cbm_mesh])
    
    cell_dim = mesh.topology.dim
    cell_tags = np.array([c.MatID for c in cbm_mesh])
    n_cells = cell_tags.size
    indices = np.arange(n_cells, dtype=np.int32)
    values = cell_tags.astype(np.int32)
    mesh.topology.create_connectivity(cell_dim, cell_dim)
    mt = meshtags(mesh, cell_dim, indices, values)
    mt.name = "cell_tags"
    with io.XDMFFile(mesh.comm, str(path), "w") as xf:
        xf.write_mesh(mesh)
        xf.write_meshtags(mt, mesh.geometry)

    rmats = {}
    for c in cbm_mesh:
        #pdb.set_trace()
        rmats[c.id-1] = b3_secfem.RegionMat(material=matLibrary[c.MatID-1])#,
                                            #alpha_deg=float(c.theta_1[0]),
                                            #beta_deg=float(c.theta_3))
    mysec = b3_secfem.SectionInput(mesh_path=path,
                                   region_materials=rmats,
                                   per_cell_beta_deg=fiber_orientations_vec,
                                   per_cell_alpha_deg=plane_orientations_vec)
        
    '''
    cell_map = mesh.topology.index_map(mesh.topology.dim)
    n_cells = cell_map.size_local + cell_map.num_ghosts
    cell_idx = np.arange(n_cells, dtype=np.int32)
    
    _ = meshtags(mesh, mesh.topology.dim, cell_idx, [matdict[cbm_mesh[0].MatID]]*n_cells)
    
    _ = meshtags(mesh, mesh.topology.dim, cell_idx, plane_orientations_vec)

    _ = meshtags(mesh, mesh.topology.dim, cell_idx, fiber_orientations_vec)

    b3_secfem.write_xdmf(path, mesh)

    rmats = {}
    for c in cbm_mesh:
        #pdb.set_trace()
        rmats[c.id-1] = b3_secfem.RegionMat(material=matLibrary[c.MatID-1],
                                            beta_deg=float(c.theta_1[0]),
                                            alpha_deg=float(c.theta_3))

    mysec = b3_secfem.SectionInput(mesh_path=path,
                                   region_materials=rmats,
                                   per_cell_beta_deg=plane_orientations_vec,
                                   per_cell_alpha_deg=fiber_orientations_vec)
    '''
    return mysec


def anbax_recovery(anba, n_el, force, moment, voigt_convention, T):
    """
    Function to recover stresses and strains from an applied loading
    Results are generated in the global and local ('M' for material coordinate system) coordinates

    INPUTS:
    anba    -   dolfin construct from anbax
    n_el    -   number of mesh elements
    force   -   Forces in anbax coordinates, [F1, F2, F3], e.g. force = [2.2, 3.4, 1.1]
    moment  -   Moments in anbax coordinates, [M1, M2, M3], e.g. moment = [4.2, 5.7, 6.2]
    voigt_convention    -   "anba" with [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy] or "paraview" with [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]
    T       -   Transformation matrix to convert results from ANBA to SONATA/VABS coordinates

    OUTPUTS:
    *_tran issues that outputs were converted to the SONATA/VABS coordinates
    remove '_tran' from the following output names when needed in ANBA coordinates

    tmp_StressF_tran    -   global stress field
    tmp_StressF_M_tran  -   local stress field
    tmp_StrainF_tran    -   global strain field
    tmp_StrainF_M_tran  -   local strain field



    """

    fields = b3_secfem.recover_strains(anba)
    tmp_StressF_vec = fields.sigma
    tmp_StrainF_vec = fields.epsilon
    
    # Stress field
    #anba.stress_field(force, moment, reference="global", voigt_convention=voigt_convention)  # get stress field in global sys
    #tmp_StressF_vec = np.array(anba.STRESS.vector().vec())  # global stress field
    #anba.stress_field(force, moment, reference="local", voigt_convention=voigt_convention)  # get stress field in local sys (material coordinates)
    #tmp_StressF_M_vec = np.array(anba.STRESS.vector().vec())  # local stress field

    # Strain field
    #anba.strain_field(force, moment, reference="global", voigt_convention=voigt_convention)  # get strain field in global sys
    #tmp_StrainF_vec = np.array(anba.STRAIN.vector().vec())  # global strain field
    #anba.strain_field(force, moment, reference="local", voigt_convention=voigt_convention)  # get strain field in local sys (material coordinates)
    #tmp_StrainF_M_vec = np.array(anba.STRAIN.vector().vec())  # local strain field

    # cd = anba.STRESS.function_space().dofmap().cell_dofs  # index numbers of cells from dolfin mesh that was used for stress recovery (each cell has 6 dofs)

    s_11 = np.zeros(n_el)
    s_22 = np.zeros(n_el)
    s_33 = np.zeros(n_el)
    s_23 = np.zeros(n_el)
    s_13 = np.zeros(n_el)
    s_12 = np.zeros(n_el)
    tmp_StressF = np.zeros((n_el, 3, 3))
    tmp_StressF_tran = np.zeros((n_el, 3, 3))

    s_11_M = np.zeros(n_el)
    s_22_M = np.zeros(n_el)
    s_33_M = np.zeros(n_el)
    s_23_M = np.zeros(n_el)
    s_13_M = np.zeros(n_el)
    s_12_M = np.zeros(n_el)
    tmp_StressF_M = np.zeros((n_el, 3, 3))
    tmp_StressF_M_tran = np.zeros((n_el, 3, 3))

    e_11 = np.zeros(n_el)
    e_22 = np.zeros(n_el)
    e_33 = np.zeros(n_el)
    e_23 = np.zeros(n_el)
    e_13 = np.zeros(n_el)
    e_12 = np.zeros(n_el)
    tmp_StrainF = np.zeros((n_el, 3, 3))
    tmp_StrainF_tran = np.zeros((n_el, 3, 3))

    e_11_M = np.zeros(n_el)
    e_22_M = np.zeros(n_el)
    e_33_M = np.zeros(n_el)
    e_23_M = np.zeros(n_el)
    e_13_M = np.zeros(n_el)
    e_12_M = np.zeros(n_el)
    tmp_StrainF_M = np.zeros((n_el, 3, 3))
    tmp_StrainF_M_tran = np.zeros((n_el, 3, 3))

    # cell_id = np.zeros((n_el, 6))

    if voigt_convention == "anba":  # [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy]
        for i in range(n_el):
            # stresses in "global" system
            s_11[i] = tmp_StressF_vec[i * 6]
            s_22[i] = tmp_StressF_vec[i * 6 + 1]
            s_33[i] = tmp_StressF_vec[i * 6 + 2]
            s_23[i] = tmp_StressF_vec[i * 6 + 3]  # equiv to s_23
            s_13[i] = tmp_StressF_vec[i * 6 + 4]  # equiv to s_31
            s_12[i] = tmp_StressF_vec[i * 6 + 5]  # equiv to s_21
            tmp_StressF[i, :, :] = np.array([[s_11[i], s_12[i], s_13[i]], [s_12[i], s_22[i], s_23[i]], [s_13[i], s_23[i], s_33[i]]])
            tmp_StressF_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StressF[i]), T)  # transform to sonata coordinate system
            # stresses in "local" system
            #s_11_M[i] = tmp_StressF_M_vec[i * 6]
            #s_22_M[i] = tmp_StressF_M_vec[i * 6 + 1]
            #s_33_M[i] = tmp_StressF_M_vec[i * 6 + 2]
            #s_23_M[i] = tmp_StressF_M_vec[i * 6 + 3]  # equiv to s_23_M
            #s_13_M[i] = tmp_StressF_M_vec[i * 6 + 4]  # equiv to s_31_M
            #s_12_M[i] = tmp_StressF_M_vec[i * 6 + 5]  # equiv to s_21_M
            #tmp_StressF_M[i, :, :] = np.array([[s_11_M[i], s_12_M[i], s_13_M[i]], [s_12_M[i], s_22_M[i], s_23_M[i]], [s_13_M[i], s_23_M[i], s_33_M[i]]])
            #tmp_StressF_M_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StressF_M[i]), T)  # transform to sonata coordinate system

            # strains in "global" system
            e_11[i] = tmp_StrainF_vec[i * 6]
            e_22[i] = tmp_StrainF_vec[i * 6 + 1]
            e_33[i] = tmp_StrainF_vec[i * 6 + 2]
            e_23[i] = tmp_StrainF_vec[i * 6 + 3]  # equiv to e_23
            e_13[i] = tmp_StrainF_vec[i * 6 + 4]  # equiv to e_31
            e_12[i] = tmp_StrainF_vec[i * 6 + 5]  # equiv to e_21
            tmp_StrainF[i, :, :] = np.array([[e_11[i], e_12[i], e_13[i]], [e_12[i], e_22[i], e_23[i]], [e_13[i], e_23[i], e_33[i]]])
            tmp_StrainF_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StrainF[i]), T)  # transform to sonata coordinate system
            # strains in "local" system
            #e_11_M[i] = tmp_StrainF_M_vec[i * 6]
            #e_22_M[i] = tmp_StrainF_M_vec[i * 6 + 1]
            #e_33_M[i] = tmp_StrainF_M_vec[i * 6 + 2]
            #e_23_M[i] = tmp_StrainF_M_vec[i * 6 + 3]  # equiv to e_23_M
            #e_13_M[i] = tmp_StrainF_M_vec[i * 6 + 4]  # equiv to e_31_M
            #e_12_M[i] = tmp_StrainF_M_vec[i * 6 + 5]  # equiv to e_21_M
            #tmp_StrainF_M[i, :, :] = np.array([[e_11_M[i], e_12_M[i], e_13_M[i]], [e_12_M[i], e_22_M[i], e_23_M[i]], [e_13_M[i], e_23_M[i], e_33_M[i]]])
            #tmp_StrainF_M_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StrainF_M[i]), T)  # transform to sonata coordinate system

            # cell_id[i, :] = cd(i)
    elif voigt_convention == "paraview":  # different ordering compared to "anba"; "paraview" ordering: [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]
        print("ToDo - Process to paraview output")

    # Export to Paraview format (to be tested!)
    # file_res = do.XDMFFile('output_filename.xdmf')
    # file_res.parameters['functions_share_mesh'] = True
    # file_res.parameters['rewrite_function_mesh'] = False
    # file_res.parameters["flush_output"] = True
    # file_res.write(anba.STRESS, t=2)  # t=unique_number


    return tmp_StressF_tran, tmp_StressF_M_tran, tmp_StrainF_tran, tmp_StrainF_M_tran



def ComputeShearCenter(stiff_matrix):  # shear center equiv. to elastic axes
    K1 = np.array([[stiff_matrix[i, j] for j in range(3)] for i in range(3)])
    K3 = np.array([[stiff_matrix[i, j+3] for j in range(3)] for i in range(3)])
    Y = np.linalg.solve(K1, -K3)
    return [-Y[2, 0], Y[1, 0]]
    # return [-Y[1,2], Y[0,2]]

def ComputeTensionCenter(stiff_matrix):  # tension center equiv. to neutral axes
    K1 = np.array([[stiff_matrix[i, j] for j in range(3)] for i in range(3)])
    K3 = np.array([[stiff_matrix[i, j+3] for j in range(3)] for i in range(3)])
    Y = np.linalg.solve(K1, -K3)
    return [Y[0, 2], -Y[0, 1]]
    # return [Y[2,1], -Y[2,0]]

def ComputeMassCenter(mass_matrix):
    M1 = np.array([[mass_matrix[i, j] for j in range(3)] for i in range(3)])
    M3 = np.array([[mass_matrix[i, j+3] for j in range(3)] for i in range(3)])
    Y = np.linalg.solve(M1, -M3)
    return [Y[0,2], Y[1,0]]
    # return [Y[2,1], -Y[2,0]]

def DecoupleStiffness(stiff_matrix):
    K = np.array([[stiff_matrix[i, j] for j in range(6)] for i in range(6)])
    K1 = np.array([[stiff_matrix[i, j] for j in range(3)] for i in range(3)])
    K3 = np.array([[stiff_matrix[i, j+3] for j in range(3)] for i in range(3)])
    Y = np.linalg.solve(K1, -K3)
    I3 = np.eye(3)
    Z3 = np.zeros((3,3))
    TL = np.block([[I3, Z3], [Y.T, I3]])
    TR = np.block([[I3, Y], [Z3, I3]])
    return TL @ K @ TR

def PrincipalAxesRotationAngle(decoupled_stiff_matrix):
    K1 = np.array([[decoupled_stiff_matrix[i, j] for j in range(3)] for i in range(3)])
    K3 = np.array([[decoupled_stiff_matrix[i+3, j+3] for j in range(3)] for i in range(3)])
    (w1, v1) = np.linalg.eig(K1)
    (w3, v3) = np.linalg.eig(K3)

    if np.abs(v3[0,0]) < np.abs(v3[0,1]):
        angle = np.arccos(v3[0,0])
    else:
        angle = -np.arcsin(v3[0,1])
    return angle
