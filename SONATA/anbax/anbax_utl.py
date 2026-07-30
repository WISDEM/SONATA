import numpy as np
from mpi4py import MPI

import b3_secfem
import basix.ufl
import ufl
from dolfinx import mesh as dmesh
from dolfinx.mesh import meshtags
from dolfinx import io

def dolfin_solve(cbm_mesh, cbm_nodes, cbm_materials):
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


    Notes
    ----------
    the cells of cbm_mesh already contain the nodes. So the information is
    currently passed twice. But consistent with the export_cells_for_vabs.

    """
    # Would like to avoid writing mesh to a file in the future, but for now just give dummy name
    path = "dolfinx_temp_mesh.xdmf"

    matLibrary = [cbm_materials[m].b3mat for m in cbm_materials]

    elem = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
    domain = ufl.Mesh(elem)

    coords = np.zeros((len(cbm_nodes), 2))
    for n in cbm_nodes:
        coords[n.id-1,:] = n.coordinates

    n_cells = len(cbm_mesh)
    tris = np.zeros((n_cells, 3), dtype=np.int32)
    fiber_orientations = np.zeros(n_cells)
    plane_orientations = np.zeros(n_cells)
    materials_vec = [None] * n_cells
    for c in cbm_mesh:
        tris[c.id-1,:] = [n.id-1 for n in c.nodes]
        plane_orientations[c.id-1] = c.theta_1[0]
        fiber_orientations[c.id-1] = c.theta_3
        materials_vec[c.id-1] = matLibrary[c.MatID-1]

    # Create the dolfinx mesh
    mesh = dmesh.create_mesh(MPI.COMM_WORLD, tris, domain, coords)

    # Any dolfinx mesh creation or write-read sequence reorders the elements,
    # so we have to keep up with the changes.
    # b3_secfem does this after it reads the mesh, so that step is taken care of, but need to do it here too
    oci = np.asarray(mesh.topology.original_cell_index)
    materials_vec = [materials_vec[m] for m in oci]
    fiber_orientations = fiber_orientations[oci]
    plane_orientations = plane_orientations[oci]

    # This meshtag doesn't do anything if not using the "region_mat" approach, but code complains if no tags are given
    indices = np.arange(n_cells, dtype=np.int32)
    values = np.array([c.id-1 for c in cbm_mesh], dtype=np.int32)
    cell_dim = mesh.topology.dim
    mesh.topology.create_connectivity(cell_dim, cell_dim)
    mt = meshtags(mesh, cell_dim, indices, values)
    mt.name = "cell_tags"

    # Write mesh to file
    with io.XDMFFile(mesh.comm, str(path), "w") as xf:
        xf.write_mesh(mesh)
        xf.write_meshtags(mt, mesh.geometry)

    # Create the b3_secfem section structure and return the solution
    mysec = b3_secfem.SectionInput(mesh_path=path,
                                   degree=1,
                                   per_cell_material=materials_vec,
                                   per_cell_alpha_deg=fiber_orientations,
                                   per_cell_beta_deg=plane_orientations)

    result = b3_secfem.solve( mysec )

    # Keep oci tracker consistent
    result.oci = oci[result.oci]
    # Store per-cell arrays (all already in dolfinx order) so recovery can
    # build per-cell material frames and detect isotropic vs orthotropic.
    result.per_cell_beta_deg  = plane_orientations
    result.per_cell_alpha_deg = fiber_orientations
    result.per_cell_material  = materials_vec

    return result


def anbax_unit_recovery(anba, T=None):
    """
    Function to recover unit stresses and strains from an applied loading
    Results are generated in global coordinates with dimensions: ([3 unit F + 3 unit M, num elements, 6 voigt])

    INPUTS:
    anba    -   dolfin construct from anbax
    T       -   Transformation matrix to convert results from ANBA to SONATA/VABS coordinates

    OUTPUTS:
    elem_stress_tran  -   global stress field in SONATA/VABS coordinates
    elem_stressM_mat  -   material-frame stress field before SONATA/VABS rotation
    elem_strain_tran  -   global strain field in SONATA/VABS coordinates
    elem_strainM_mat  -   material-frame strain field before SONATA/VABS rotation
    """

    fields  = b3_secfem.recover_unit_load_strains(anba) #unit_load_
    stress  = fields.sigma.copy()   # dim: ([3 unit F + 3 unitM, num elements, 6 voigt])
    strain  = fields.epsilon.copy() # dim: ([3 unit F + 3 unitM, num elements, 6 voigt])

    if T is None:
        T = np.eye(3)

    n_el = stress.shape[1]

    # Per-cell plane-orientation angles (theta_1) in dolfinx order.
    beta_arr  = anba.per_cell_beta_deg  if anba.per_cell_beta_deg  is not None else np.zeros(n_el)
    mat_list  = anba.per_cell_material  if anba.per_cell_material  is not None else [None] * n_el

    def _Q_sonata_material(t_deg):
        """Rotation Q whose columns are the SONATA material axes expressed in
        b3_secfem global (x, y, z):
          col 0 -> fiber axis    = z  (beam direction)
          col 1 -> perimeter     = (cos t, sin t, 0)
          col 2 -> through-thick = (-sin t, cos t, 0)
        so that  sigma_SONATA_material = Q.T @ sigma_b3global @ Q
        gives stress in (fiber, perimeter, through-thickness) convention.
        """
        t = np.deg2rad(t_deg)
        c, s = np.cos(t), np.sin(t)
        return np.array([[0.,  c, -s],
                         [0.,  s,  c],
                         [1.,  0.,  0.]])

    # Initialize the outputs:
    elem_stress_tran  = np.zeros(stress.shape)
    elem_strain_tran  = np.zeros(strain.shape)
    elem_stressM_mat  = np.zeros(stress.shape)
    elem_strainM_mat  = np.zeros(strain.shape)

    for k in range(n_el):
        mat = mat_list[k]
        is_isotropic = (mat is None or isinstance(mat, b3_secfem.IsotropicMaterial))

        if not is_isotropic:
            Q = _Q_sonata_material(float(beta_arr[k]))

        for ii in range(6):
            s = stress[ii, k, :]
            e = strain[ii, k, :]

            # Unpack Voigt -> 3x3 symmetric (b3_secfem Voigt: s11,s22,s33,s23,s13,s12)
            sg3 = np.array([[s[0], s[5], s[4]],
                            [s[5], s[1], s[3]],
                            [s[4], s[3], s[2]]])
            eg3 = np.array([[e[0], e[5], e[4]],
                            [e[5], e[1], e[3]],
                            [e[4], e[3], e[2]]])

            # Global stress in SONATA/VABS axes (T permutes b3_secfem x,y,z -> SONATA 2,3,1)
            istress = T.T @ sg3 @ T
            istrain = T.T @ eg3 @ T
            elem_stress_tran[ii, k, :] = np.r_[np.diag(istress), istress[1,2], istress[0,2], istress[0,1]]
            elem_strain_tran[ii, k, :] = np.r_[np.diag(istrain), istrain[1,2], istrain[0,2], istrain[0,1]]

            if is_isotropic:
                # Isotropic: no preferred material direction -> material frame = global frame
                sm3 = istress
                em3 = istrain
            else:
                # Orthotropic: rotate to SONATA material frame (fiber, perimeter, through-thickness)
                # Q.T @ sg3_b3global @ Q gives stress in SONATA material convention directly
                sm3 = Q.T @ sg3 @ Q
                em3 = Q.T @ eg3 @ Q

            elem_stressM_mat[ii, k, :] = np.r_[np.diag(sm3), sm3[1,2], sm3[0,2], sm3[0,1]]
            elem_strainM_mat[ii, k, :] = np.r_[np.diag(em3), em3[1,2], em3[0,2], em3[0,1]]

    oci2orig = np.argsort(anba.oci)
    elem_stress_tran = elem_stress_tran[:, oci2orig, :]
    elem_stressM_mat = elem_stressM_mat[:, oci2orig, :]
    elem_strain_tran = elem_strain_tran[:, oci2orig, :]
    elem_strainM_mat = elem_strainM_mat[:, oci2orig, :]

    return elem_stress_tran, elem_stressM_mat, elem_strain_tran, elem_strainM_mat


def anbax_recovery(anba, force, moment, T=None):
    """
    Function to recover total stresses and strains from an applied loading
    Results are generated in global coordinates with dimensions: ([num elements, 6 voigt])

    INPUTS:
    anba    -   dolfin construct from anbax
    force   -   Forces in anbax coordinates, [F1, F2, F3], e.g. force = [2.2, 3.4, 1.1]
    moment  -   Moments in anbax coordinates, [M1, M2, M3], e.g. moment = [4.2, 5.7, 6.2]
    T       -   Transformation matrix to convert results from ANBA to SONATA/VABS coordinates

    OUTPUTS:
    stress_sum  -   global total stress field in SONATA/VABS coordinates
    strain_sum  -   global total strain field in SONATA/VABS coordinates
    """

    stress, stressM, strain, strainM = anbax_unit_recovery(anba, T=T)
    n_el = stress.shape[1]

    # Rotate forces and moments
    if True or T is None:
        myforce, mymoment = force, moment
    else:
        myforce  = T.T @ force @ T
        mymoment = T.T @ moment @ T

    # Add up total stress and strain through super-position linear combo of unit forces and moments
    stress_sum  = np.zeros((n_el,6))
    stressM_sum = np.zeros((n_el,6))
    strain_sum  = np.zeros((n_el,6))
    strainM_sum = np.zeros((n_el,6))
    for k in range(3):
        stress_sum += myforce[ k] * stress[  k,:,:]
        stress_sum += mymoment[k] * stress[3+k,:,:]

        stressM_sum += myforce[ k] * stressM[  k,:,:]
        stressM_sum += mymoment[k] * stressM[3+k,:,:]

        strain_sum += myforce[ k] * strain[  k,:,:]
        strain_sum += mymoment[k] * strain[3+k,:,:]

        strainM_sum += myforce[ k] * strainM[  k,:,:]
        strainM_sum += mymoment[k] * strainM[3+k,:,:]

    return stress_sum, stressM_sum, strain_sum, strainM_sum


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
