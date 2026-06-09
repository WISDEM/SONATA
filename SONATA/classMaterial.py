#------------------------------------------------------------------------------------
# Import Libraries
#------------------------------------------------------------------------------------

# Standard library
from collections import OrderedDict

# Third party modules
import numpy as np

import b3_secfem


#------------------------------------------------------------------------------------
# Material class
#------------------------------------------------------------------------------------
class Material:
    """
    general material class

    Attributes
    ----------
    id : int
        material identifier

    name : str
        short material name

    description : str
        description of the material, e.g.: unidirektional ht-carbon fiber
        composite with epoxy matrix (FVC of 60%)

    source : str
        source of the material properties e.g.: elasitc properties derived
        from Schuermann, (p.184, 417) with a semi-empiric Puck approach

    orth : int
        orth is the flag to indicate whether the material is isotropic (0),
        orthotropic (1) or general anisotropic (2) in consitency with VABS
        Manual for Users (2011)

    rho : float
        density in kg/m**3

    """

    __slots__ = ("id", "name", "description", "source", "orth", "rho", "b3mat")

    def __init__(self, ID="NOID", name="noname", description="nodescription", source="nosource", orth=None, rho=0, **kwargs):
        self.id = ID
        self.name = name
        self.description = description
        self.source = source
        self.orth = orth
        self.rho = float(rho)
        self.b3mat = None

    def __repr__(self):
        if self.orth == 0:
            return str("%s: IsotropicMaterial: %s" % (self.id, self.name))
        elif self.orth == 1:
            return str("%s: OrthotropicMaterial: %s" % (self.id, self.name))
        elif self.orth == 2:
            return str("%s: OrthotropicMaterial: %s" % (self.id, self.name))
        else:
            return str("%s: UndefinedMaterial: %s" % (self.id, self.name))

    def constitutive_tensor(self):
        """
        Calculate the local consitutive tensor for the material.

        This is just a template function to be overridden by specific
        material types.

        Returns
        -------
        constitutive_tensor : (6,6) numpy.ndarray
            Local constitutive tensor. Template just returns zeros.

        Notes
        -----

        For consistency, shear strain components (gamma_ij) are the engineering
        shear strain
        that is twice the components of the elasticity tensor shear strains
        (e.g., gamma_ij = 2*eps_ij = eps_ij + eps_ji) for use with this tensor.

        """
        if self.b3mat is not None:
            return self.b3mat.C_local()
        else:
            return np.zeros((6,6))

    def rotated_constitutive_tensor(self, plane_orientation, fiber_orientation):
        """
        Transforms the material constitutive tensor from material to global
        coordinate system.

        Parameters
        ----------
        plane_orientation : float
            Rotation parameter for in plane orientation of the fiber.
            Units are degrees.
        fiber_orientation : float
            Rotation parameter for fiber orientation (has not been
            significantly used/tested).
            Units are degrees.

        Returns
        -------
        tensor_global : (6,6) numpy.ndarray
            Elasticity constitutive tensor for converting from strains
            to stresses. See Notes for details.
            This is calculated in the global coordinates.

        Notes
        -----

        The returned tensor is size 6x6 and converts from strains in the form:
        [eps11, eps22, eps33, gamma23, gamma13, gamma12]
        to stresses of
        [sigma11, sigma22, sigma33, sigma23, sigma13, sigma12].

        The shear strain components (gamma_ij) are the engineering shear strain
        that is twice the components of the elasticity tensor shear strains
        (e.g., gamma_ij = 2*eps_ij = eps_ij + eps_ji)

        The implementation is heavily copied from:
            anba4/anba4/material/material.cpp/TransformationMatrix
        This is added to the python implementation for easier access when doing
        extra calculations for viscoelastic materials.

        """

        # material coordinate tensor
        tensor_material = self.constitutive_tensor()

        # Angles
        pi180 = np.pi / 180.

        # alpha->fiber plane oriention; beta->fiber oriention.
        alpha = plane_orientation
        beta = fiber_orientation

        # calculate rotation matrix, this is heavily copied from the c++ code.
        sn_a = -np.sin(alpha*pi180)
        cn_a =  np.cos(alpha*pi180)
        sn_b = -np.sin(beta*pi180)
        cn_b =  np.cos(beta*pi180)

        transformMatrix = np.zeros((6,6))

        transformMatrix[0, 0] = cn_a * cn_a * cn_b * cn_b
        transformMatrix[0, 1] = sn_a * sn_a
        transformMatrix[0, 2] = cn_a * cn_a * sn_b * sn_b
        transformMatrix[0, 3] = -2.0 * cn_a * sn_a * sn_b
        transformMatrix[0, 4] = -2.0 * cn_a * cn_a * sn_b * cn_b
        transformMatrix[0, 5] = 2.0 * cn_a * sn_a * cn_b

        transformMatrix[1, 0] = sn_a * sn_a * cn_b * cn_b
        transformMatrix[1, 1] = cn_a * cn_a
        transformMatrix[1, 2] = sn_a * sn_a * sn_b * sn_b
        transformMatrix[1, 3] = 2.0 * cn_a * sn_a * sn_b
        transformMatrix[1, 4] = -2.0 * sn_a * sn_a * sn_b * cn_b
        transformMatrix[1, 5] = -2.0 * cn_a * sn_a * cn_b

        transformMatrix[2, 0] = sn_b * sn_b
        transformMatrix[2, 2] = cn_b * cn_b
        transformMatrix[2, 4] = 2.0 * cn_b * sn_b

        transformMatrix[3, 0] = -sn_a * sn_b * cn_b
        transformMatrix[3, 2] = sn_a * sn_b * cn_b
        transformMatrix[3, 3] = cn_a * cn_b
        transformMatrix[3, 4] = -sn_a * cn_b * cn_b + sn_a * sn_b * sn_b
        transformMatrix[3, 5] = cn_a * sn_b

        transformMatrix[4, 0] = cn_a * sn_b * cn_b
        transformMatrix[4, 2] = -cn_a * sn_b * cn_b
        transformMatrix[4, 3] = sn_a * cn_b
        transformMatrix[4, 4] = -cn_a * sn_b * sn_b + cn_a* cn_b * cn_b
        transformMatrix[4, 5] = sn_a * sn_b

        transformMatrix[5, 0] = -sn_a * cn_a * cn_b * cn_b
        transformMatrix[5, 1] = cn_a * sn_a
        transformMatrix[5, 2] = -sn_a * cn_a * sn_b * sn_b
        transformMatrix[5, 3] = -cn_a * cn_a * sn_b + sn_a *sn_a * sn_b
        transformMatrix[5, 4] = 2.0 * sn_a * sn_b * cn_a * cn_b
        transformMatrix[5, 5] = cn_a * cn_a * cn_b - sn_a * sn_a * cn_b


        # transform from local -> global
        tensor_global = transformMatrix @ tensor_material @ transformMatrix.T

        return tensor_global


#------------------------------------------------------------------------------------
# IsotropicMaterial class
#------------------------------------------------------------------------------------
class IsotropicMaterial(Material):
    """
    Isotropic Material


    Attributes
    ----------
    E : float
        in GPa; Young's modulus

    nu : float
        nondimensional; Poisson's ratio

    alpha : np.ndarray
        in 1/K; coefficient of thermal expansion in direction

    YS :  float
        in N/m**2; Yield Strenth (Streckgrenze)
    UTS : float
        in N/m**2; Ultimate Tensile Strenght (Zugfestigkeit)

    viscoelastic :
        Properties associated with different time scales of the material to
        define viscoelasticity.
    """

    __slots__ = ("E", "nu", "alpha", "YS", "UTS", "viscoelastic")

    def __init__(self, **kw):
        kw["orth"] = 0
        super().__init__(**kw)
        self.E = None
        self.nu = None
        self.alpha = None
        self.YS = None
        self.UTS = None
        
        if kw.get("E") is not None:
            self.E = float(kw.get("E"))

        if kw.get("nu") is not None:
            self.nu = float(kw.get("nu"))

        if kw.get("alpha") is not None:
            self.alpha = float(kw.get("alpha"))

        if kw.get("YS") is not None:
            self.YS = float(kw.get("YS"))

        if kw.get("UTS") is not None:
            self.UTS = float(kw.get("UTS"))

        self.viscoelastic = {}
        self.b3mat = b3_secfem.IsotropicMaterial(E=self.E, nu=self.nu, rho=self.rho, name=self.name)

    '''
    # Using the abstract version that calls the b3_secfem materials class instead
    def constitutive_tensor(self):
        """
        Calculate the local consitutive tensor for the material.

        Returns
        -------
        constitutive_tensor : (6,6) numpy.ndarray
            Elasticity constitutive tensor for converting from strains
            to stresses. See Notes for details.
            This is calculated in the local material coordinates.

        Notes
        -----

        The returned tensor is size 6x6 and converts from strains in the form:
        [eps11, eps22, eps33, gamma23, gamma13, gamma12]
        to stresses of
        [sigma11, sigma22, sigma33, sigma23, sigma13, sigma12].

        The shear strain components (gamma_ij) are the engineering shear strain
        that is twice the components of the elasticity tensor shear strains
        (e.g., gamma_ij = 2*eps_ij = eps_ij + eps_ji)

        The implementation is heavily copied from:
            anba4/anba4/material/material.cpp/IsotropicMaterial
        This is added to the python implementation for easier access when doing
        extra calculations for viscoelastic materials.

        """

        E = self.E
        nu = self.nu
        G = E / (2 * (1 + nu))

        delta = E / (1. + nu) / (1 - 2.*nu)
        diag = (1. - nu) * delta
        off_diag = nu * delta

        constitutive_tensor = np.zeros((6, 6))

        constitutive_tensor[0, 0] = diag
        constitutive_tensor[0, 1] = off_diag
        constitutive_tensor[0, 2] = off_diag

        constitutive_tensor[1, 0] = off_diag
        constitutive_tensor[1, 1] = diag
        constitutive_tensor[1, 2] = off_diag

        constitutive_tensor[2, 0] = off_diag
        constitutive_tensor[2, 1] = off_diag
        constitutive_tensor[2, 2] = diag

        constitutive_tensor[3, 3] = G
        constitutive_tensor[4, 4] = G
        constitutive_tensor[5, 5] = G

        return constitutive_tensor
    '''

#------------------------------------------------------------------------------------
# OrthotropicMaterial class
#------------------------------------------------------------------------------------
class OrthotropicMaterial(Material):
    """
    Orthotropic Material


    Attributes
    ----------
    E : np.ndarray
        in GPa; [E_1, E_2, E_3], with E_i: axial tensile modules in direction i
        and E_2 and E_3 the transverse tensile modules respectively

    G : np.ndarray
        in GPa; [G_12, G_13, G_23], with G_ij, is the shear modulus in
        direction j on the plane whose normal is in direction  i; for
        transversal insotropic materials G_13 = G_12

    nu : np.ndarray
        nondimensional; [nu12, nu_13, nu_23], nu_ij is the Poisson's ratio that
        corresponds to a contraction in direction j when an extension is
        applied in direction i.

    alpha_11 : np.ndarray
        in 1/K; [alpha_11, alpha_22, alpha_33], alpha_ii is the coefficient of
        thermal expansion in direction ii

    Xt :  float
        in N/m**2; 0° tensile strenght
    Xc : float
        in N/m**2; 0° compressive strenght
    Yt : float
        in N/m**2; 90° tensile strenght
    Yc : float
        in N/m**2; 90° compressive strenght
    S21 :
        in N/m**2; in-/out of plane shear strength

    viscoelastic :
        Properties associated with different time scales of the material to
        define viscoelasticity.
    """

    __slots__ = ('E', 'G', 'nu', 'alpha', 'Xt', 'Xc', 'Yt', 'Yc', 'S21', 'S23',
                 'viscoelastic')

    def __init__(self, flag_mat, **kw):
        kw['orth'] = 1
        Material.__init__(self, **kw)
        self.E = None
        self.G = None
        self.nu = None
        self.alpha = None
        self.Xt = None
        self.Xc = None
        self.Yt = None
        self.Yc = None
        self.S21 = None

        if kw.get('E') is not None:
            self.E = np.asarray(kw.get('E')).astype(float)

        if kw.get('G') is not None:
            self.G = np.asarray(kw.get('G')).astype(float)

        if kw.get('nu') is not None:
            self.nu = np.asarray(kw.get('nu')).astype(float)

        if kw.get('alpha') is not None:
            self.alpha = np.asarray(kw.get('alpha')).astype(float)

        if all(k in kw for k in ('E_1', 'E_2', 'E_3')):
            self.E = np.array([kw.get('E_1'), kw.get('E_2'), kw.get('E_3')]).astype(float)

        if all(k in kw for k in ('G_12', 'G_13', 'G_23')):
            self.G = np.array([kw.get('G_12'), kw.get('G_13'), kw.get('G_23')]).astype(float)

        if all(k in kw for k in ('nu_12', 'nu_13', 'nu_23')):
            self.nu = np.array([kw.get('nu_12'), kw.get('nu_13'), kw.get('nu_23')]).astype(float)

        if all(k in kw for k in ('alpha_11', 'alpha_22', 'alpha_33')):
            self.alpha = np.array([kw.get('alpha_11'), kw.get('alpha_22'), kw.get('alpha_33')]).astype(float)

        if flag_mat:  # wisdem includes vectors for the following material properties that are to be converted in order to comply with SONATA and VABS/anbax
            if kw.get('Xt') is not None:
                self.Xt = float(kw.get('Xt')[0])  # retrieve axial tensile strength in [MPa] from provided 3D vector

            if kw.get('Xc') is not None:
                self.Xc = float(kw.get('Xc')[0])  # retrieve axial compression strength in [MPa] from provided 3D vector

            if kw.get('Yt') is not None:
                self.Yt = float(kw.get('Xt')[1])  # retrieve transverse tensile strength in [MPa] from provided 3D vector

            if kw.get('Yc') is not None:
                self.Yc = float(kw.get('Xc')[1])  # retrieve transverse compression strength in [MPa] from provided 3D vector

            if kw.get('S') is not None:
                self.S21 = float(kw.get('S')[0])  # retrieve in-/out of plane shear strength [MPa] in [MPa] from provided 3D vector

        else:
            if kw.get('Xt') is not None:
                self.Xt = float(kw.get('Xt'))  # Axial Tensile Strength in [MPa]

            if kw.get('Xc') is not None:
                self.Xc = float(kw.get('Xc'))  # Axial Compression Strength  [MPa]

            if kw.get('Yt') is not None:
                self.Yt = float(kw.get('Yt'))  # Transverse Tensile strenght  [MPa]

            if kw.get('Yc') is not None:
                self.Yc = float(kw.get('Yc'))  # Transverse  Compression strenght  [Mpa]

            if kw.get('S21') is not None:
                self.S21 = float(kw.get('S21'))  # in-/out of plane shear strength [MPa]

        # self.S23 = float(kw.get('S23'))

        self.viscoelastic = {}
        
        # Default b3_secfem coordinate system with fiber axis on E1
        #E1, E2, E3 = self.E #[0], m.E[1], m.E[2]
        #G12, G13, G23 = self.G #[0], m.G[1], m.G[2]
        #nu12, nu13, nu23 = self.nu #[0], m.nu[1], m.nu[2]

        # If trying to trick b3_secfem to give results like Anba that has fiber axis on E3:
        E1, E2, E3 = self.E[1], self.E[2], self.E[0]
        G12, G13, G23 = self.G[2], self.G[0], self.G[1]
        nu12, nu31, nu32 = self.nu[2], self.nu[0], self.nu[1]
        nu13 = nu31*E1/E3
        nu23 = nu32*E2/E3
        
        self.b3mat = b3_secfem.OrthotropicMaterial(E1=E1, # Young's modulus, fibre [Pa] (Ezz along beam axis)
                                            E2=E2, # Young's modulus, transverse-2 [Pa]
                                            E3=E3, # Young's modulus, transverse-3 [Pa]
                                            G12=G12, G13=G13, G23=G23,
                                            nu12=nu12, nu13=nu13, nu23=nu23,
                                            rho=self.rho, name=self.name)

    '''
    # Using the abstract version that calls the b3_secfem materials class instead
    def constitutive_tensor(self):
        """
        Calculate the local consitutive tensor for the material.

        Returns
        -------
        constitutive_tensor : (6,6) numpy.ndarray
            Elasticity constitutive tensor for converting from strains
            to stresses. See Notes for details.
            This is calculated in the local material coordinates.

        Notes
        -----

        The returned tensor is size 6x6 and converts from strains in the form:
        [eps11, eps22, eps33, gamma23, gamma13, gamma12]
        to stresses of
        [sigma11, sigma22, sigma33, sigma23, sigma13, sigma12].

        The shear strain components (gamma_ij) are the engineering shear strain
        that is twice the components of the elasticity tensor shear strains
        (e.g., gamma_ij = 2*eps_ij = eps_ij + eps_ji)

        The implementation is heavily copied from:
            anba4/anba4/material/material.cpp/OrthotropicMaterial
        This is added to the python implementation for easier access when doing
        extra calculations for viscoelastic materials.

        Directions should be consistent with:
            https://ieawindsystems.github.io/windIO/main/source/detailed_turbine_documentation.html#materials

        """

        # using tensor direction indices z=1, x=2, y=3
        # (anba uses this direction ordering to return stress/strain)
        # Implementation should match anbax_util.py > build_mat_library.
        #
        # This means elastic modulus order should be:
        #   [Along Beam, Along Perimeter, Through Thickness]
        # Assuming no fiber orientation rotation angle.

        # If hacking b3_secfem to look
        e_xx = self.E[1]
        e_yy = self.E[2]
        e_zz = self.E[0]
        g_yz = self.G[1] # [G12, G13, G23][1] = G13 = yz
        g_xz = self.G[0] # [G12, G13, G23][0] = G12 = xz
        g_xy = self.G[2] # [G12, G13, G23][2] = G23 = xy
        nu_zy = self.nu[1] # [nu12, nu13, nu23][1] = nu13 = zy =/= yz
        nu_zx = self.nu[0] # [nu12, nu13, nu23][0] = nu12 = zx =/= xz
        nu_xy = self.nu[2] # [nu12, nu13, nu23][2] = nu23 = xy =/= yx

        # Calculate the other 3 poisson ratios.
        nu_yx = e_yy * nu_xy / e_xx
        nu_xz = e_xx * nu_zx / e_zz
        nu_yz = e_yy * nu_zy / e_zz

        constitutive_tensor = np.zeros((6, 6))

        delta = (1.0
                 - nu_xy * nu_yx
                 - nu_yz * nu_zy
                 - nu_xz * nu_zx
                 -2.0 * nu_yx * nu_zy * nu_xz) / (e_xx * e_yy * e_zz)

        constitutive_tensor[0, 0] = (1.0-nu_yz*nu_zy)/(e_yy*e_zz*delta)
        constitutive_tensor[0, 1] = (nu_xy+nu_zy*nu_xz)/(e_xx*e_zz*delta)
        constitutive_tensor[0, 2] = (nu_xz+nu_xy*nu_yz)/(e_xx*e_yy*delta)

        constitutive_tensor[1, 0] = constitutive_tensor[0, 1]
        constitutive_tensor[1, 1] = (1-nu_xz*nu_zx)/(e_xx*e_zz*delta)
        constitutive_tensor[1, 2] = (nu_yz+nu_yx*nu_xz)/(e_xx*e_yy*delta)

        constitutive_tensor[2, 0] = constitutive_tensor[0, 2]
        constitutive_tensor[2, 1] = constitutive_tensor[1, 2]
        constitutive_tensor[2, 2] = (1-nu_xy*nu_yx)/(e_xx*e_yy*delta)

        constitutive_tensor[3, 3] = g_yz
        constitutive_tensor[4, 4] = g_xz
        constitutive_tensor[5, 5] = g_xy

        return constitutive_tensor
    '''


#------------------------------------------------------------------------------------
# Utility functions
#------------------------------------------------------------------------------------
def read_materials(yml, viscoelastic_yaml=None):
    """
    Read material definitions from YAML data and create material objects.

    Parameters
    ----------
    yml : list
        List of dictionaries containing material definitions from YAML file.
        Each dictionary should contain material properties and optionally an 'id' field.

    Returns
    -------
    OrderedDict
        Dictionary of material objects keyed by material ID, sorted by ID.
        Contains IsotropicMaterial or OrthotropicMaterial instances based on 'orth' flag.

    Notes
    -----
    If no 'id' is provided in the material definition, materials are automatically
    assigned sequential IDs starting from 1. The flag_materials_vector_input is set
    to True in this case to handle WISDEM-specific material definitions.
    """
    materials = OrderedDict()
    for i, mat in enumerate(yml):
        if 'id' in mat:
            ID = mat['id']
            # If no ID is issued for materials, automatically define the following flag to allocate wisdem specfic material definitions
            flag_materials_vector_input = False
        else:
            ID = i + 1
            # use this in case no mat ID was issues
            # TODO: better create a user defined flag
            flag_materials_vector_input = True

        if mat.get('orth') == 0:
            materials[ID] = IsotropicMaterial(ID=ID, **mat)
        elif mat.get('orth') == 1:
            materials[ID] = OrthotropicMaterial(ID=ID, flag_mat=flag_materials_vector_input, **mat)

        if viscoelastic_yaml is not None:
            # read viscoelastic properties from separate yaml file and assign to material
            import yaml
            with open(viscoelastic_yaml, 'r') as f:
                viscoelastic_data = yaml.safe_load(f)

            for ve_mat in viscoelastic_data.get('materials', []):
                ve_mat_name = ve_mat.get('name', '').lower()
                if ve_mat_name == materials[ID].name.lower():
                    # assign viscoelastic properties to the material
                    viscoelastic_mat_keys = ['time_scales_v', 'E_v',
                                             'E_1_v', 'E_2_v', 'E_3_v',
                                             'G_12_v', 'G_13_v', 'G_23_v']

                    for k in viscoelastic_mat_keys:
                        if ve_mat.get(k) is not None:
                            materials[ID].viscoelastic[k] = [float(tmp)
                                                     for tmp in ve_mat.get(k)]

    materials = OrderedDict(sorted(materials.items()))

    return materials


def find_material(materials, attr, value):
    """
    Find a material object based on a given attribute and value.

    Parameters
    ----------
    materials : OrderedDict
        Dictionary of material objects keyed by material ID, sorted by ID.
    attr : str
        Attribute of the material object to search for.
    value : str
        Value of the attribute to search for.

    Returns
    -------
    Material
        Material object if found, otherwise None.
    """
    return next((x for x in materials.values() if getattr(x, attr).lower() == value.lower()), None)
