#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:59:15 2019

@author: gu32kij
"""
# Core Library modules
import os

# Third party modules
import numpy as np

# First party modules

if __name__ == "__main__":
    os.chdir("../..")



class BeamSectionalProps(object):
    """
    this class stores the beam cross-sectional data and has methods to 
    read the result files from vabs
    
    The coordinate system follows the definition of VABS and DYMORE. 
    For more information see the the vabs and dymore user manual.
    the x/x1 axis is in the direction of the beam. the y/x2 axis points towards
    the leading edge, and z/x3 points accordingly upwards.
    
    The 6x6 sectional stiffness matrix, TS (Timoshenko Stiffness Matrix)
    (1-extension; 2,3-shear, 4-twist; 5,6-bending) relates the sectional axial 
    strain, epsilon1, transverse shearing strains, epsilon2 and epsilon3, 
    twisting curvatures, kappa1 and two bending curvatures, kappa2 and kappa3, 
    to the axial force, F1, transverse shear forces, F2 and F3, twisting 
    moment, M1, and two bending moments, M2 and M3. The relationship between 
    these sectional strains and sectional stress resultants takes the form of a 
    symmetric, 6x6 matrix
    | F1 |   | k11 k12 k13 k14 k15 k16 |   | epsilon1 | 
    | F2 |   | k12 k22 k23 k24 k25 k26 |   | epsilon2 |
    | F3 |   | k13 k23 k33 k34 k35 k36 |   | epsilon3 |
    | M1 | = | k14 k24 k34 k44 k45 k46 | * | kappa1   |
    | M2 |   | k15 k25 k35 k45 k55 k56 |   | kappa2   |
    | M3 |   | k16 k26 k36 k46 k56 k66 |   | kappa3   |

    The 6x6 sectional mass matrix, MM, relates the sectional linear velocities, 
    denoted v1, v2 and v3, and angular velocities, denoted omega1, omega2 and 
    omgea3, to the sectional linear momenta, denoted p1, p2 and p3, and angular
    momenta, denoted h1, h2 and h3. The relationship between these sectional 
    velocities and sectional momenta takes the form of a symmetric, 6x6 matrix:
    | p1 |   | M11 M12 M13 M14 M15 M16 |   | v1     | 
    | p2 |   | M12 M22 M23 M24 M25 M26 |   | v2     |
    | p3 |   | M13 M23 M33 M34 M35 M36 |   | v3     |
    | h1 | = | M14 M24 M34 M44 M45 M46 | * | omega1 |
    | h2 |   | M15 M25 M35 M45 M55 M56 |   | omega2 |
    | h3 |   | M16 M26 M36 M46 M56 M66 |   | omega3 |
    
    Due to the nature of the problem, many of these coefficients vanish, and 
    the remaining entries are written as:
    | p1 |   | m00        0        0         0  m00*Xm3  −m00*Xm2 |   | v1    | 
    | p2 |   | 0        m00        0  −m00*Xm3        0         0 |   | v2    |
    | p3 |   | 0          0      m00   m00*Xm2        0         0 |   | v3    |
    | h1 | = | 0   −m00*Xm3  m00*Xm2       m11        0         0 | * | omega1|
    | h2 |   | m00*Xm3    0        0        0       m22      −m23 |   | omega2|
    | h3 |   | −m00*Xm2   0        0        0      −m23       m33 |   | omega3|
    
    where m00 is the sectional mass per unit span and Xm2 and Xm3 the 
    coordinates of the sectional center of mass. m22 and m33 are the section 
    mass moments of inertia per unit span with respect to unit vectors i2 and 
    i3, respectively, and m23 the sectional cross-product of inertia per unit 
    span. The polar moment of inertia per unit span, denoted m11, is given by
    m11 = m22 + m33.


    Attributes
    ----------
    TS : ndarray
        Timoshenko Stiffness Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
    MM : ndarray
        6x6 Mass Matrix    
    Xg : ndarray
        Geometric Center of the Cross Section
    Xt : ndarray
        Neutral Axes (or Tension Center) of the Cross Section
    Xs : ndarray
        The Generalized Shear Center of the Cross Section
    VS: ndarray
        Vlasov Stiffness Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)
    Ag: ndarray
        Trapeze Effects Related Matrices
    Bk: ndarray
        Trapeze Effects Related Matrices
    Ck: ndarray
        Trapeze Effects Related Matrices
    Dk: ndarray
        Trapeze Effects Related Matrices   


    Properties
    ----------
    CS : ndarray
        Classical Stiffness Matrix (1-extension; 2-twist; 3,4-bending)
    CF : ndarray
        Classical Flexibility Matrix (1-extension; 2-twist; 3,4-bending)
    TF : ndarray
        Timoshenko Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
    MMatMC : ndarray
        6x6 Mass Matrix at Mass Center
    Xm : ndarray
        Center of Gravity (Mass Center)
    m00 : float
        mass per unit span
    m11 : float
        polar moment of inertia per unit span
    m22 : float
        mass moments of inertia per unit span with respect to unit vectors i2 
    m33 : float
        mass moments of inertia per unit span with respect to unit vectors i3 
    VF: ndarray
        Vlasov Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)

    ToDo
    ----------
    - Implement function to rotate the Vlasov Stiffness matrix and the Trapez 
        Effect related terms
    - Make one unit convertion method that converts the complete instance
    """    
    
    __slots__ = ('TS', 'MM', 'Xg', 'Xt', 'Xs', 'Delta', 'PIA', 'VS', 'Ag', \
                 'Bk', 'Ck', 'Dk', 'ELE', 'U')
    

    def __init__(self, fname=None):
        self.TS = None
        self.MM = None
        self.Xg = None
        self.Xt = None
        self.Xs = None
        self.Delta = None
        self.PIA = None
        self.VS = None
        self.Ag = None
        self.Bk = None
        self.Ck = None
        self.Dk = None

        self.ELE = None
        self.U = None

        if fname:
            self.read_vabs_K(fname)

    @property
    def m00(self):
        return self.MM[0, 0]

    @property
    def m11(self):
        return self.MM[3, 3]

    @property
    def m22(self):
        return self.MM[4, 4]

    @property
    def m33(self):
        return self.MM[5, 5]

    @property
    def m23(self):
        return -self.MM[4, 5]

    @property
    def Xm(self):
        return np.array([-self.MM[0, 5] / self.m00, self.MM[0, 4] / self.m00])

    @property
    def rg(self):
        return np.sqrt(self.m11 / self.m00)

    @property
    def MMatMC(self):
        mask = np.zeros(self.MM.shape, dtype=bool)
        np.fill_diagonal(mask, True)
        mask[4, 5] = True
        mask[5, 4] = True
        return self.MM * mask

    #    @property
    #    def CS(self):
    #        tmp = np.delete(self.TS,[1,2],axis=0)
    #        return np.delete(tmp,[1,2],axis=1)

    @property
    def TF(self):
        return np.linalg.inv(self.TS)

    #    @property
    #    def CF(self):
    #        return np.linalg.inv(self.CS)

    @property
    def VF(self):
        return np.linalg.inv(self.VS)

# ======================================================
#       MAIN
# ======================================================
if __name__ == "__main__":
    filename = "cbm_noname_2019-04-02_145252576501.vab.K"
    bp1 = BeamSectionalProps(filename)
    bp2 = bp1.rotate(np.pi / 8, copy=True)
