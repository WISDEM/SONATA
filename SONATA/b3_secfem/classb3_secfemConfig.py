"""
Created on Wed Apr 01 14:15:44 2020

@author: RFeil
"""
# Third party modules
import numpy as np


class b3_secfemConfig(object):

    """
    this class contains the Configuration for a b3_secfem

    Attributes:
    ----------

    F : nparray([F1, F2, F3]])
        F1 is the sectional aial force, F2 and F3 are the sectional transverse
        shear forces along x2 and x3 respectively

    M : nparray([M1, M2, M3]])
        M1 is the sectional torque, M2 is the setional bending moment around x2
        and M3 is the sectional bending moment around x3.


    """

    def __init__(self, **kw):
        self.recover_flag = 0
        self.ref_sys = "global"
        self.voigt_convention = "b3_secfem"

        if "recover_flag" in kw:
            self.recover_flag = kw["recover_flag"]

        self.F = np.array([0, 0, 0])  # Timoshenko_flag == 0 -> F = [0]
        self.M = np.array([0, 0, 0])



if __name__ == "__main__":
    test = b3_secfemConfig(recover_flag=1)
