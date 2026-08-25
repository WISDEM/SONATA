# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:30:28 2017

@author: TPflumm
"""
# Core Library modules
import math
import os

# Third party modules
import numpy as np

# First party modules
from SONATA.cbm.classTensorException import TensorException

if __name__ == "__main__":
    os.chdir("../..")



class Stress(object):
    """Object to represent the stress tensor"""

    __slots__ = "__tensor"

    def __init__(self, Vec=None):
        self.__tensor = np.zeros((3, 3))

        if Vec is not None:
            self.__tensor = np.zeros((3, 3))
            idx = np.triu_indices_from(self.__tensor)
            self.__tensor[idx] = Vec
            self.__tensor[(idx[1], idx[0])] = Vec

    def __getTensor(self):
        return self.__tensor

    def __setTensor(self, x):
        if isinstance(x, np.ndarray) and x.shape == (3, 3):
            self.__tensor = x
        else:
            raise TensorException

    tensor = property(__getTensor, __setTensor)

    @property
    def sigma11(self):
        return self.__tensor[0, 0]

    @property
    def sigma12(self):
        return self.__tensor[0, 1]

    @property
    def sigma13(self):
        return self.__tensor[0, 2]

    @property
    def sigma22(self):
        return self.__tensor[1, 1]

    @property
    def sigma23(self):
        return self.__tensor[1, 2]

    @property
    def sigma33(self):
        return self.__tensor[2, 2]

    @property
    def sigma_vM(self):
        """Calculates the von_Mises stresses"""
        sigma_vM = math.sqrt(
            self.sigma11 ** 2
            + self.sigma22 ** 2
            + self.sigma33 ** 2
            - self.sigma11 * self.sigma22
            - self.sigma11 * self.sigma33
            - self.sigma22 * self.sigma33
            + 3 * (self.sigma12 ** 2 + self.sigma13 ** 2 + self.sigma23 ** 2)
        )
        return sigma_vM


if __name__ == "__main__":
    s = Stress([11, 12, 13, 22, 23, 33])
    s.tensor = np.array([[11, 3, 2], [23, 23, 2], [2, 1, 2]])
    print(s.sigma_vM)
