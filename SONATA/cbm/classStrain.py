# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:25:49 2017

@author: TPflumm
"""
# Core Library modules
import os

# Third party modules
import numpy as np

# First party modules
from SONATA.cbm.classTensorException import TensorException

if __name__ == "__main__":
    os.chdir("../..")



class Strain(object):
    """object to represent the strain tensor"""

    __slots__ = "__tensor"

    def __init__(self, Vec=None):
        self.__tensor = np.zeros((3, 3))

        if Vec is not None:
            self.__tensor = np.zeros((3, 3))
            idx = np.triu_indices_from(self.__tensor)
            self.__tensor[idx] = Vec
            self.__tensor[(idx[1], idx[0])] = Vec
            self.__tensor[0, 1] = self.__tensor[0, 1] / 2
            self.__tensor[0, 2] = self.__tensor[0, 2] / 2
            self.__tensor[1, 2] = self.__tensor[1, 2] / 2

    def __getTensor(self):
        return self.__tensor

    def __setTensor(self, x):
        if isinstance(x, np.ndarray) and x.shape == (3, 3):
            self.__tensor = x
        else:
            raise TensorException

    tensor = property(__getTensor, __setTensor)

    @property  # better method with the use of decorators
    def epsilon11(self):
        return self.tensor[0, 0]

    @property
    def epsilon12(self):
        return self.tensor[0, 1]

    @property
    def gamma12(self):  # 2*epsilon12
        return 2 * self.tensor[0, 1]

    @property
    def epsilon13(self):
        return self.tensor[0, 2]

    @property
    def gamma13(self):  # 2*epsilon13
        return 2 * self.tensor[0, 2]

    @property
    def epsilon22(self):
        return self.tensor[1, 1]

    @property
    def epsilon23(self):
        return self.tensor[1, 2]

    @property
    def gamma23(self):  # 2*epsilon23
        return 2 * self.tensor[1, 2]

    @property
    def epsilon33(self):
        return self.tensor[2, 2]


if __name__ == "__main__":
    s = Strain([11, 12, 13, 22, 23, 33])
    s.tensor = np.array([[11, 3, 2], [23, 23, 2], [2, 1, 2]])
