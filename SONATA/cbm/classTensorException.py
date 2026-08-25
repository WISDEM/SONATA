# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:58:55 2018

@author: TPflumm
"""


class TensorException(Exception):
    def __init__(self, typ=None):
        Exception.__init__(self)
        self.typ = typ

    def __str__(self):
        return "the tensor has to be of type np.array with the size (3,3)"
