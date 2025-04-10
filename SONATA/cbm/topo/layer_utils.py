# -*- coding: utf-8 -*-
"""
Created on Tue Jan 02 10:52:53 2018

@author: TPflumm
"""


def get_layer(lid, SegmentLst):
    if lid < 0:
        print("ERROR:\t lid<0 -> refering to a web and not a layer or a segment")
        return None

    else:
        segid = int(lid / 1000)
        for x in SegmentLst[segid].LayerLst:
            if x.ID == lid:
                # print "i found the lid:", x.ID
                break

            else:
                x = None

        if x == None:
            print("ERROR:\t the lid:", lid, " was not found")

        return x