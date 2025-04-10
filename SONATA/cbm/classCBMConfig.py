# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:53:09 2018

@author: TPflumm
"""
# Core Library modules
import os
from collections import OrderedDict

# Third party modules
import numpy as np
import yaml

# First party modules
from SONATA.classMaterial import find_material, read_materials
from SONATA.anbax.classANBAXConfig import ANBAXConfig

if __name__ == "__main__":
    os.chdir("/media/gu32kij/work/TPflumm/SONATA")


class CBMConfig(object):
    """
    Second Generation Configuration Class for the SONATA_CBM Disciplin  
    
    Attributes
    ----------
    filename : str
        filename, when config filename is given.
        
    setup: dict
        contains the following fields: datasource, material_db, 
        radial_station, Theta, scale_factor, BalanceWeight and mesh_resolution
    
    webs: OrderedDict
        Ordered Dictionary that contains subdictionaries of with Pos1 and Pos2 
        keys
    
    segments: OrderedDict
         Ordered Dictionary that contains subdictionaries with Corematerial, 
         Layup and Layup_names as keys
    
    bw: dict
    
    flags: dict

    vabs_cfg: VABSConfig
    anbax_cfg: ANBAXConfig

    
    """

    __slots__ = ("filename", "setup", "webs", "segments", "bw", "flags", "vabs_cfg", "anbax_cfg")

    def __init__(self, inputdata=None, materials=None):
        self.setup, self.webs, self.segments, self.bw = {}, {}, {}, {}
        self.filename = ""

        if isinstance(inputdata, dict):
            yml = inputdata
            self.read_yaml_cbm(yml, materials)
        else:
            print("Input data is not a dictionary. Check yaml input file.")

        self.anbax_cfg = ANBAXConfig()
        self.flags = {"mesh_core": True}

    def read_yaml_cbm(self, yml, materials):
        """ read the yml dictionary of a cross section and assign class
        attributes to this configuration object
        
        Parameters:
        ----------
        yml : dict
           dictionary of the yaml cross section beam model (cbm) input
        """

        # Setup:
        self.setup = {}
        self.setup["datasource"] = None
        self.setup["material_db"] = None
        self.setup["radial_station"] = None
        self.setup["Theta"] = None  # GET FROM BLADE DEFINITION
        self.setup["scale_factor"] = 1
        self.setup["BalanceWeight"] = False
        self.setup["mesh_resolution"] = 280  # GET FROM function

        if yml.get("mesh_resolution"):
            self.setup["mesh_resolution"] = yml.get("mesh_resolution")

        # Webs:
        foo = {}
        if yml.get("webs"):
            for item in yml.get("webs"):
                w = {}
                w["id"] = item.get("id")
                w["Pos1"] = item["position"][0]
                w["Pos2"] = item["position"][1]
                w["curvature"] = item.get("curvature")
                foo[item.get("id")] = w
            self.webs = OrderedDict(sorted(foo.items(), key=lambda x: x[1]["id"]))

        # Segments
        self.segments = OrderedDict()
        for i, s in enumerate(yml.get("segments")):
            d = {}
            key = s.get("id")
            if s.get("filler") == None:
                d["CoreMaterial"] = 0

            elif isinstance(s.get("filler"), int):
                d["CoreMaterial"] = materials[s.get("filler")].id

            else:
                d["CoreMaterial"] = find_material(materials, "name", s.get("filler")).id

            layerlst = s.get("layup")
            if layerlst and all(isinstance(l, list) for l in layerlst):
                layerlst = s.get("layup")
                d["Layup_names"] = np.asarray(layerlst)[:, 5].tolist()
                d["Layup"] = np.asarray(layerlst)[:, :5].astype(float)               
            elif layerlst and all(isinstance(l, dict) and l for l in layerlst):
                d["Layup"] = np.asarray([[l.get("start"), l.get("end"), l.get("thickness"), l.get("orientation"), find_material(materials, "name", l.get("material_name")).id] for l in layerlst])
                d["Layup_names"] = [l.get("name") for l in layerlst]

            else:
                d["Layup"] = np.empty((0, 0))
                d["Layup_names"] = np.empty((0, 0))

            self.segments[key] = d

        # BalanceWeight
        if yml.get("trim_mass"):
            self.bw = yml.get("trim_mass")  #
            self.setup["BalanceWeight"] = True


if __name__ == "__main__":

    # classic configuration file:
    os.chdir("/media/gu32kij/work/TPflumm/SONATA")
    # fname = 'jobs/VariSpeed/advanced/sec_config.yml'
    fname = "examples/sec_config.yml"
    cfg = CBMConfig(fname)

    with open("jobs/VariSpeed/UH-60A_adv.yml", "r") as myfile:
        inputs = myfile.read()

    yml = yaml.load(inputs)
    materials = read_materials(yml.get("materials"))
    yml = yml.get("components").get("blade").get("2d_fem").get("sections")[0]
    wt_cfg = CBMConfig(yml, materials)
