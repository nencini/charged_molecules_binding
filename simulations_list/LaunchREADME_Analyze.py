import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import fnmatch
import yaml
import gc
import math
import warnings
import re
import time

from optparse import OptionParser
from collections import OrderedDict
from operator import add
from linecache import getline


#sys.path.insert(1, '/home/nenciric/Documents/git/charged_molecules_binding/simulations_list/')

import AnalysisToolbox as AT

from datetime import date
today = date.today()

gc.collect()

folder_path="/media/nenciric/Ricky20201/2020/simulations/"
systems=["etidocaine","TPP","SMS","dibucaine"]


#folder_path="/home/ricky/Documents/from_work/MD/simulations/production_run/"
#systems=["testing_density_center"]


folder_path="/DATA/hector/200POPC_PN_NAs150mM/"
systems=["PN_m"]


for file in os.listdir(folder_path):
    input_corr_file = folder_path+os.fsdecode(file)
    for system in systems:
        if fnmatch.fnmatch(os.fsdecode(file), "*"+system+"*"):
            newcomer=AT.AnalysisToolbox(folder_path,os.fsdecode(file),system,["ORDER_PARAMETER","BOX_DIMENSIONS"])
            #newcomer=AT.AnalysisToolbox(folder_path,os.fsdecode(file),"etidocaine",["ORDER_PARAMETER","BOX_DIMENSIONS"])
            newcomer.add_new_folders()
            newcomer.analysis_module()
