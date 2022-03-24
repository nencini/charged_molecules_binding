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


sys.path.insert(1, '/home/nenciric/Documents/git/charged_molecules_binding/simulations_list/')
import AnalysisToolbox as AT

from datetime import date
today = date.today()

gc.collect()

folder_path="/media/nenciric/Ricky2020/2020/simulations/"
systems=["etidocaine","TPP","SMS","dibucaine"]


for file in os.listdir(folder_path):
    input_corr_file = folder_path+os.fsdecode(file)
    for system in systems:
        if fnmatch.fnmatch(os.fsdecode(file), "*"+system+"*"):
            AT.AnalysisToolbox(folder_path,os.fsdecode(file),system,["ORDER_PARAMETER","BOX_DIMENSIONS"])


folder_path="/media/nenciric/Ricky20201/simulations/"
systems=["etidocaine","TPP","SMS","dibucaine"]


for file in os.listdir(folder_path):
    input_corr_file = folder_path+os.fsdecode(file)
    for system in systems:
        if fnmatch.fnmatch(os.fsdecode(file), "*"+system+"*"):
            AT.AnalysisToolbox(folder_path,os.fsdecode(file),system,["ORDER_PARAMETER","BOX_DIMENSIONS"])
                                        
