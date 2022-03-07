#!/usr/bin/env python

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

from optparse import OptionParser
from collections import OrderedDict
from operator import add
from linecache import getline

from datetime import date
today = date.today()

gc.collect()

folder_path="/media/nenciric/Ricky2020/simulations/"
systems=["dibucaine","etidocaine","SMS","TPP"]

#class by Joe Melcr
#upgraded by Anne
#from https://github.com/akiirik/Databank/blob/main/Scripts/BuildDatabank/OrderParameter.py
class OrderParameter:
    """
    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """
    def __init__(self, name, resname, atom_A_name, atom_B_name, *args):  #removed name, resname from arguments
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        #        self.name = name             # name of the order parameter, a label
        self.resname = resname       # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name

        
        self.name = name
        #M_atom_A_name + " " + M_atom_B_name # generic name of atom A and atom B
        for field in self.__dict__:
            if not isinstance(field, str):
                warnings.warn("provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.")#.format(field)
            else:
                if not field.strip():
                    raise RuntimeError("provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.")#.format(field)
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 0:
            self.avg = None
            self.std = None
            self.stem = None
        elif len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
            self.stem = None
        else:
            warnings.warn("Number of optional positional arguments is {len}, not 2 or 0. Args: {args}\nWrong file format?")
        self.traj = []  # for storing OPs
        self.selection = []
        

    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        # print(atoms)
        bond_len_max=1.5  # in A, max distance between atoms for reasonable OP calculation
        bond_len_max_sq=bond_len_max**2
        
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()

        if d2>bond_len_max_sq:
            at1=atoms[0].name
            at2=atoms[1].name
            resnr=atoms[0].resid
            d=math.sqrt(d2)
            warnings.warn("Atomic distance for atoms \
            {at1} and {at2} in residue no. {resnr} is suspiciously \
            long: {d}!\nPBC removed???")
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        return S
        


    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))
    @property
    def get_avg_std_stem_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        std=np.std(self.traj)
        # convert to numpy array
        return (np.mean(self.traj),std,std/np.sqrt(len(self.traj)-1))
    @property
    def get_op_res(self):
        """
        Provides average and stddev of all OPs in self.traj
        """

        # convert to numpy array
        return self.traj

def go_through_simulation(path,name,system):
    topology = path+name+"/"+name+".gro"
    trajectory = path+name+"/"+name+".xtc"
    top_file = path+name+"/"+name+".top"
    topology_tpr = path+name+"/"+name+".tpr"
    readme = path+name+ "/README.yaml"
    today = str(date.today())
    
    u = mda.Universe(topology,trajectory)
    
    if not os.path.isfile(readme):
        sim={}
    else:
        with open(readme) as yaml_file:
            sim = yaml.load(yaml_file, Loader=yaml.FullLoader)
       
    begin_time=u.trajectory.time

    if not  'TRAJECTORY_SIZE' in sim:   
        sim['TRAJECTORY_SIZE'] = os.path.getsize(trajectory)/1048576

    Nframes=len(u.trajectory)
    timestep = u.trajectory.dt
    
    if not 'TIMESTEP' in sim: 
        sim['TIMESTEP'] = timestep
    
    trj_length = Nframes * timestep
    
    if not 'TRJLENGTH' in sim:
        sim['TRJLENGTH'] = trj_length
        
    if not 'TRJBEGIN' in sim:
        sim['TRJBEGIN'] = begin_time
        
    if not 'TRJEND' in sim:
        sim['TRJEND'] = begin_time+trj_length
        
    if not 'BINDINGEQ' in sim:
        sim['BINDINGEQ'] = input("Biding of {} eqilibrated after [ps] \n".format(output))

        
    if not 'TEMPERATURE' in sim:
        #get temperature from tpr; taken from AddData.py by Anne Kiirikki
        file1 =  'temporary_tpr.txt'

        print("Exporting information with gmx dump")  

        os.system('echo System | gmx dump -s '+ topology_tpr + ' > '+file1)

        with open(file1, 'rt') as tpr_info:
            for line in tpr_info:
                if 'ref-t' in line:
                    sim['TEMPERATURE']=float(line.split()[1])
  
    
    if not 'COMPOSITION' in sim:
        sim["COMPOSITION"]={}
        with open(top_file,"r") as f:
            molecules_list = False
            for line in f.readlines():
                if molecules_list:
                    if not line.startswith(";"):
                        items = line.split()
                        if len(items)==2:
                            sim["COMPOSITION"][items[0]]=items[1]
                elif line.startswith("[ molecules ]"):
                    molecules_list = True

        


    with open(readme, 'w') as f:
        yaml.dump(sim,f, sort_keys=False)
        
        
    with open(system +"_"+ today+ '.out', 'a') as f:
        try:
            f.write("{:80} {:>5} {:>6} {:>4} {:>4} ".format(name,int(timestep),sim['TRJLENGTH']/1000,int(int(sim['BINDINGEQ'])/1000),sim['TEMPERATURE']))
        except:
            f.write("{:80} {:>5} {:>6} {:>4} {:>4} ".format(name,int(timestep),sim['TRJLENGTH']/1000,sim['BINDINGEQ'],sim['TEMPERATURE']))

        for key in sim["COMPOSITION"]:
            f.write("{:>5}: {:>6}, ".format(key, sim["COMPOSITION"][key]))
        f.write("\n")
        
    print("great success!!!")
    
def initialize_output(system,folder_path):
    today = str(date.today())
    with open(system +"_"+ today+ '.out', 'w') as f:
        f.write("# List of simulations for: {} \n".format(system))
        f.write("# Path: {} \n".format(folder_path))
        f.write("# Create on: {} \n \n".format(today))
        f.write("#     1 - System \n")
        f.write("#     2 - Saving frequency [ps]  \n")
        f.write("#     3 - Total length [ns] \n")
        f.write("#     4 - Binding equilibrated after [ns]  \n")
        f.write("#     5 - Temperature [K] \n")
        f.write("#     6 - Composition \n \n")

def initialize_output_analyze(system,folder_path):  
    today = str(date.today())
    column=5
    with open(system + '_analysis_'+ today +'.out', 'w') as f:
            f.write("# List of simulations for: {} \n".format(system))
            f.write("# Path: {} \n".format(folder_path))
            f.write("# Create on: {} \n \n".format(today))
            f.write("#     1 - System \n")
            f.write("#     2 - Start time [ns] \n")
            f.write("#     3 - End time [ns] \n")
            f.write("#     4 - Saving frequency [ps]  \n")
          
    
    ordPars={}
    with open("OP_def.def","r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    items = line.split()
                    ordPars[items[0]] = OrderParameter(*items)
    
    with open(system + '_analysis_'+ today +'.out', 'a') as f:
        f.write("#     {}-{} - ".format(column,column+8))
        for op in ordPars:
            f.write("{} , error {}, ".format(op, op))
        f.write(" \n")
    column+=9
        
    with open(system + '_analysis_'+ today +'.out', 'a') as f:
        f.write("#     {}-{} - Box - Z,  APL \n".format(column,column+1))
    column+=2


            
            
"""Analysis of different stuff - to be developed"""
class AnalysisToolbox:
    def __init__(self,path,name,system,analysis):
        self.topology = path+name+"/"+name+".gro"
        self.topology_tpr = path+name+"/"+name+".tpr"
        self.trajectory = path+name+"/"+name+".xtc"
        self.mol = mda.Universe(self.topology,self.trajectory)
        print(self.mol)
        self.output = name
        self.today = str(date.today())

        
        readme = path+name+ "/README.yaml"
        with open(readme) as yaml_file:
            content = yaml.load(yaml_file, Loader=yaml.FullLoader)
        
        self.readme=content
       
        try:
            self.Nframes=len(self.mol.trajectory)-(int(self.readme['BINDINGEQ'])/int(self.readme['TIMESTEP']))
        except:
            self.Nframes=len(self.mol.trajectory)

        self.system=system
        choose_function = {"box_dimensions": [self.ini_box_dimensions,self.box_dimensions,self.fin_box_dimensions],
                           "order_parameter": [self.ini_order_parameter,self.order_parameter,self.fin_order_parameter]}
        
        print(name)

        if not 'ANALYSIS' in self.readme:
            self.readme['ANALYSIS']={}

        
        self.column=5
        try:
            for analyze in analysis:
                choose_function[analyze][0]()

            
            with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
                f.write(" \n")

            begin_analysis=int(int(self.readme['BINDINGEQ'])/int(self.readme['TIMESTEP']))
            print("going throught the trajectory")
            for self.frame in self.mol.trajectory[begin_analysis:]:
                last_frame=self.frame.time
                #print(last_frame)
                for analyze in analysis:
                    choose_function[analyze][1]()
                
            print("exiting trajectory")
            with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
                f.write("{:75} {:>10} {:>10} {:>10} ".format(
                    name,int(self.readme['BINDINGEQ'])/1000,last_frame/1000,int(self.readme['TIMESTEP'])))
        
            for analyze in analysis:
                choose_function[analyze][2]()
        except:
            print("some trouble")
            with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
                f.write("{:75}  Trouble ".format(name))
            try:
                with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
                    f.write("eqil time: {}, total time: {} \n".format(self.readme['BINDINGEQ'],self.readme['TRJLENGTH']))
            except:
                pass
            
            try:
                os.system('rm ' + self.output + '.xtc' ) 
            except:
                pass
        
             
        with open(readme, 'w') as f:
            yaml.dump(self.readme,f, sort_keys=False)
                    
        

    def ini_box_dimensions(self):
        """At the moment assumes 100 lipids per leaflet"""

        self.box_sizes=[]
        

        
        
    def ini_order_parameter(self):
        
        """
        parses input file with Order Parameter definitions
        file format is as follows:
        OP_name    resname    atom1    atom2  +extra: OP_mean  OP_std
        (flexible cols)
        fname : string
        input file name
        returns : dictionary
        with OrderParameters class instances
        """
        
        print("Make molecules whole in the trajectory")
        os.system('echo System | gmx trjconv -f ' + self.trajectory + ' -s ' + self.topology_tpr + ' -o ' + self.output + ' -pbc mol ' ) # -b ' + str(EQtime))
        self.mol = mda.Universe(self.topology,self.output+'.xtc')
        self.Nframes=len(self.mol.trajectory)-(int(self.readme['BINDINGEQ'])/int(self.readme['TIMESTEP']))
        
        self.ordPars = {}

        with open("OP_def.def","r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    items = line.split()
                    self.ordPars[items[0]] = OrderParameter(*items)
        
        
        
    
        # make atom selections for each OP and store it as its attribute for later use in trajectory
        for self.op in self.ordPars.values():
            # selection = pairs of atoms, split-by residues
            selection = self.mol.select_atoms("resname {rnm} and name {atA} {atB}".format(
                                    rnm=self.op.resname, atA=self.op.atAname, atB=self.op.atBname)
                                    ).atoms.split("residue")

            for res in selection:
                # check if we have only 2 atoms (A & B) selected
                if res.n_atoms != 2:
                    for atom in res.atoms:
                        atA=self.op.atAname
                        atB=self.op.atBname
                        nat=res.n_atoms
                        warnings.warn("Selection >> name {atA} {atB} << \
                           contains {nat} atoms, but should contain exactly 2!")
            self.op.selection = selection

        
            self.Nres=len(self.op.selection)
            self.op.traj= [0]*self.Nres

    
    def box_dimensions(self):
        current_box_z = self.frame.dimensions[2]/10
        current_box_x = self.frame.dimensions[0]/10
          
        current_time = self.frame.time
           
        self.box_sizes.append([current_time,current_box_z,current_box_x])

        
        
    def order_parameter(self):
        for self.op in self.ordPars.values():     
            for i in range(0,self.Nres):
                residue=self.op.selection[i]
                S = self.op.calc_OP(residue)

                self.op.traj[i]=self.op.traj[i]+S/self.Nframes

                
    def fin_order_parameter(self):
        if not 'ORDER_PARAMETER' in self.readme['ANALYSIS']:
                self.readme['ANALYSIS']['ORDER_PARAMETER']={}
        with open(self.output+ "_order_parameter_" + str(self.today),"w") as f:
            f.write("Order parameter analysis \n")
            f.write("Analyzed from {} to {} ns. With saving frequency {} ps.".format(int(self.readme['BINDINGEQ'])/1000,int(self.readme['TRJLENGTH'])/1000,int(self.readme['TIMESTEP'])))
            f.write("# OP_name    resname    atom1    atom2    OP_mean   OP_stddev   OP_err.est. \n\
#--------------------------------------------------------------------------------------------\n" )
            for self.op in self.ordPars.values():
                (self.op.mean, self.op.std, self.op.stem) = self.op.get_avg_std_stem_OP
                f.write( "{:12s} {:10s} {:8s} {:7s} {: 2.5f}  {: 2.5f}    {: 2.5f} \n".format(
                         self.op.name, self.op.resname, self.op.atAname, self.op.atBname,
                         self.op.mean, self.op.std, self.op.stem)
                       )
                with open(self.system + '_analysis_'+ self.today +'.out', 'a') as g:
                    g.write("{: 2.5f} {: 2.5f}  ".format(
                        self.op.mean, self.op.stem))
                try:
                    self.readme['ANALYSIS']['ORDER_PARAMETER'][self.op.name]=float(self.op.mean)
                    self.readme['ANALYSIS']['ORDER_PARAMETER'][self.op.name+"_error"]=float(self.op.stem)
                except:
                    print("problem with writing OPs to readme")
        os.system('rm ' + self.output + '.xtc' ) 
        
        
        


    def fin_box_dimensions(self):
        average_dimensions=np.average(self.box_sizes,axis=0)
        with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
            f.write("{: 3.2f} {: 3.2f}  ".format( average_dimensions[1], average_dimensions[2]**2/100))
        with open(self.output+'_boxSizes_' +self.today+ '.out', 'w') as f:
            np.savetxt(f, self.box_sizes,fmt='%8.4f  %.8f %.8f')
            
        if not 'BOX_DIMENSIONS' in self.readme['ANALYSIS']:
            self.readme['ANALYSIS']['BOX_DIMENSIONS']={}
        

        self.readme['ANALYSIS']['BOX_DIMENSIONS']['REPEAT_DISTANCE']=float(average_dimensions[1])
        self.readme['ANALYSIS']['BOX_DIMENSIONS']['AREA_PER_LIPID']=float(average_dimensions[2]**2/100) 
       

        
def box_dimensions(topology,trajectory,output,system):
    
    u = mda.Universe(topology,trajectory)
    box_sizes=[]
    
    """Go through the simulation - do I want this???, gets APL as a function of time"""
    for ts in u.trajectory:
        #count the index of the frame, numbered from 0, 
        frame = ts.frame  

            
        #reads the dimension (in Angstroms)
        current_box_z = ts.dimensions[2]
        current_box_x = ts.dimensions[0]
          
        last_time = ts.time

            
        box_sizes.append([last_time,current_box_x/10,current_box_z/10])

        
    with open(output+'_boxSizes.out', 'w') as f:
        np.savetxt(f, box_sizes,fmt='%8.4f  %.8f %.8f')
        

"""Go through all simulations and calculate OP and box dimentions"""
for file in os.listdir(folder_path):
    for system in systems:
        initialize_output_analyze(system,folder_path)


for file in os.listdir(folder_path):
    input_corr_file = folder_path+os.fsdecode(file)
    for system in systems:
        if fnmatch.fnmatch(os.fsdecode(file), "*"+system+"*"):
            AnalysisToolbox(folder_path,os.fsdecode(file),system,["order_parameter","box_dimensions"])

