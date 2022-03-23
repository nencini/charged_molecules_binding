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
import time

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
        self.path=path
        self.name=name
        
        readme = path+name+ "/README.yaml"
        with open(readme) as yaml_file:
            content = yaml.load(yaml_file, Loader=yaml.FullLoader)
        
        self.readme=content
        self.system=system
        
        self.check_for_latest_files()
        
        choose_function = {"BOX_DIMENSIONS": [self.ini_box_dimensions,self.box_dimensions,self.fin_box_dimensions],
                           "ORDER_PARAMETER": [self.ini_order_parameter,self.order_parameter,self.fin_order_parameter],
                          "binding_coefficient": [self.ini_binding_coefficient,self.binding_coefficient,self.fin_binding_coefficient]}
                                                 
        
        print(name)

        if not 'ANALYSIS' in self.readme:
            self.readme['ANALYSIS']={}

        
        self.column=5
        try:
            #probably for OP, should work without it here
            #self.Nframes=len(self.mol.trajectory)-(int(self.readme['BINDINGEQ'])/int(self.readme['TRAJECTORY']['TIMESTEP']))
        
            """For whathever reason does not work withou the outer loop. IT SHOULD though!!"""
            
            if not 'BINDINGEQ' in self.readme:
                analysis=[]
            else:
                try:
                    int(self.readme['BINDINGEQ'])
                except:
                    analysis=[]
            
            for i in range(0,len(analysis)):
                for analyze in analysis:
                    print("check if analyzed: {}".format(analyze))
                    if analyze in self.readme["ANALYSIS"]:
                        if self.readme["ANALYSIS"][analyze]["FROM_XTC"]==self.readme["FILES"]["xtc"]["MODIFIED"]:                
                            analysis.remove(analyze)
           
            if analysis==[]:
                pass
            else:
                for analyze in analysis:
                    choose_function[analyze][0]()

                begin_analysis=int(int(self.readme['BINDINGEQ'])/int(self.readme["FILES"]['xtc']['TIMESTEP']))
                print("going throught the trajectory")
                self.frames=0
                
                for self.frame in self.mol.trajectory[begin_analysis:]:
                    last_frame=self.frame.time

                    for analyze in analysis:
                        choose_function[analyze][1]()
                    self.frames+=1
                
                print("exiting trajectory, writing output ...")
            
                for analyze in analysis:
                    choose_function[analyze][2]()
                print("exit output")
        except Exception as e:
            print(e)
            print("some trouble")

            with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
                f.write("{:75}  Trouble ".format(name))
            try:
                with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
                    f.write("eqil time: {}, total time: {} \n".format(self.readme['BINDINGEQ'],self.readme["FILES"]['xtc']['LENGTH']))
            except:
                pass
            
            try:
                os.system('rm ' + self.output + '.xtc' ) 
            except:
                pass
        
             
        with open(readme, 'w') as f:
            yaml.dump(self.readme,f, sort_keys=False)
    
    
    ###############
    # MODIFICATIONS
    ###############
    
    def check_for_latest_files(self):
        files_to_consider=["xtc","edr","tpr","top","mdp","ndx","gro","cpt","log"]
        sim=self.readme
        if not "FILES" in sim:
            sim["FILES"]={}
        for file in files_to_consider:
            if not file in sim["FILES"]:
                sim["FILES"][file]={}
    
            file_adress = self.path+self.name+"/"+self.name+"."+file
            sim["FILES"][file]["NAME"] = self.name+"."+file
    
            try:
                timepre=os.path.getmtime(file_adress)
                file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
                sim["FILES"][file]["SIZE"]=os.path.getsize(file_adress)/1048576
                sim["FILES"][file]["MODIFIED"] = file_mod
            except:
                sim["FILES"][file]["SIZE"]= "none"
                sim["FILES"][file]["MODIFIED"] = "none"
        
        Nframes=len(self.mol.trajectory)
        timestep = self.mol.trajectory.dt
        trj_length = Nframes * timestep
        
        sim["FILES"]["xtc"]['TIMESTEP'] = timestep
        sim["FILES"]['xtc']['LENGTH'] = trj_length
                
                
    ############
    # BINDING COEFFICIENT
    ############
    
    def ini_binding_coefficient(self):
        """
        similar initiotion to density
        collects density data along the z coordinate
        on the range of the first frame box size.
        
        The final result is cut-off to the smallest box size
        """
        
        self.nbin = 200
        self.c = self.mol.select_atoms('resname POPC')      #selection for centering the profile   
        box_z = self.mol.dimensions[2] # takes the size of the first frame and uses it for the density calculation
        self.nbin=int(box_z*3)
        self.bsize=box_z/10   # convert to [nm]
        self.min_z = box_z   # used to search for the min box size to cut off the profile on the edges
        self.bin_width =  2*self.bsize/ self.nbin #     # bin width, used to be d
        self.boxH = self.bsize # self.bsize used instead of boxH
        self.z_coor = np.linspace(-self.bsize,self.bsize,self.nbin+1)[:-1] + self.bin_width/2 # x changed to z_coor
        
        self.fz_water = np.zeros(self.nbin)
        self.fz_lipid = np.zeros(self.nbin)        
        self.fz_particle = np.zeros(self.nbin)
    
        self.water = self.mol.select_atoms('resname TIP3')
        self.lipid = self.mol.select_atoms('resname POPC')
        
        resnames={"etidocaine":"resname ETI","dibucaine":"resname DIB", "SMS": "resname SMS", "TPP":"resname TPA"}
        
        
        self.particleName = resnames[self.system]
        
        self.particle = self.mol.select_atoms(self.particleName)
    
        self.wght_water=np.ones(self.water.atoms.names.shape[0])
        self.wght_lipid=np.ones(self.lipid.atoms.names.shape[0])   
        self.wght_particle=np.ones(self.particle.atoms.names.shape[0])

        
    def binding_coefficient(self):
        self.current_dimentions = self.mol.trajectory.ts.dimensions
        self.box_z = self.current_dimentions[2]
        if self.box_z/10<self.min_z:
            self.min_z=self.box_z/10
                

        crds = self.mol.atoms.positions
        ctom = self.c.atoms.center_of_mass()[2]
        crds[:,2] += self.box_z/2 - ctom
        self.mol.atoms.positions = crds
        self.mol.atoms.pack_into_box()

        self.fz_water += self.density_cal(self.water,self.wght_water)
        self.fz_lipid +=  self.density_cal(self.lipid,self.wght_lipid)
        self.fz_particle +=  self.density_cal(self.particle,self.wght_particle)
        
    def fin_binding_coefficient(self):
        self.fz_water /= self.frames
        self.fz_lipid /=  self.frames
        self.fz_particle /=  self.frames
        
        """Get the indexes of the final density data where all the time steps contribute
        In other words, take the coordinates of the smalest box from the simulation"""
        final_FF_start=int(np.round(self.nbin/2-self.min_z/self.bin_width/2))+1
        final_FF_end=int(np.round(self.nbin/2+self.min_z/self.bin_width/2))-1
        
        data = np.vstack((self.z_coor,self.fz_water,self.fz_lipid,self.fz_particle)).transpose()
        density = data[final_FF_start+1:final_FF_end-1,:]
        
    
        
        minus_water=0
        minus_lipid=0
        plus_water=0
        plus_lipid=0
        count_minus_water=0
        count_minus_lipid=0
        count_plus_water=0
        count_plus_lipid=0

        for density_slice in density:
            if density_slice[0]<0:
                if density_slice[2]<density_slice[1]:
                    minus_water+=density_slice[3]/density[0,3]
                    count_minus_water+=1
                else:
                    minus_lipid+=density_slice[3]/density[0,3]
                    count_minus_lipid+=1
            else:
                if density_slice[2]<density_slice[1]:
                    plus_water+=density_slice[3]/density[0,3]
                    count_plus_water+=1
                else:
                    plus_lipid+=density_slice[3]/density[0,3]
                    count_plus_lipid+=1

                    
        if not 'BINDING_COEFFICIENT' in self.readme['ANALYSIS']:
            self.readme['ANALYSIS']['BINDING_COEFFICIENT']={}
        
        self.readme['ANALYSIS']['BINDING_COEFFICIENT']['DENSITY_CROSS']={}
        self.readme['ANALYSIS']['BINDING_COEFFICIENT']['DENSITY_CROSS']['FROM_XTC']=self.readme["FILES"]["xtc"]["MODIFIED"]
        self.readme['ANALYSIS']['BINDING_COEFFICIENT']['DENSITY_CROSS']['ANALYZED']=self.today
        self.readme['ANALYSIS']['BINDING_COEFFICIENT']['DENSITY_CROSS']['MINUS_VALUE']=float(minus_lipid/minus_water*count_minus_water/count_minus_lipid)
        self.readme['ANALYSIS']['BINDING_COEFFICIENT']['DENSITY_CROSS']['PLUS_VALUE']=float(plus_lipid/plus_water*count_plus_water/count_plus_lipid)
        
        with open(self.output+"_density.out", 'wb') as f:
            np.savetxt(f, density,fmt='%8.4f  %.8f %.8f %.8f')
      

    ############
    # BOX DIMENSIONS
    ############       
    
    def ini_box_dimensions(self):
        """At the moment assumes 100 lipids per leaflet"""
        print("Box dimensions will be analyzed")

        self.box_sizes=[]
        
    def box_dimensions(self):
        current_box_z = self.frame.dimensions[2]/10
        current_box_x = self.frame.dimensions[0]/10
          
        current_time = self.frame.time
           
        self.box_sizes.append([current_time,current_box_z,current_box_x])

   
    def fin_box_dimensions(self):
        average_dimensions=np.average(self.box_sizes,axis=0)
        with open(self.system + '_analysis_'+ self.today +'.out', 'a') as f:
            f.write("{: 3.2f} {: 3.2f}  ".format( average_dimensions[1], average_dimensions[2]**2/100))
        with open(self.output+'_boxSizes_' +self.today+ '.out', 'w') as f:
            np.savetxt(f, self.box_sizes,fmt='%8.4f  %.8f %.8f')
            
        if not 'BOX_DIMENSIONS' in self.readme['ANALYSIS']:
            self.readme['ANALYSIS']['BOX_DIMENSIONS']={}
        
        self.readme['ANALYSIS']['BOX_DIMENSIONS']['FROM_XTC']=self.readme["FILES"]["xtc"]["MODIFIED"]
        self.readme['ANALYSIS']['BOX_DIMENSIONS']['ANALYZED']=self.today
        self.readme['ANALYSIS']['BOX_DIMENSIONS']['REPEAT_DISTANCE']=float(average_dimensions[1])
        self.readme['ANALYSIS']['BOX_DIMENSIONS']['AREA_PER_LIPID']=float(average_dimensions[2]**2/100) 
       
    ############
    # ORDER PARAMETER
    ############         
        
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
        
        print("Order parameter will be analyzed")
        print("Make molecules whole in the trajectory")
        os.system('echo POPC | gmx trjconv -f ' + self.trajectory + ' -s ' + self.topology_tpr + ' -o lipids.gro -pbc mol -b 0 -e 0 ' ) # -b ' + str(EQtime))
        os.system('echo POPC | gmx trjconv -f ' + self.trajectory + ' -s ' + self.topology_tpr + ' -o ' + self.output + ' -pbc mol ' ) # -b ' + str(EQtime))
        self.mol = mda.Universe('lipids.gro',self.output+'.xtc')
        self.Nframes=len(self.mol.trajectory)-(int(self.readme['BINDINGEQ'])/int(self.readme["FILES"]['xtc']['TIMESTEP']))
        
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
        
        print("Ini OP sucessfuly done")

    

        
        
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
            f.write("Analyzed from {} to {} ns. With saving frequency {} ps.".format(int(self.readme['BINDINGEQ'])/1000,int(self.readme["FILES"]['xtc']['LENGTH'])/1000,int(self.readme["FILES"]['xtc']['TIMESTEP'])))
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
                    self.readme['ANALYSIS']['ORDER_PARAMETER']['FROM_XTC']=self.readme["FILES"]["xtc"]["MODIFIED"]
                    self.readme['ANALYSIS']['ORDER_PARAMETER']['ANALYZED']=self.today
                    self.readme['ANALYSIS']['ORDER_PARAMETER'][self.op.name]=float(self.op.mean)
                    self.readme['ANALYSIS']['ORDER_PARAMETER'][self.op.name+"_error"]=float(self.op.stem)
                except:
                    print("problem with writing OPs to readme")
        os.system('rm ' + self.output + '.xtc' ) 
        os.system('rm lipids.gro' ) 
        
        
    ############
    # DENSITY
    ############      
    def ini_density(self):
        pass
        
    def density(self):
        pass
        
    def fin_density(self):
        pass
        
    ############
    # OTHER FUNCTIONS
    ############   

    def density_cal(self,atoms_for_density,atom_weights):
    
        crds_of_atoms = (atoms_for_density.atoms.positions[:,2] - self.box_z/2)/10                          
        vbin = self.bin_width*np.prod(self.current_dimentions[:2])/100      
        return np.histogram(crds_of_atoms,bins=self.nbin,range=(-self.boxH,self.boxH),weights=atom_weights/vbin)[0]
     

        
        
        
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
#for file in os.listdir(folder_path):
#    for system in systems:
#        initialize_output_analyze(system,folder_path)


#for file in os.listdir(folder_path):
#    input_corr_file = folder_path+os.fsdecode(file)
#    for system in systems:
#        if fnmatch.fnmatch(os.fsdecode(file), "*"+system+"*"):
#            AnalysisToolbox(folder_path,os.fsdecode(file),system,["order_parameter","box_dimensions"])

