import os
import yaml
import numpy as np
import fnmatch
from datetime import date
today = date.today()


concentrations=["30mM","70mM","140mM","280mM","420mM","1968mM"]
folder_path="/media/nenciric/Ricky20201/simulations/"
waters=["5200w","5985w","7581w","8800w","9975w","11252w","19900w","20000w","177600w"]


for conc in concentrations:
    with open(conc+"_SWISS.dat","w") as f:
        f.write("# Created by R.Nencini on: {} \n \n".format(today))
        f.write("#1 System \n#2 Real system conc. = eti/waters * 55.5  \n#3-#10 alpha1, alp1 err - beta2, beta2 err \n#11 repeat distance \n \n")

for conc in concentrations:
    for water in waters:
        for file in os.listdir(folder_path):
            input_file=folder_path+os.fsdecode(file)
            if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
                if not fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
                    if fnmatch.fnmatch(os.fsdecode(file),"*"+conc+"*"):
                        if fnmatch.fnmatch(os.fsdecode(file),"*"+water+"*"):
                            #print(os.fsdecode(file))
                            readme=folder_path+os.fsdecode(file)+"/README.yaml"
                            with open(readme) as yaml_file:
                                content=yaml.load(yaml_file, Loader=yaml.FullLoader)
                            if "ORDER_PARAMETER" in content["ANALYSIS"]:                                
                                print("{} {:2.1f} {:.4f} {:.1f}".format(
                                    os.fsdecode(file), int(content["COMPOSITION"]["etidocaine"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                    content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))
                                
                                with open(conc+"_SWISS.dat","a") as f:
                                    f.write("{} {:2.1f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}  {:.1f} \n".format(
                                    os.fsdecode(file), int(content["COMPOSITION"]["etidocaine"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                    content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1_error"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2_error"],
                                    content["ANALYSIS"]["ORDER_PARAMETER"]["beta1"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta1_error"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta2"],
                                    content["ANALYSIS"]["ORDER_PARAMETER"]["beta2_error"],
                                    content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))


for conc in concentrations:
    with open(conc+"_PARAMCHEM.dat","w") as f:
        f.write("# Created by R.Nencini on: {} \n \n".format(today))
        f.write("#1 System \n#2 Real system conc. = eti/waters * 55.5  \n#3-#10 alpha1, alp1 err - beta2, beta2 err \n#11 repeat distance \n \n")

for conc in concentrations:
    for water in waters:
        for file in os.listdir(folder_path):
            input_file=folder_path+os.fsdecode(file)
            if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
                if fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
                    if not fnmatch.fnmatch(os.fsdecode(file),"*ECC*"):
                        if fnmatch.fnmatch(os.fsdecode(file),"*"+conc+"*"):
                            if fnmatch.fnmatch(os.fsdecode(file),"*"+water+"*"):
                                #print(os.fsdecode(file))
                                readme=folder_path+os.fsdecode(file)+"/README.yaml"
                                with open(readme) as yaml_file:
                                    content=yaml.load(yaml_file, Loader=yaml.FullLoader)
                                if "ORDER_PARAMETER" in content["ANALYSIS"]:
                                    print("{} {:2.1f} {:.4f} {:.1f}".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))

                                    with open(conc+"_PARAMCHEM.dat","a") as f:
                                        f.write("{} {:2.1f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}  {:.1f} \n".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1_error"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2_error"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta1"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta1_error"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta2"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta2_error"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))


for conc in concentrations:
    with open(conc+"_ECC.dat","w") as f:
        f.write("# Created by R.Nencini on: {} \n \n".format(today))
        f.write("#1 System \n#2 Real system conc. = eti/waters * 55.5  \n#3-#10 alpha1, alp1 err - beta2, beta2 err \n#11 repeat distance \n \n")

for conc in concentrations:
    for water in waters:
        for file in os.listdir(folder_path):
            input_file=folder_path+os.fsdecode(file)
            if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
                if fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
                    if fnmatch.fnmatch(os.fsdecode(file),"*ECC*"):
                        if fnmatch.fnmatch(os.fsdecode(file),"*"+conc+"*"):
                            #print(os.fsdecode(file))
                            if fnmatch.fnmatch(os.fsdecode(file),"*"+water+"*"):
                                print(os.fsdecode(file))
                                readme=folder_path+os.fsdecode(file)+"/README.yaml"
                                with open(readme) as yaml_file:
                                    content=yaml.load(yaml_file, Loader=yaml.FullLoader)
                                if "ORDER_PARAMETER" in content["ANALYSIS"]:
                                    print("{} {:2.1f} {:.4f} {:.1f}".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))

                                    with open(conc+"_ECC.dat","a") as f:
                                        f.write("{} {:2.1f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}  {:.1f} \n".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1_error"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2_error"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta1"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta1_error"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta2"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta2_error"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))

folder_path="/media/nenciric/Ricky2020/2020/simulations/"
for conc in concentrations:
    with open(conc+"_SWISS_disk2.dat","w") as f:
        f.write("# Created by R.Nencini on: {} \n \n".format(today))
        f.write("#1 System \n#2 Real system conc. = eti/waters * 55.5  \n#3-#10 alpha1, alp1 err - beta2, beta2 err \n#11 repeat distance \n \n")

for conc in concentrations:
    for water in waters:
        for file in os.listdir(folder_path):
            input_file=folder_path+os.fsdecode(file)
            if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
                if not fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
                    if fnmatch.fnmatch(os.fsdecode(file),"*"+conc+"*"):
                        if fnmatch.fnmatch(os.fsdecode(file),"*"+water+"*"):
                            #print(os.fsdecode(file))
                            readme=folder_path+os.fsdecode(file)+"/README.yaml"
                            with open(readme) as yaml_file:
                                content=yaml.load(yaml_file, Loader=yaml.FullLoader)
                            if "ORDER_PARAMETER" in content["ANALYSIS"]:
                                print("{} {:2.1f} {:.4f} {:.1f}".format(
                                    os.fsdecode(file), int(content["COMPOSITION"]["etidocaine"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                    content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))

                                with open(conc+"_SWISS_disk2.dat","a") as f:
                                    f.write("{} {:2.1f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}  {:.1f} \n".format(
                                    os.fsdecode(file), int(content["COMPOSITION"]["etidocaine"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                    content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1_error"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2_error"],
                                    content["ANALYSIS"]["ORDER_PARAMETER"]["beta1"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta1_error"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta2"],
                                    content["ANALYSIS"]["ORDER_PARAMETER"]["beta2_error"],
                                    content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))




for conc in concentrations:
    with open(conc+"_PARAMCHEM_disk2.dat","w") as f:
        f.write("# Created by R.Nencini on: {} \n \n".format(today))
        f.write("#1 System \n#2 Real system conc. = eti/waters * 55.5  \n#3-#10 alpha1, alp1 err - beta2, beta2 err \n#11 repeat distance \n \n")

for conc in concentrations:
    for water in waters:
        for file in os.listdir(folder_path):
            input_file=folder_path+os.fsdecode(file)
            if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
                if fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
                    if not fnmatch.fnmatch(os.fsdecode(file),"*ECC*"):
                        if fnmatch.fnmatch(os.fsdecode(file),"*"+conc+"*"):
                            if fnmatch.fnmatch(os.fsdecode(file),"*"+water+"*"):
                                print(os.fsdecode(file))
                                readme=folder_path+os.fsdecode(file)+"/README.yaml"
                                with open(readme) as yaml_file:
                                    content=yaml.load(yaml_file, Loader=yaml.FullLoader)
                                if "ORDER_PARAMETER" in content["ANALYSIS"]:
                                    print("{} {:2.1f} {:.4f} {:.1f}".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))

                                    with open(conc+"_PARAMCHEM_disk2.dat","a") as f:
                                        f.write("{} {:2.1f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}  {:.1f} \n".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1_error"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2_error"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta1"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta1_error"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta2"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta2_error"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))



for conc in concentrations:
    with open(conc+"_ECC_disk2.dat","w") as f:
        f.write("# Created by R.Nencini on: {} \n \n".format(today))
        f.write("#1 System \n#2 Real system conc. = eti/waters * 55.5  \n#3-#10 alpha1, alp1 err - beta2, beta2 err \n#11 repeat distance \n \n")

for conc in concentrations:
    for water in waters:
        for file in os.listdir(folder_path):
            input_file=folder_path+os.fsdecode(file)
            if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
                if fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
                    if fnmatch.fnmatch(os.fsdecode(file),"*ECC*"):
                        if fnmatch.fnmatch(os.fsdecode(file),"*"+conc+"*"):
                            #print(os.fsdecode(file))
                            if fnmatch.fnmatch(os.fsdecode(file),"*"+water+"*"):
                                print(os.fsdecode(file))
                                readme=folder_path+os.fsdecode(file)+"/README.yaml"
                                with open(readme) as yaml_file:
                                    content=yaml.load(yaml_file, Loader=yaml.FullLoader)
                                if "ORDER_PARAMETER" in content["ANALYSIS"]:
                                    print("{} {:2.1f} {:.4f} {:.1f}".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))

                                    with open(conc+"_ECC_disk2.dat","a") as f:
                                        f.write("{} {:2.1f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}  {:.1f} \n".format(
                                        os.fsdecode(file), int(content["COMPOSITION"]["ETI"])/int(content["COMPOSITION"]["TIP3"])*55.5*1000, content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["alpha1_error"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2"], content["ANALYSIS"]["ORDER_PARAMETER"]["alpha2_error"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta1"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta1_error"],content["ANALYSIS"]["ORDER_PARAMETER"]["beta2"],
                                        content["ANALYSIS"]["ORDER_PARAMETER"]["beta2_error"],
                                        content["ANALYSIS"]["BOX_DIMENSIONS"]["REPEAT_DISTANCE"]))


print(".... ")
for file in os.listdir(folder_path):
    input_file=folder_path+os.fsdecode(file)
    if fnmatch.fnmatch(os.fsdecode(file),"*etidocaine*"):
        if not  fnmatch.fnmatch(os.fsdecode(file),"*paramchem*"):
            if not fnmatch.fnmatch(os.fsdecode(file),"*ECC*"):
                print(os.fsdecode(file))

