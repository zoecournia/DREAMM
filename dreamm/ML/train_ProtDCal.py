# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:24:27 2019

@author: alexis
"""

import os, sys
from moleculekit.molecule import Molecule

file = sys.argv[1]

filename = sys.argv[2]

mol = Molecule(filename)

#ProtDcal
filename4 = os.path.join(os.path.dirname(sys.argv[0]), '..', 'outputs/prepared/fixed/' + file + '_Prot.pdb')
mol.renumberResidues()
mol.write(filename4)

command = "mkdir " + os.path.join(os.path.dirname(sys.argv[0]), "ProtDCal_v4.5/Datasets/PMPs/" + file)
os.system(command)

command = "cp " + filename4 + " " + os.path.join(os.path.dirname(sys.argv[0]), "ProtDCal_v4.5/Datasets/PMPs/" + file)
os.system(command)

with open(os.path.join(os.path.dirname(sys.argv[0]), 'ProtDCal_v4.5/Projects/dreamm.proj'), 'r') as file2:
    # read a list of lines into data
    data = file2.readlines()

data[1] = data[1].rstrip() + file + "\n"

command = "mkdir " + os.path.join(os.path.dirname(sys.argv[0]), "ProtDCal_v4.5/Projects/" + file)
os.system(command)
with open(os.path.join(os.path.dirname(sys.argv[0]), 'ProtDCal_v4.5/Projects/' + file + '/' + file + '.proj'), 'w') as file2:
    file2.writelines(data)

os.chdir(os.path.join(os.path.dirname(sys.argv[0]), "ProtDCal_v4.5"))

command = "java -Xms1g -Xmx4g -jar ProtDCal.jar -p Projects/" + file + " -o Outputs/"
os.system(command)

command = "rm -r " + "Projects/" + file
os.system(command)
command = "rm -r " + "Datasets/PMPs/" + file
os.system(command)

with open('Outputs/' + file + '/' + file + '_AA.txt', 'r') as file2:
    data = file2.readlines()

data = data[2:]

with open('Outputs/' + file + '/'+ file + '.csv' , 'w') as file2:
    file2.writelines(data)
    
