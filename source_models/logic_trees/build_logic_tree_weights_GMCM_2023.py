import os, sys
from os.path import join
import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

gmcm_expert_weights = '/home/547/jdg547/Expert_elicitation/GMCM_expert_weights.csv'
cratonic_target_file = '/home/547/jdg547/Expert_elicitation/GMCM Target Questions for NSHA23 - Cratonic sources.csv'
noncratonic_target_file = '/home/547/jdg547/Expert_elicitation/GMCM Target Questions for NSHA23 - Non-cratonic sources.csv'
plate_boundary_target_file = '/home/547/jdg547/Expert_elicitation/GMCM Target Questions for NSHA23 - plate-margin sources and other questions.csv'
gmcm_model_weights_filename = 'GMCM_model_weights_2023.csv'
gmcm_logic_tree_file = 'GMCM_NSHA23_logic_tree.xml'

def largest_remainder(weights, expected_sum=1, precision=0):
    """Use largest remainder method to round weights such that                                                                                                             
    the sum to 1.                                                                                                                                                          
    params weights: The raw weights                                                                                                                                        
    params expected_sum: The total the weights should sum to                                                                                                               
    params precision: Number of decimal places to round to                                                                                                                 
    """
    total_number = expected_sum*np.power(10,precision)
    weights = weights*np.power(10,precision)
    updated_weights = np.floor(weights)
    remainders = weights - updated_weights
    unallocated_places = total_number - np.sum(updated_weights)
    for i in range(int(unallocated_places)):
        max_remainder_index = np.argmax(remainders)
        updated_weights[max_remainder_index] = updated_weights[max_remainder_index] + 1
        remainders[max_remainder_index] = 0
    updated_weights = updated_weights/np.power(10, precision)
    return updated_weights

exp_weights = {}
f = open(gmcm_expert_weights, 'r')
lines = f.readlines()[1:]
f.close()
for line in lines:
    row = line.rstrip('\n').split(',')
    exp_weights[row[0]] = float(row[1])
print(exp_weights)

#Plate Boundary sources
pb_labels = []
pb_weights = []

f = open(plate_boundary_target_file)
lines = f.readlines()
f.close()
for line in lines:
    print(line)
    row = line.rstrip('\n').split(',')
    exp_ind = row[0]
