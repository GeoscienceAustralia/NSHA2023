import os, sys
from os.path import join
import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ssc_weights = '/home/547/jdg547/Expert_elicitation/SSC_expert_weights.csv'

ssc_target_file = '/home/547/jdg547/Expert_elicitation/SSC Target Questions for NSHA23 Responses.csv'

exp_weights = {}
f = open(ssc_weights, 'r')
lines = f.readlines()[1:]
f.close()
for line in lines:
    row = line.rstrip('\n').split(',')
    exp_weights[row[0]] = float(row[1])
print(exp_weights)

# For adding catalogue weights to
mw_cat_w = 0
ml_cat_w = 0
combo_cat_w = 0

mw_cat_no_weighting = 0 #Record raw responses without expert calibration
combo_cat_no_weighting = 0

cut_off_w = 0 # Drop models with weights less that this
combo_cat_mw_w = 0 # Weight for mw catalogue if both used
combo_cat_ml_w = 0 # Weight for ml catalogue if both used

decluster_trunc_true = 0
decluster_trunc_false = 0
decluster_trunc_true_no_w = 0 # Without expert weighting
decluster_trunc_false_no_w = 0
decluster_trunc_value = 0

f = open(ssc_target_file)
lines = f.readlines()#[70:]
f.close()
for line in lines:
    print(line)
    row = line.rstrip('\n').split(',')
    exp_ind = row[0]
    print(exp_ind)
    cut_off_w += float(row[3])*exp_weights[exp_ind]
    if row[5] == 'Only the homgenised moment magnitude catalogue developed by Geoscience Australia':
        mw_cat_w += exp_weights[exp_ind]
        mw_cat_no_weighting += 1
    elif row[5] == 'A weighted combination of both catalogues':
        combo_cat_w += exp_weights[exp_ind]
        combo_cat_no_weighting += 1
    else:
        msg = 'Unknown magnitude option'
        raise ValueError(msg, row[5])
    combo_cat_mw_w += float(row[7])*exp_weights[exp_ind]
    combo_cat_ml_w += float(row[10])*exp_weights[exp_ind]
    if row[12] == 'Yes':
        decluster_trunc_true += exp_weights[exp_ind]
        decluster_trunc_true_no_w += 1
        decluster_trunc_value += float(row[14])*exp_weights[exp_ind]
    elif row[12] == 'No':
        decluster_trunc_false += exp_weights[exp_ind]
        decluster_trunc_false_no_w += 1
        decluster_trunc_value += float(row[14])*exp_weights[exp_ind]
    else:
        if row[12] == '':
            print('Skipping blank response for declustering truncation from expert %s' % exp_ind)
        else:
            msg = 'Unknown decluster truncation option'
            raise ValueError(msg, row[12])
    
print('cut_off_w', cut_off_w)
print('mw_cat_w', mw_cat_w)
print('combo_cat_w', combo_cat_w)
print('combo_cat_mw_w', combo_cat_mw_w)
print('combo_cat_ml_w', combo_cat_ml_w)
print('decluster_trunc_true', decluster_trunc_true)
print('decluster_trunc_false', decluster_trunc_false)
print('decluster_trunc_value', decluster_trunc_value)
