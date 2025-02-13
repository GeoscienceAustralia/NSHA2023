import os, sys
from os.path import join
import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ssc_expert_weights = '/home/547/jdg547/Expert_elicitation/SSC_expert_weights.csv'
ssc_target_file = '/home/547/jdg547/Expert_elicitation/SSC_Target_Questions_revised.csv' # Weights after second round of EE
#ssc_target_file = '/home/547/jdg547/Expert_elicitation/SSC Target Questions for NSHA23 Responses.csv'
ssc_model_weights_filename = 'SSC_model_weights_2023.csv'

exp_weights = {}
f = open(ssc_expert_weights, 'r')
lines = f.readlines()[1:]
f.close()
for line in lines:
    print(line)
    row = line.rstrip('\n').split(',')
    exp_weights[row[0]] = float(row[1])
print(exp_weights)

labels = []
values = []

# For adding catalogue weights to
mw_cat_w = 0
ml_cat_w = 0
combo_cat_w = 0

mw_cat_no_weighting = 0 #Record raw responses without expert calibration
combo_cat_no_weighting = 0

cut_off_w = 0 # Drop models with weights less that this
cut_off_w_equal_experts = 0 # See what it is if we weight everyone equally
combo_cat_mw_w = 0 # Weight for mw catalogue if both used
combo_cat_ml_w = 0 # Weight for ml catalogue if both used

decluster_trunc_true = 0
decluster_trunc_false = 0
decluster_trunc_true_no_w = 0 # Without expert weighting
decluster_trunc_false_no_w = 0
decluster_trunc_value = 0

background_b_individ_sz = 0
background_b_aggregate_neodom = 0
background_b_centroid_neodom = 0
reg_b_individ_sz = 0
reg_b_aggregate_neodom = 0
reg_b_centroid_neodom = 0
ss_b_neodom = 0
ss_b_superdom = 0
mmin_4p0 = 0
mmin_4p5 = 0
mmin_5p0 = 0
bg_w = 0
reg_w = 0
seismo_w = 0
ss_w = 0
ss_fsm_w = 0

adaptive_ss_w = 0
fixed_ss_w = 0

f = open(ssc_target_file)
lines = f.readlines()
f.close()
for line in lines:
    print(line)
    row = line.rstrip('\n').split(',')
    exp_ind = row[0]
    print(exp_ind)
    cut_off_w += float(row[3])*exp_weights[exp_ind]
    cut_off_w_equal_experts += float(row[3])
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

    # Check weights sum to one - background model b-values
    bg_b_sum = float(row[17])+ float(row[20]) + float(row[23])
    if np.allclose(bg_b_sum, 1.0, rtol=1e-05):
        background_b_individ_sz += float(row[17])*exp_weights[exp_ind]
        background_b_aggregate_neodom += float(row[20])*exp_weights[exp_ind]
        background_b_centroid_neodom += float(row[23])*exp_weights[exp_ind]
    else:
        msg = 'Background zone b-value weights (%s, %s, and %s) do not sum to one for expert %s' \
            % (row[17], row[20], row[23], exp_ind)
        print(msg)
        print('Renormalising to sum to one')
        #        raise ValueError(msg)
        background_b_individ_sz += float(row[17])/bg_b_sum*exp_weights[exp_ind]
        background_b_aggregate_neodom += float(row[20])/bg_b_sum*exp_weights[exp_ind]
        background_b_centroid_neodom += float(row[23])/bg_b_sum*exp_weights[exp_ind]
        bg_b_sum_normalised = (float(row[17])/bg_b_sum + float(row[20])/bg_b_sum +\
                               float(row[23])/bg_b_sum)
        print(bg_b_sum_normalised)

    # Check weights sum to one - regional model b-values
    reg_b_sum = float(row[26])+ float(row[29]) + float(row[32])
    if np.allclose(reg_b_sum, 1.0, rtol=1e-05):
        reg_b_individ_sz += float(row[26])*exp_weights[exp_ind]
        reg_b_aggregate_neodom += float(row[29])*exp_weights[exp_ind]
        reg_b_centroid_neodom += float(row[32])*exp_weights[exp_ind]
    else:
        msg = 'Regional zone b-value weights (%s, %s, and %s) do not sum to one for expert %s' \
            % (row[26], row[29], row[32], exp_ind)
        print(msg)
        print('Renormalising to sum to one')
        #        raise ValueError(msg)
        reg_b_individ_sz += float(row[26])/reg_b_sum*exp_weights[exp_ind]
        reg_b_aggregate_neodom += float(row[29])/reg_b_sum*exp_weights[exp_ind]
        reg_b_centroid_neodom += float(row[32])/reg_b_sum*exp_weights[exp_ind]
        reg_b_sum_normalised = (float(row[26])/reg_b_sum + float(row[29])/reg_b_sum +\
                               float(row[32])/reg_b_sum)
        print(reg_b_sum_normalised)

    # Check weights sum to one - smoothed seismicity model b-values
    ss_b_sum = float(row[35])+ float(row[38])
    if np.allclose(ss_b_sum, 1.0, rtol=1e-05):
        ss_b_neodom += float(row[35])*exp_weights[exp_ind]
        ss_b_superdom += float(row[38])*exp_weights[exp_ind]
    else:
        msg = 'Smoothed seismicity b-value weights (%s, and %s) do not sum to one for expert %s' \
            % (row[35], row[38], exp_ind)
        print(msg)
        print('Renormalising to sum to one')
        ss_b_neodom += float(row[35])/ss_b_sum*exp_weights[exp_ind]
        ss_b_superdom += float(row[38])/ss_b_sum*exp_weights[exp_ind]
        ss_b_sum_normalised = float(row[35])/ss_b_sum + float(row[38])/ss_b_sum
        print(ss_b_sum_normalised)

    # Check weights sum to one - mmin choice
    mmin_sum = float(row[41])+ float(row[44]) + float(row[47])
    if np.allclose(mmin_sum, 1.0, rtol=1e-05):
        mmin_4p0 += float(row[41])*exp_weights[exp_ind]
        mmin_4p5 += float(row[44])*exp_weights[exp_ind]
        mmin_5p0 += float(row[47])*exp_weights[exp_ind]
    else:
        msg = 'MMin weights (%s, %s, and %s) do not sum to one for expert %s' \
            % (row[41], row[44], row[47], exp_ind)
        print(msg)
        print('Renormalising to sum to one')
        mmin_4p0 += float(row[41])/mmin_sum*exp_weights[exp_ind]
        mmin_4p5 += float(row[44])/mmin_sum*exp_weights[exp_ind]
        mmin_5p0 += float(row[47])/mmin_sum*exp_weights[exp_ind]
        mmin_sum_normalised = float(row[41])/mmin_sum + float(row[44])/mmin_sum +\
            float(row[47])
        print(mmin_sum_normalised)

    # Check weights sum to one - source model interclass weights
    sm_sum = float(row[50])+ float(row[53]) + float(row[56]) + float(row[59]) +float(row[62])
    if np.allclose(sm_sum, 1.0, rtol=1e-05):
        bg_w += float(row[50])*exp_weights[exp_ind]
        reg_w += float(row[53])*exp_weights[exp_ind]
        seismo_w += float(row[56])*exp_weights[exp_ind]
        ss_w += float(row[59])*exp_weights[exp_ind]
        ss_fsm_w += float(row[62])*exp_weights[exp_ind]
    else:
        msg = 'Source model class weights (%s, %s, %s, %s, and %s) do not sum to one for expert %s' \
            % (row[50], row[53], row[56], row[59],row[62], exp_ind)
        print(msg)
        print('Renormalising to sum to one')
        bg_w += float(row[50])/sm_sum*exp_weights[exp_ind]
        reg_w += float(row[53])/sm_sum*exp_weights[exp_ind]
        seismo_w += float(row[56])/sm_sum*exp_weights[exp_ind]
        ss_w += float(row[59])/sm_sum*exp_weights[exp_ind]
        ss_fsm_w += float(row[62])/sm_sum*exp_weights[exp_ind]
        sm_sum_normalised = float(row[50])/sm_sum + float(row[53])/sm_sum + \
            float(row[56])/sm_sum + float(row[59])/sm_sum + float(row[62])/sm_sum
        print(sm_sum_normalised)
        
    # Check weights all sum to one - smoothed seismicity model types
    ss_type_sum = float(row[65])+ float(row[68])
    if np.allclose(ss_type_sum, 1.0, rtol=1e-05):
        fixed_ss_w += float(row[65])*exp_weights[exp_ind]
        adaptive_ss_w += float(row[68])*exp_weights[exp_ind]
    else:
        msg = 'Smoothed seismicity model type weights (%s, and %s) do not sum to one for expert %s' \
            % (row[65],row[68], exp_ind)
        print(msg)
        print('Renormalising to sum to one')
        fixed_ss_w += float(row[65])/ss_type_sum*exp_weights[exp_ind]
        adaptive_ss_w += float(row[68])/ss_type_sum*exp_weights[exp_ind]
        ss_type_sum_normalised = float(row[65])/ss_type_sum + float(row[68])/ss_type_sum
        print(ss_type_sum_normalised)
        
# Write out raw source weights
sm_weights = np.array([bg_w, reg_w, seismo_w, ss_w, ss_fsm_w])
labels = ['Background class weight', 'Regional class weight', \
    'Seismotectonic class weight', 'Smoothed seismicity class weight', \
    'Smoothed seismicity FSM class weight']
f_out = open('source_class_weights_raw.csv', 'w')
for i, label in enumerate(labels):
    outline = label + ',%.10f\n' % sm_weights[i]
    f_out.write(outline)
f_out.close()

# Remove source weight class with total weight < cut-off value                                                                       
print(bg_w)
min_sm_weight_ind = np.argmin(sm_weights)
print(sm_weights)
print(min_sm_weight_ind)
print(sm_weights[min_sm_weight_ind])
if sm_weights[min_sm_weight_ind] < cut_off_w:
    sm_weights[min_sm_weight_ind] = 0
    tmp_sum = np.sum(sm_weights)
    # Renormalise other weights                                                                                                                                                                  
    sm_weights = sm_weights/tmp_sum
    bg_w = sm_weights[0]
    reg_w = sm_weights[1]
    seismo_w = sm_weights[2]
    ss_w = sm_weights[3]
    ss_fsm_w = sm_weights[4]



print('cut_off_w', cut_off_w)
print('cut_off_w_equal_weighting', cut_off_w_equal_experts/14)
print('mw_cat_w', mw_cat_w)
print('combo_cat_w', combo_cat_w)
print('mw_cat_equal_weights', mw_cat_no_weighting/14)
print('combo_cat_equal_weights', combo_cat_no_weighting/14)
print('combo_cat_mw_w', combo_cat_mw_w)
print('combo_cat_ml_w', combo_cat_ml_w)
print('decluster_trunc_true', decluster_trunc_true)
print('decluster_trunc_false', decluster_trunc_false)
print('decluster_trunc_value', decluster_trunc_value)
print('background_b_individ_sz', background_b_individ_sz)
print('background_b_aggregate_neodom', background_b_aggregate_neodom)
print('background_b_centroid_neodom', background_b_centroid_neodom)
print('reg_b_individ_sz', reg_b_individ_sz)
print('reg_b_aggregate_neodom', reg_b_aggregate_neodom)
print('reg_b_centroid_neodom', reg_b_centroid_neodom)
print('ss_b_neodom', ss_b_neodom)
print('ss_b_superdom', ss_b_superdom)
print('mmin_4p0', mmin_4p0)
print('mmin_4p5', mmin_4p5)
print('mmin_5p0', mmin_5p0)
print('bg_w', bg_w)
print('reg_w', reg_w)
print('seismo_w', seismo_w)
print('ss_w', ss_w)
print('ss_fsm_w ', ss_fsm_w)
print('adaptive_ss_w', adaptive_ss_w)
print('fixed_ss_w', fixed_ss_w)

# Write out to file
labels = ['Cut off weight', 'Mw catalogue only', 'Combined catalogues', 'Combined catalogues Mw weight', \
          'Combined catalogues Ml weight', 'Decluter truncation = True', 'Decluster truncation = False', \
          'Decluster truncation value', 'Bg b individual source zone', 'Bg b aggregated domains', \
          'Bg b centroid domains', 'Regional b individual source zone', 'Regional b aggregated domains', \
          'Regional b centroid domains', 'Smoothed seismicity b neo domains', 'Smoothed seismicity b neo superdomains', \
          'Mmin=4.0', 'Mmin=4.5', 'Mmin=5.0', 'Background class weight', 'Regional class weight', \
          'Seismotectonic class weight', 'Smoothed seismicity class weight', 'Smoothed seismicity FSM class weight', \
          'Adaptive SS weight', 'Fixed SS weight']
values = [cut_off_w, mw_cat_w, combo_cat_w, combo_cat_mw_w, combo_cat_ml_w, decluster_trunc_true, decluster_trunc_false, decluster_trunc_value,\
          background_b_individ_sz, background_b_aggregate_neodom, background_b_centroid_neodom, reg_b_individ_sz, \
          reg_b_aggregate_neodom, reg_b_centroid_neodom, ss_b_neodom, ss_b_superdom, mmin_4p0, mmin_4p5, mmin_5p0, \
          bg_w, reg_w, seismo_w, ss_w, ss_fsm_w, adaptive_ss_w, fixed_ss_w]

print(labels, len(labels))
print(values, len(values))
f_out = open(ssc_model_weights_filename, 'w')
for i, label in enumerate(labels):
    outline = label + ',%.2f\n' % values[i]
    f_out.write(outline)
f_out.close()

# Now do some plots
# Catalogue type
if not os.path.exists('plots'):
    os.mkdir('plots')
width = 0.2
x_vals = np.array([0.2, 0.8])
fig, ax = plt.subplots()
n_experts = 14
ax.bar(x_vals, [mw_cat_no_weighting/n_experts, combo_cat_no_weighting/n_experts],
       width, fill=False, hatch='///', label='Equal expert weighting',
       color = '0.7')
ax.bar(x_vals+width, [mw_cat_w, combo_cat_w] , width, color='0.5', label='Calibrated expert weighting')
ax.set_xlim([0,1.2])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')#, linewidth=0.5)
ax.text(0.95, 0.125, 'Cut-off weight')
ax.plot([0,2],[cut_off_w_equal_experts, cut_off_w_equal_experts], linestyle = 'dotted', c = 'b')
plt.xticks(x_vals + width / 2,
           ('Mw catalogue only', 'Weighted combination\n of Mw and Ml\n catalogues'),
           fontsize=14)
# Finding the best position for legends and putting it
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/catalogue_selection.png')

# Catalogue weights
plt.clf()
x_vals = np.array([0.5])
fig, ax = plt.subplots()
ax.bar(x_vals, [combo_cat_ml_w],
       width, fill=False, hatch='///', label='Ml catalogue_weight',
       color = '0.7')
ax.bar(x_vals+width, [combo_cat_mw_w], width, color='0.5', label='Mw catalogue weight')
ax.set_xlim([0,1.2])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')#, linewidth=0.5)                                              
ax.text(0.95, 0.125, 'Cut-off weight')
#plt.xticks(x_vals + width / 2,
#           ('Catalogue choice'),
#           fontsize=14)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/catalogue_weights.png')

# b-value calculation methods - background zones
plt.clf()
x_vals = np.array([0.2, 0.6, 1.0])
fig, ax = plt.subplots()
#ax.bar(x_vals, [,
#       width, fill=False, hatch='///', label='Ml catalogue_weight',
#       color = '0.7')
ax.bar(x_vals,
       [background_b_individ_sz, background_b_aggregate_neodom, background_b_centroid_neodom],
       width, color='0.5', label='b-value method,\nbackground zones')
ax.set_xlim([0,1.2])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')#, linewidth=0.5)                                              
ax.text(0.95, 0.125, 'Cut-off weight')
plt.xticks(x_vals,
           ('Individual\nsource zone', 'Aggregated by\nneotectonic domains',
            'Centroid relative to\nneotectonic domains'),
           fontsize=11)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/bvalue_background_weights.png')

# b-value calculation methods - regional zones                                                                                    
plt.clf()
x_vals = np.array([0.2, 0.6, 1.0])
fig, ax = plt.subplots()
ax.bar(x_vals,
       [reg_b_individ_sz, reg_b_aggregate_neodom, reg_b_centroid_neodom],
       width, color='0.5', label='b-value method,\nregional zones')
ax.set_xlim([0,1.2])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')#, linewidth=0.5)                                              
ax.text(0.95, 0.125, 'Cut-off weight')
plt.xticks(x_vals,
           ('Individual\nsource zone', 'Aggregated by\nneotectonic domains',
            'Centroid relative to\nneotectonic domains'),
           fontsize=11)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/bvalue_regional_weights.png')

# b-value calculation methods - smoothed seismicity models
plt.clf()
x_vals = np.array([0.2, 0.6])
fig, ax = plt.subplots()
ax.bar(x_vals,
       [ss_b_neodom, ss_b_superdom],
       width, color='0.5', label='b-value method,\nsmoothed seismicity')
ax.set_xlim([0,1.0])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')#, linewidth=0.5)                                                                           
ax.text(0.75, 0.125, 'Cut-off weight')
plt.xticks(x_vals,
           ('Neotectonic domains', 'Neotectonic super domains'),
           fontsize=11)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/bvalue_smoothed seismicity_weights.png')

# Smoothed seismicity models method
plt.clf()
x_vals = np.array([0.2, 0.6])
fig, ax = plt.subplots()
ax.bar(x_vals,
       [adaptive_ss_w, fixed_ss_w],
       width, color='0.5', label='Smoothed seismiciy\nmethod')
ax.set_xlim([0,1.0])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')#, linewidth=0.5)                                                                         
ax.text(0.75, 0.125, 'Cut-off weight')
plt.xticks(x_vals,
           ('Adaptive', 'Fixed'),
           fontsize=11)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/smoothed seismicity_method_weights.png')

# Hazard calculation Mmin                                                                                                                           
plt.clf()
x_vals = np.array([0.2, 0.6, 1.0])
fig, ax = plt.subplots()
ax.bar(x_vals,
       [mmin_4p0, mmin_4p5, mmin_5p0],
       width, color='0.5', label='Minimum magnitude\n(for hazard calculations)')
ax.set_xlim([0,1.2])
ax.set_ylim([0, 1.0])
ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')
ax.text(0.95, 0.125, 'Cut-off weight')
plt.xticks(x_vals,
           ('4.0', '4.5', '5.0'),
           fontsize=11)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/minimum_magnitude_weights.png')

#Plot expert weights
plt.clf()
ews = []
keys = []
print(exp_weights)
for key, value in exp_weights.items():
    keys.append(key)
    ews.append(float(value))
exp_weights = np.array(ews)
print(exp_weights)
x_vals = np.arange(len(exp_weights))
fig, ax = plt.subplots()
ax.bar(x_vals,
       exp_weights, 
       width, color='0.5', label='Expert weights')
#ax.set_xlim([0,1.2])
ax.set_ylim([0, 1.0])
#ax.plot([0,2],[cut_off_w, cut_off_w], linestyle = 'dotted', c = 'r')
#ax.text(0.95, 0.125, 'Cut-off weight')
#plt.xticks(x_vals,
#           ('4.0', '4.5', '5.0'),
#           fontsize=11)
plt.legend(loc='best')
plt.ylabel('Weight', fontsize=14)
plt.tight_layout()
plt.savefig('plots/SSC_expert_weights.png')
exp_weights
