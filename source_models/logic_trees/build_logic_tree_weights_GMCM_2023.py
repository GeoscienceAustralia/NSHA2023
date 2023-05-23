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
#sys.exit()
#Plate Boundary sources
num_pb_gmms = 7 # Number of plate boundary gmm options
pb_labels = ['Allen2012_SS14', 'Allen2022', 'AtkinsonBoore2006',
             'KuehnEtAl2020SSlab', 'KuehnEtAl2020SInter', 'NGAEastGMPE',
             'SomervilleEtAl2009NonCratonic_SS14']
nc_c_labels = ['AbrahamsonEtAl2014', 'Allen2012_SS14', 'AtkinsonBoore2006',
               'BooreEtAl2014', 'CampbellBozorgnia2014', 'ChiouYoungs2008SWISS06',
               'ChiouYoungs2014', 'DrouetBrazil2015', 'ESHM20Craton', 'NGAEastGMPE',
               'RietbrockEdwards2019Mean', 'SomervilleEtAl2009NonCratonic_SS14',
               'SomervilleEtAl2009YilgarnCraton_SS14', 'TromansEtAl2019',
               'ZhaoEtAl2006AscSWISS08']
num_nc_c_gmms = len(nc_c_labels)
pb_weights_equal = np.zeros(num_pb_gmms)
pb_weights = np.zeros(num_pb_gmms)
nc_weights = np.zeros(num_nc_c_gmms)
nc_weights_equal = np.zeros(num_nc_c_gmms)
c_weights = np.zeros(num_nc_c_gmms)
c_weights_equal = np.zeros(num_nc_c_gmms)

# Other parameters
gmm_sigma = np.zeros(4) # Sigma truncation options
gmm_sigma_equal = np.zeros(4) # Sigma truncation options      
gmm_cutoff_w = 0
sigma_truncation_values = [3, 4, 5, 6]
pb_model_weight_indices = [3, 6, 9, 12, 15, 18, 21] # plate boundary
nc_c_model_weight_indices = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, \
                             36, 39, 42, 45] # Non-cratonic and cratonic (same model choices
sigma_weight_indices = [24, 27, 30, 33]
gmm_cutoff_index = 36

f = open(plate_boundary_target_file)
lines = f.readlines()
f.close()

for i, line in enumerate(lines):
    row = line.rstrip('\n').split(',')
    exp_ind = row[0]
    expert_i_weights = []
    expert_i_sigma_weights = []
    for j,ind in enumerate(pb_model_weight_indices):
        try:
            w = float(row[ind])
            expert_i_weights.append(w)
        except:
            msg = 'Weight for model %s for expert %s is not defined, setting to zero' \
                % (pb_labels[j], exp_ind)
            expert_i_weights.append(0.0)
    expert_i_weights = np.array(expert_i_weights)
    if np.allclose(sum(expert_i_weights), 1, rtol=1e-5):
        pass
    else:
        msg = 'Weights for expert %s do not sum to one, renormalised' % exp_ind
        print(msg)
        expert_i_weights = expert_i_weights/sum(expert_i_weights)
    for j,ind in enumerate(pb_model_weight_indices):
        pb_weights_equal[j] = pb_weights_equal[j] + expert_i_weights[j]/len(lines)
        pb_weights[j] = pb_weights[j] + expert_i_weights[j]*exp_weights[exp_ind]

    # Now get sigma truncation and cutoff weight
    for j,ind in enumerate(sigma_weight_indices):
        try:
            w = float(row[ind])
            expert_i_sigma_weights.append(w)
        except:
            msg = 'Sigma truncation weights for expert %s is not defined, setting to zero' \
                % (exp_ind)
            expert_i_sigma_weights.append(0.0)
    expert_i_sigma_weights = np.array(expert_i_sigma_weights)
    if np.allclose(sum(expert_i_sigma_weights), 1, rtol=1e-5):
        pass
    else:
        msg = 'Sigma truncation weights (%s) for expert %s do not sum to one, renormalised' \
            % (expert_i_sigma_weights, exp_ind)
        print(msg)
        expert_i_sigma_weights = expert_i_sigma_weights/sum(expert_i_sigma_weights)
    for j,ind in enumerate(sigma_weight_indices):
        gmm_sigma_equal[j] = gmm_sigma_equal[j] + expert_i_sigma_weights[j]/len(lines)
        gmm_sigma[j] = gmm_sigma[j] + expert_i_sigma_weights[j]*exp_weights[exp_ind]

    # Get cut-off value for lowly weighted models
    gmm_cutoff_w += float(row[gmm_cutoff_index])*exp_weights[exp_ind]

# Now do non-cratonic gmms
f = open(noncratonic_target_file)
lines = f.readlines()
f.close()
for i, line in enumerate(lines):
    row = line.rstrip('\n').split(',')
    exp_ind = row[0]
    expert_i_weights = []
    expert_i_sigma_weights = []
    for j,ind in enumerate(nc_c_model_weight_indices):
        try:
            w = float(row[ind])
            expert_i_weights.append(w)
        except:
            msg = 'Weight for model %s for expert %s is not defined, setting to zero' \
                % (nc_c_labels[j], exp_ind)
            expert_i_weights.append(0.0)
    expert_i_weights = np.array(expert_i_weights)
    if np.allclose(sum(expert_i_weights), 1, rtol=1e-5):
        pass
    else:
        msg = 'Weights for expert %s do not sum to one, renormalised' % exp_ind
        print(msg)
        expert_i_weights = expert_i_weights/sum(expert_i_weights)
    for j,ind in enumerate(nc_c_model_weight_indices):
        nc_weights_equal[j] = nc_weights_equal[j] + expert_i_weights[j]/len(lines)
        nc_weights[j] = nc_weights[j] + expert_i_weights[j]*exp_weights[exp_ind]

# Now do cratonic gmms
f = open(cratonic_target_file)
lines = f.readlines()
f.close()
for i, line in enumerate(lines):
    row = line.rstrip('\n').split(',')
    exp_ind = row[0]
    expert_i_weights = []
    expert_i_sigma_weights = []
    for j,ind in enumerate(nc_c_model_weight_indices):
        try:
            w = float(row[ind])
            expert_i_weights.append(w)
        except:
            msg = 'Weight for model %s for expert %s is not defined, setting to zero' \
                % (nc_c_labels[j], exp_ind)
            expert_i_weights.append(0.0)
    expert_i_weights = np.array(expert_i_weights)
    if np.allclose(sum(expert_i_weights), 1, rtol=1e-5):
        pass
    else:
        msg = 'Weights for expert %s do not sum to one, renormalised' % exp_ind
        print(msg)
        expert_i_weights = expert_i_weights/sum(expert_i_weights)
    for j,ind in enumerate(nc_c_model_weight_indices):
        c_weights_equal[j] = c_weights_equal[j] + expert_i_weights[j]/len(lines)
        c_weights[j] = c_weights[j] + expert_i_weights[j]*exp_weights[exp_ind]

# Write out raw weights before removing lowly weighted models
f_out = open('./plate_boundary_weights_raw.csv', 'w')
for i, label in enumerate(pb_labels):
    print(label, pb_weights[i])
    outline = label +',%.10f\n' % pb_weights[i]
    f_out.write(outline)
f_out.close

f_out = open('./non_cratonic_weights_raw.csv', 'w')
for i, label in enumerate(nc_c_labels):
    print(label, nc_weights[i])
    outline = label +',%.10f\n' % nc_weights[i]
    f_out.write(outline)
f_out.close

f_out = open('./cratonic_weights_raw.csv', 'w')
for i, label in enumerate(nc_c_labels):
    print(label, c_weights[i])
    outline = label +',%.10f\n' % c_weights[i]
    f_out.write(outline)
f_out.close

# Now we want to remove models with weights less than the cutoff weight
# This is done iteratively: The lowest weighted model below the cutoff weight
# is removed and then the weights of the remaining models normalised to sum
# to one. This is repeated until all remaining models have a weight above
# the cut-off value
# Plate-boundary models
min_w = min(w for w in pb_weights if w > 0)
while min_w < gmm_cutoff_w:
    ind = np.argwhere(pb_weights == min_w).flatten()
    print(ind)
    # Set lowest non-zero weight to zero
    for i in ind:
        pb_weights[i] = 0.0
    # Renormalise    
    pb_weights = pb_weights/np.sum(pb_weights)
    min_w = min(w for w in pb_weights if w > 0)    
# Now round weights and ensure they sum to one using largest
# remainder method
print(pb_weights)
pb_weights = largest_remainder(pb_weights, expected_sum=1, precision=2)
print(pb_weights)

# Non-cratonic models
min_w = min(w for w in nc_weights if w > 0)
while min_w < gmm_cutoff_w:
    ind = np.argwhere(nc_weights == min_w).flatten()
    print(ind)
    # Set lowest non-zero weight to zero                                                                                                    
    for i in ind:
        nc_weights[i] = 0.0
    # Renormalise                                                                                                                           
    nc_weights = nc_weights/np.sum(nc_weights)
    min_w = min(w for w in nc_weights if w > 0)
# Now round weights and ensure they sum to one using largest                                                                              
# remainder method                                                                                                                         
print(nc_weights)
nc_weights = largest_remainder(nc_weights, expected_sum=1, precision=2)
print(nc_weights)

# Cratonic models
min_w = min(w for w in c_weights if w > 0)
while min_w < gmm_cutoff_w:
    ind = np.argwhere(c_weights == min_w).flatten()
    print(ind)
    # Set lowest non-zero weight to zero                                                                                                    
    for i in ind:
        c_weights[i] = 0.0
    # Renormalise                                                                                                                           
    c_weights = c_weights/np.sum(c_weights)
    min_w = min(w for w in c_weights if w > 0)
# Now round weights and ensure they sum to one using largest                                                                               
# remainder method                                                                                                                         
print(c_weights)
c_weights = largest_remainder(c_weights, expected_sum=1, precision=2)
print(c_weights)

# Write out results 
f_out = open('./gmm_sigma_trunc_and_cutoff.csv', 'w')
print('\nGMM cutoff weights', gmm_cutoff_w)
outline = 'GMM cutoff weight,%.10f\n' % gmm_cutoff_w
print(outline)
f_out.write(outline)
print('\nSigma truncation weights')
for i, value in enumerate(sigma_truncation_values):
    print(value, gmm_sigma[i])
    outline = 'Sigma truncation level:%i,%.10f\n' % (value, gmm_sigma[i])
    print(outline)
    f_out.write(outline)
f_out.close()

print('\nPlate boundary model weights')
f_out =	open('./plate_boundary_weights.csv', 'w')
for i, label in enumerate(pb_labels):
    print(label, pb_weights[i])
    outline = label +',%.2f\n' % pb_weights[i]
    f_out.write(outline)
f_out.close()
print(sum(pb_weights))
for i, label in enumerate(pb_labels):
    print(label, pb_weights_equal[i])
print(sum(pb_weights_equal))

print('\nNon-cratonic model weights')
f_out = open('./noncratonic_weights.csv', 'w')
for i, label in enumerate(nc_c_labels):
    print(label, nc_weights[i])
    outline = label +',%.2f\n' % nc_weights[i]
    f_out.write(outline)
f_out.close()
print(sum(nc_weights))
print('\nNon-cratonic model weights - equal expert weighting')
for i, label in enumerate(nc_c_labels):
    print(label, nc_weights_equal[i])
print(sum(nc_weights_equal))

print('\nCratonic model weights')
f_out = open('./cratonic_weights.csv', 'w')
for i, label in enumerate(nc_c_labels):
    print(label, c_weights[i])
    outline = label +',%.2f\n' % c_weights[i]
    f_out.write(outline)
f_out.close()
print(sum(c_weights))
print('\nCratonic model weights - equal expert weighting')
for i, label in enumerate(nc_c_labels):
    print(label, c_weights_equal[i])
print(sum(c_weights_equal))

# Now we write-out to nrml file
outline = '<?xml version="1.0" encoding="UTF-8"?>\n'
outline += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
outline += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
outline += '    <logicTree logicTreeID="lt1">\n\n'
outline += '<!-- GMM selection and weights defined through expert elicitation process, using calibration\n' \
    'power of 0.4 and a 75th percentile cut-off to remove lowly weighted branches -->\n\n'
outline += '        <logicTreeBranchingLevel branchingLevelID="bl1">\n\n'
outline += '            <logicTreeBranchSet uncertaintyType="gmpeModel" branchSetID="bs1"\n'
outline += '                    applyToTectonicRegionType="Non_cratonic">\n\n'
for i, gmm_name in enumerate(nc_c_labels):
    weight = nc_weights[i]
    if weight > 0:
        outline += '                <logicTreeBranch branchID="%s">\n' % gmm_name
        outline += '                    <uncertaintyModel>%s</uncertaintyModel>\n' % gmm_name
        outline += '                    <uncertaintyWeight>%.2f</uncertaintyWeight>\n' % weight
        outline += '                </logicTreeBranch>\n\n'
outline += '           </logicTreeBranchSet>\n\n'
outline += '        </logicTreeBranchingLevel>\n\n'
outline += '        <logicTreeBranchingLevel branchingLevelID="bl2">\n\n'
outline += '            <logicTreeBranchSet uncertaintyType="gmpeModel" branchSetID="bs2"\n'
outline += '                    applyToTectonicRegionType="Cratonic">\n\n'
for i, gmm_name in enumerate(nc_c_labels):
    weight = c_weights[i]
    if weight > 0:
        outline += '                <logicTreeBranch branchID="%s">\n' % gmm_name
        outline += '                    <uncertaintyModel>%s</uncertaintyModel>\n' % gmm_name
        outline += '                    <uncertaintyWeight>%.2f</uncertaintyWeight>\n' % weight
        outline += '                </logicTreeBranch>\n\n'
outline += '           </logicTreeBranchSet>\n\n'
outline += '        </logicTreeBranchingLevel>\n\n'
outline += '        <logicTreeBranchingLevel branchingLevelID="bl3">\n\n'
outline += '            <logicTreeBranchSet uncertaintyType="gmpeModel" branchSetID="bs3"\n'
outline += '                    applyToTectonicRegionType="Subduction">\n\n'
for i, gmm_name in enumerate(pb_labels):
    weight = pb_weights[i]
    if weight > 0:
        outline += '                <logicTreeBranch branchID="%s">\n' % gmm_name
        outline += '                    <uncertaintyModel>%s</uncertaintyModel>\n' % gmm_name
        outline += '                    <uncertaintyWeight>%.2f</uncertaintyWeight>\n' % weight
        outline += '                </logicTreeBranch>\n\n'
            
outline += '           </logicTreeBranchSet>\n\n'
outline += '        </logicTreeBranchingLevel>\n\n'
outline += '    </logicTree>\n'
outline += '</nrml>'

# gmm_path = join(target_path, 'ground_motion_results', 'gmm_logic_tree')
xml_filename = 'NSHA23_Aus_GMM_logic_tree.xml'
f_out = open(xml_filename, 'w')
f_out.write(outline)
f_out.close()

