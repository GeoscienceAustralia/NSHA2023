# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:41:07 2017

@author: u56903
"""
from os import path, getcwd, walk
from tools.make_nsha_oq_inputs import make_logic_tree
#from logic_tree import LogicTree
from source_models.utils.utils import largest_remainder
from numpy import array, hstack, unique
from shutil import copyfile

xmllist = []
weighted_smoothing = False # for weighting two adaptive with different Mmin

###############################################################################
# copy regional source models
###############################################################################

relpath = path.join('..', 'zones', '2023_mw')
faultpath = path.join('..', 'faults')
destinationPath = '2023_final'

#print '\n!!! Using Original Magnitudes !!!\n'
#relpath = path.join('..', 'zones', '2012_mx_ge_4.0')
"""
# copy NSHA13
sourceXML = path.join(relpath, 'NSHA13', 'input', 'collapsed', 'NSHA13_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy AUS6
sourceXML = path.join(relpath, 'AUS6', 'input', 'collapsed', 'AUS6_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy DIMAUS
sourceXML = path.join(relpath, 'DIMAUS', 'input', 'collapsed', 'DIMAUS_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])
"""
###############################################################################
# copy background source models
###############################################################################

# copy NSHA13_Background
sourceXML = path.join(relpath, 'NSHA13_Background', 'input', 'collapsed', 'NSHA13_BACKGROUND_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy ARUP
sourceXML = path.join(relpath, 'ARUP', 'input', 'collapsed', 'ARUP_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy ARUP_Background
sourceXML = path.join(relpath, 'ARUP_Background', 'input', 'collapsed', 'ARUP_Background_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy Domains
sourceXML = path.join(relpath, 'Domains_multi_mc', 'input', 'collapsed', 'Domains_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy Leonard08
sourceXML = path.join(relpath, 'Leonard2008', 'input', 'collapsed', 'Leonard2008_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])
'''
# copy SinMcC2016
sourceXML = path.join(relpath, 'SinMcC2016', 'input', 'collapsed', 'SIN_MCC_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])
'''
###############################################################################
# copy seismotectonic source models
###############################################################################

# copy NSHA13
sourceXML = path.join(relpath, 'NSHA13', 'input', 'seismo_collapsed', 'NSHA13_collapsed_NFSM.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy AUS6
sourceXML = path.join(relpath, 'AUS6', 'input', 'seismo_collapsed', 'AUS6_collapsed_NFSM.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy DIMAUS
sourceXML = path.join(relpath, 'DIMAUS', 'input', 'seismo_collapsed', 'DIMAUS_collapsed_NFSM.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy smoothed seismicity source models
###############################################################################

# copy GA adaptive
if weighted_smoothing == True:
    print('not relevant')

# Use four smoothed models
else:
    '''
    sourceXML = path.join('..', 'smoothed_seismicity', 'Cuthbertson2018', 'cuthbertson2018_source_model_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'cuthbertson2018_source_model_banda.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    sourceXML = path.join('..', 'smoothed_seismicity', 'Hall2007_2023', 'Hall2007_2023_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'Hall2007_2023_banda.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #GA adaptive
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness', \
                          'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #GA fixed kernel 
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness', \
                          'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    '''
    print('Skipping Smoothed...')

###############################################################################
# copy smoothed seismicity source models with faults
###############################################################################

# copy GA adaptive
if weighted_smoothing == True:
    print('not relevant')

# Use four smoothed models with faults
else:
    '''
    sourceXML = path.join('..', 'smoothed_seismicity', 'Cuthbertson2018', 'cuthbertson2018_source_model_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'cuthbertson2018_source_model_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    '''
    sourceXML = path.join('..', 'smoothed_seismicity', 'Hall2007_2023', 'Hall2007_2023_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'Hall2007_2023_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    sourceXML = path.join('..', 'smoothed_seismicity', 'Hall2007_2023', 'Hall2007_2023_trunc_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'Hall2007_2023_trunc_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
        
    #GA adaptive
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness', \
                          'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_trunc_decluster', \
                          'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_trunc_decluster_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_trunc_decluster_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #GA fixed kernel 
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness', \
                          'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_trunc_declustered', \
                          'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_trunc_declustered_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_trunc_declustered_banda_nfsm.xml')
    copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    
###############################################################################
# parse weights file
###############################################################################

lines = open('../../shared/2023_seismic_source_model_weights_draft.csv')

class_wgts = []
mod_wgts = []

for line in lines:
    dat = line.strip().split(',')
    if dat[1].startswith('Source_type'):
        class_wgts.append({'class':dat[3], 'wgt':float(dat[4])})
    elif dat[1].startswith('Source_model'):
        mod_wgts.append({'class':dat[2], 'name':dat[3], 'wgt':float(dat[4])})


# set meta dict
# set up metadata dictionary
modelPath = getcwd() # path where source logic tree is to be saved
meta = {'modelPath': modelPath, 'modelFile':'nsha23_source_model_logic_tree.xml', 
        'splitXMLPath': True} # assume source files in job dir 


###############################################################################
# recalibrate source type weights
###############################################################################
'''
#src_wts[0] = 0. # smoothed
#src_wts[1] = 0. # smoothed faults

# rescale source types
src_wts = array(src_wts)/sum(array(src_wts))
'''
###############################################################################
# get weights
###############################################################################
mod_dict = []
src_wts = []
branch_xml = []
# loop throgh source model types
for cls in class_wgts:
    print('\n'+cls['class'])
    src_type_wts = []
    
    # now loop through models within source type
    for mod in mod_wgts:
        if mod['class'] == cls['class']:
            print('    '+mod['name'])
            
            # loop through xml files and find matches
            for xl in xmllist:
                
                if xl.upper().startswith(mod['name'].upper()):
                
                    # do a couple of checks
                    #print mod, xl, mw
                    if mod['name'] == 'NSHA13' and xl.upper().startswith('NSHA13_BACKGROUND'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                    
                    elif mod['name'] == 'NSHA13' and mod['class'] == 'Seismotectonic' and xl.endswith('collapsed.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'DIMAUS' and mod['class'] == 'Seismotectonic' and xl.endswith('collapsed.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'DIMAUS' and mod['class'] == 'Regional' and xl.endswith('NFSM.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'AUS6' and mod['class'] == 'Seismotectonic' and xl.endswith('collapsed.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'AUS6' and mod['class'] == 'Regional' and xl.endswith('NFSM.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'NSHA13' and mod['class'] == 'Regional' and xl.endswith('NFSM.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'GA_adaptive' and mod['class'] == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'GA_adaptive' and mod['class'] == 'Smoothed_faults' and xl.endswith('banda.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'GA_fixed' and mod['class'] == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'GA_fixed' and mod['class'] == 'Smoothed_faults' and xl.endswith('banda.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'Cuthbertson' and mod['class'] == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                    
                    elif mod['name'] == 'Cuthbertson' and mod['class'] == 'Smoothed_faults' and xl.endswith('banda.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                        
                    elif mod['name'] == 'Hall' and mod['class'] == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                    
                    elif mod['name'] == 'Hall' and mod['class'] == 'Smoothed_faults' and xl.endswith('banda.xml'):
                        print('Not adding '+xl+' to '+mod['class']+' set')
                    
                    # else, add file to list
                    else:
                        mod_wt = 1.0
                        
                        # multiply ARUP models by 0.5
                        if mod['name'].startswith('ARUP'):
                           mod_wt = 0.5
                           print('        Modifying ARUP model weight')
                           
                        if mod['name'].startswith('Hall'):
                           mod_wt = 0.5
                           print('        Modifying Risk Frontiers model weight')
                           
                        if mod['name'].startswith('GA_adaptive'):
                           mod_wt = 0.5
                           print('        Modifying GA_adaptive model weight')
                           
                        if mod['name'].startswith('GA_fixed'):
                           mod_wt = 0.5
                           print('        Modifying GA_fixed model weight')
                           
                        
                        # append weights within source type
                        src_wts.append(mod_wt * mod['wgt'] * cls['wgt'])
                        print(mod_wt * mod['wgt'] * cls['wgt'])
                        
                        # append branch file
                        branch_xml.append(xl)
                         
                        #print xl, mod, st, mw
                        
                        # get models actually added
                        #mod_dict.append({'xml':xl, 'model':mod, 'model_wt':mw, 'src_type':st, 'src_wt':sw, 'cml_wt':mw*sw})
                        
# re-normalise source type weights if within type neq 1.0
print(src_wts)
#src_wts = array(src_type_wts) / sum(src_type_wts)
    
                
#print branch_wts

# check weights sum to one!
if not sum(src_wts) == 1.0:
    print('\nWeights do not sum to 1.0!:',sum(src_wts))
    
    # assume testing, so rescale weights
    if sum(src_wts) < 0.95:
        print(src_wts)
        print('\nAre you testing?  Rescaling weights!!!\n')
        #branch_wts = branch_wts / sum(branch_wts)
        
       
# do largest remainder method to make sure numlers 
updated_weights = largest_remainder(src_wts, expected_sum=1.0, precision=4)
#updated_weights = branch_wts

# write source model logic tree
make_logic_tree(branch_xml, updated_weights, meta)

############################################################################################
# plot logic tree
############################################################################################
"""
import matplotlib.pyplot as plt
from misc_tools import dict2array
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import matplotlib.patheffects as PathEffects
#pe = [PathEffects.withStroke(linewidth=2, foreground="w")]
fig = plt.figure(1, figsize=(10, 10))

fig_src_ty = dict2array(mod_dict, 'src_type')
unq_src_ty = unique(fig_src_ty)

def xypts(x, y):
	ypad = 0.05
	xpad = 0.1
	xd = [x-xpad, x+xpad, x+xpad, x-xpad, x-xpad]
	yd = [y+ypad, y+ypad, y-ypad, y-ypad, y+ypad]
	
	return xd, yd, xpad, ypad
	
def isodd(num):
   return num % 2 != 0

# first plot source model
x, y = (0, 0.5)
xd, yd, xpad, ypad = xypts(x, y)
plt.plot(xd, yd, c='k', lw=0.5)
plt.text(x, y, 'Source Model', va='center', ha='center', size=14)

#plt.xlim([-0.2, 1])
#plt.ylim([0, 1])

# draw nwxt level
xc = 0.5
if isodd(len(unq_src_ty)) == False:
	yc = 1. - 1/(2*len(unq_src_ty))
	yinc = 1./(len(unq_src_ty))
else:
	yc = 1. - 1/(len(unq_src_ty)+1)
	yinc = 1./(len(unq_src_ty)-1)
	
y1 = 0.5
for st in unq_src_ty:
	xd, yd, xpad, ypad = xypts(xc, yc)
	plt.plot(xd, yd, c='k', lw=0.5)
	plt.text(xc, yc, st, va='center', ha='center', size=14)	
	
	# draw connecting lines
	x_connect = [xpad, xc/2., xc-xpad]
	y_connect = [0.5, yc, yc]
	
	# draw connection line
	plt.plot(x_connect, y_connect, 'k', lw=0.5)
	
	yc -= yinc
	
plt.show()
	
"""