# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 12:36:19 2023

@author: u56903
"""
import pickle
import pandas as pd
from os import remove, path
from numpy import arange, array, delete, isnan, nan, where, loadtxt
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from mag_tools import get_au_ml_zone
#from misc_tools import dictlist2array
import matplotlib as mpl
mpl.style.use('classic')

###############################################################################
# load data
###############################################################################
# first delete output pkl as doesn't seem to overwrite
if path.isfile('merged_cat_pref_mags.pkl'):
    print('Removing old pkl file...')
    remove('merged_cat_pref_mags.pkl')
else:
    print('pkl file does not exist...')

# second, load revise ml pkl
mcdat = pickle.load(open('merged_cat_revised_ml.pkl', 'rb'))

###############################################################################
# set conversion equations
###############################################################################

# conversion based on W-A 2800 - simulated data
def nsha23_ml2mw(ml):
    c = loadtxt('mw-ml_coeffs_2800.csv', delimiter=',')    
    return c[0] * ml**2 + c[1] * ml + c[2]

# piecewise linear-quadratic based on W-A 2800 - empirical data
def nsha23_piecewise_ml2mw(ml):
    c = loadtxt('mw-ml_lin2quad_coeffs.csv', delimiter=',')  
    
    if ml <= c[4]:
        mw = c[0] * ml + c[1]
    else:
        yhinge = c[0] * c[4] + c[1]
        mw = c[2]*(ml-c[4]) + c[3]*(ml-c[4])**2 + yhinge
    
    return mw

def nsha23_ms2mw(ms):
    c = loadtxt('mw-ms_empirical_coeffs.csv', delimiter=',')    
    return c[0] * ms**2 + c[1]

def nsha23_mb2mw(mb):
    c = loadtxt('mw-mb_empirical_coeffs.csv', delimiter=',')    
    
    if mb < c[3]:
        mw = c[1] + c[0] * mb
    else:
        mhinge = c[1] + c[0] * c[3]
        mw = c[2] * (mb-c[3]) + mhinge        
    return mw

###############################################################################
# first add MS from ISC
###############################################################################

lines =  open('2010-2021_isc_ms.csv').readlines()[1:]
isc_datetime = []
isc_ms = []
for line in lines:
    dat = line.strip().split(',')
    isc_datetime.append(UTCDateTime(dat[0]))
    isc_ms.append(float(dat[1]))
    
# now merge with mcat
print('Adding ISC MS magnitudes')
for i, mc in enumerate(mcdat):
    for j, idt in enumerate(isc_datetime):
        if idt == mc['DATETIME']:
            #if isnan(mc['PREFMS']):
            mcdat[i]['PREFMS'] = isc_ms[j]
            mcdat[i]['PREFMSSRC'] = 'ISC'
            print(idt)

###############################################################################
# now make location adjustments
###############################################################################

lines =  open('event_loaction_adjustment.csv').readlines()[1:]
evdt = []
lons = []
lats = []
deps = []
for line in lines:
    dat = line.strip().split(',')
    evdt.append(UTCDateTime(dat[0]))
    lons.append(float(dat[1]))
    lats.append(float(dat[2]))
    deps.append(float(dat[3]))
    
# now merge with mcat
print('\nAdjusting incorrect earthquake locations')
for i, mc in enumerate(mcdat):
    for j, dt in enumerate(evdt):
        if dt == mc['DATETIME']:
            mcdat[i]['LON'] = lons[j]
            mcdat[i]['LAT'] = lats[j]
            mcdat[i]['DEP'] = deps[j]
            print(dt)

###############################################################################
# remove events not declusterd
###############################################################################

lines =  open('manual_event_removal.csv').readlines()[1:]
evdt = []
for line in lines:
    dat = line.strip().split(',')
    evdt.append(UTCDateTime(dat[0]))
    
# now merge with mcat
print('\nManually removing aftershocks')
mcdat_old = mcdat.copy()
del mcdat
mcdat = []

for i, mc in enumerate(mcdat_old):
    add_event = True
    for j, dt in enumerate(evdt):
        if dt == mc['DATETIME']:
            print(dt)
            add_event = False
            
    if add_event == True:
        mcdat.append(mc)
            
###############################################################################
# set pref ML logic
###############################################################################
ignore_phils_mags = ['ga2013lgzaru', 'ga2021ptpeur', 'ga2021ptxjrd']
cnt = 0
for i, mc in enumerate(mcdat):
    ignore = False
    
    for ipm in ignore_phils_mags:
        if mc['GA_EventID'] == ipm:
            ignore = True
            
    # also ignore if difference too big
    if abs(mc['MLa075_net'] - mc['SeiscompML']) > 1.25:
        ignore = True
    
    #if not available, get ML region
    ml_zone = get_au_ml_zone([mcdat[i]['LON']], [mcdat[i]['LAT']])
    if ml_zone[0].endswith('Australia') or ml_zone[0].endswith('A'):
        mcdat[i]['MLREGION'] = ml_zone[0].strip('ustralia')
        if mcdat[i]['MLREGION'] == 'SWA':
            mcdat[i]['MLREGION'] = 'WA'
    else:
        mcdat[i]['MLREGION'] = 'Other'
        
    
    if mc['PREFMLSRC'] == 'Allen (unpublished)':
        mcdat[i]['PREFML_2023'] = mc['REVML_2023'] # adjusted for W-A sensitivity
        mcdat[i]['PREFMLSRC_2023'] = 'Allen (unpublished-adjust)'
        
    # if available - use recalculated/filtered ML
    elif isnan(mc['MLa05_net']) == False or isnan(mc['MLa075_net']) == False:
        add_GG91_HV_corr = True
        # use corner filter of 0.75 Hz less than ML 4.0
        if mc['MLa075_net'] < 4.0 and mc['MLa075_nob'] >= 3 and ignore == False:
            mcdat[i]['PREFML_2023'] = mc['MLa075_net']
            mcdat[i]['PREFMLSRC_2023'] = 'Cummins (recalc)'
        
        # then use 0.5 Hz filter
        elif mc['MLa05_net'] < 6.0 and mc['MLa05_nob'] >= 3 and ignore == False:
            mcdat[i]['PREFML_2023'] = mc['MLa05_net']
            mcdat[i]['PREFMLSRC_2023'] = 'Cummins (recalc)'
            
        # then use 0.1 Hz filter
        elif mc['MLa01_nob'] >= 3 and ignore == False:
            mcdat[i]['PREFML_2023'] = mc['MLa01_net']
            mcdat[i]['PREFMLSRC_2023'] = 'Cummins (recalc)'
            
        # this was missed in the initial NSHA23 catalogue
        elif isnan(mc['REVML_2023']) == False:
            mcdat[i]['PREFML_2023'] = mc['REVML_2023']
            mcdat[i]['PREFMLSRC_2023'] = 'REV_ML'
            add_GG91_HV_corr = False
            cnt += 1
            
        # take GA mag - some weirdness in phils results
        else:
            mcdat[i]['PREFML_2023'] = mc['PREFML']
            mcdat[i]['PREFMLSRC_2023'] = mc['PREFMLSRC']
            
        # add V-H correction to WA records - this is not double counting - Phil does not do this
        if mcdat[i]['MLREGION'] == 'WA' or mcdat[i]['MLREGION'] == 'WCA':
            if add_GG91_HV_corr == True:
                mcdat[i]['PREFML_2023'] += 0.13
        
                
    # now get next preferred ML
    elif isnan(mc['REVML_2023']) == False:
        mcdat[i]['PREFML_2023'] = mc['REVML_2023']
        mcdat[i]['PREFMLSRC_2023'] = 'REV_ML'
    
    # else use original ML
    elif isnan(mc['PREFML']) == False:
        mcdat[i]['PREFML_2023'] = mc['PREFML']
        mcdat[i]['PREFMLSRC_2023'] = mc['PREFMLSRC']
        
    # else use SCP3 ML
    elif isnan(mc['SeiscompML']) == False:
        mcdat[i]['PREFML_2023'] = mc['SeiscompML']
        mcdat[i]['PREFMLSRC_2023'] = 'SeiscompML'
        
    # else, set to nan
    else:
        mcdat[i]['PREFML_2023'] = nan
        mcdat[i]['PREFMLSRC_2023'] = ''
        
    ###############################################################################
    # now, get pref mags
    ###############################################################################
    # if Mw available
    if isnan(mc['PREFMW']) == False:
        mcdat[i]['PREFMX_2023'] = mc['PREFMW']
        mcdat[i]['PREFMXTYPE_2023'] = 'MW'
        mcdat[i]['PREFMXSRC_2023'] = mc['PREFMWSRC']
        
        mcdat[i]['PREFMW_2023'] = mc['PREFMW']
        mcdat[i]['PREFMWSRC_2023'] = mc['PREFMWSRC']
    
    # if MS > 5.75
    elif mc['PREFMS'] > 5.75:
        mcdat[i]['PREFMX_2023'] = mc['PREFMS']
        mcdat[i]['PREFMXTYPE_2023'] = 'MS'
        mcdat[i]['PREFMXSRC_2023'] = mc['PREFMSSRC']
        
        mcdat[i]['PREFMW_2023'] = nsha23_ms2mw(mc['PREFMS'])
        mcdat[i]['PREFMWSRC_2023'] = 'MS2MW'
        
    # if MS > 5.25 and < 1960
    elif mc['PREFMS'] > 5.25 and mc['DATETIME'] < UTCDateTime(1960,1,1):
        mcdat[i]['PREFMX_2023'] = mc['PREFMS']
        mcdat[i]['PREFMXTYPE_2023'] = 'MS'
        mcdat[i]['PREFMXSRC_2023'] = mc['PREFMSSRC']
        
        mcdat[i]['PREFMW_2023'] = nsha23_ms2mw(mc['PREFMS'])
        mcdat[i]['PREFMWSRC_2023'] = 'MS2MW'
        
    # take any ML
    elif isnan(mc['PREFML_2023']) == False:
        mcdat[i]['PREFMX_2023'] = mc['PREFML_2023']
        mcdat[i]['PREFMXTYPE_2023'] = 'ML'
        mcdat[i]['PREFMXSRC_2023'] = mc['PREFMLSRC_2023']
        
        #mcdat[i]['PREFMW_2023'] = nsha23_piecewise_ml2mw(mc['PREFML_2023']) # using empirical conversion
        mcdat[i]['PREFMW_2023'] = nsha23_ml2mw(mc['PREFML_2023']) # using simulated data conversion
        mcdat[i]['PREFMWSRC_2023'] = 'ML2MW'
        
        #print(nsha23_ml2mw(mc['PREFML_2023']), nsha23_piecewise_ml2mw(mc['PREFML_2023']))
        
    # take larger of mb/MS
    elif isnan(mc['PREFMS']) == False or isnan(mc['PREFmb']) == False:
        
        # choose larger of MS/mb
        if isnan(mc['PREFMS']) == False and isnan(mc['PREFmb']) == False:
            if mc['PREFMS'] >= mc['PREFmb']:
                mcdat[i]['PREFMX_2023'] = mc['PREFMS']
                mcdat[i]['PREFMXTYPE_2023'] = 'MS'
                mcdat[i]['PREFMXSRC_2023'] = mc['PREFMSSRC']
                
                mcdat[i]['PREFMW_2023'] = nsha23_ms2mw(mc['PREFMS'])
                mcdat[i]['PREFMWSRC_2023'] = 'MS2MW'
            else:
                mcdat[i]['PREFMX_2023'] = mc['PREFmb']
                mcdat[i]['PREFMXTYPE_2023'] = 'mb'
                mcdat[i]['PREFMXSRC_2023'] = mc['PREFmbSRC']
                
                mcdat[i]['PREFMW_2023'] = nsha23_mb2mw(mc['PREFmb'])
                mcdat[i]['PREFMWSRC_2023'] = 'mb2MW'
        
        # choose MS
        elif isnan(mc['PREFMS']) == False:
            mcdat[i]['PREFMX_2023'] = mc['PREFMS']
            mcdat[i]['PREFMXTYPE_2023'] = 'MS'
            mcdat[i]['PREFMXSRC_2023'] = mc['PREFMSSRC']
            
            mcdat[i]['PREFMW_2023'] = nsha23_ms2mw(mc['PREFMS'])
            mcdat[i]['PREFMWSRC_2023'] = 'MS2MW'
        
        # choose mb
        else:
            mcdat[i]['PREFMX_2023'] = mc['PREFmb']
            mcdat[i]['PREFMXTYPE_2023'] = 'mb'
            mcdat[i]['PREFMXSRC_2023'] = mc['PREFmbSRC']
            
            mcdat[i]['PREFMW_2023'] = nsha23_mb2mw(mc['PREFmb'])
            mcdat[i]['PREFMWSRC_2023'] = 'mb2MW'
            

#####################################################################
# now dump all data to pkl
#####################################################################        

pklfile = open('merged_cat_pref_mags_post_publication.pkl', 'wb')
pickle.dump(mcdat, pklfile, protocol=-1)
pklfile.close()

#####################################################################
# now convert to pandas and write to csv
#####################################################################        

# convert dicts back to dataframe
mcdf = pd.DataFrame(mcdat)
mcdf.to_csv('merged_cat_pref_mags_post_publication.csv', sep=',')

# number of bad MLs
print(cnt)





