# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:36:57 2023

@author: u56903
"""
import pandas as pd
import pickle
from numpy import array, where, isnan
from io_catalogues import parse_ga_event_query
from misc_tools import dictlist2array
from mag_tools import get_au_ml_zone

###############################################################################
# parse Phil's merged catalogue with revised MLs
###############################################################################

mcdf = pd.read_csv('Merged_Catalogue.csv')
#mcdf.sort_values('DATESTR',inplace=True)

# get GA IDs
mc_gaid = array(mcdf.GA_EventID)

# convert dataframe to list of dictionaries
mcdat = mcdf.to_dict('records')

# if needed, convert back to dataframe
#mcdf2 = pd.DataFrame(mcdat)

###############################################################################
# parse GA event query for merging
###############################################################################

gadat = parse_ga_event_query('earthquake_query_2017-2022.csv')

gaids = dictlist2array(gadat, 'event_id')

###############################################################################
# loop through merged cat and fill in gaps
###############################################################################

for i, mc in enumerate(mcdat):
    idx = where(mc['GA_EventID'] == gaids)[0]
    
    # if index found, merge
    if len(idx) > 0:
        print(str(i) + ': Merging Event: ' + gaids[idx[0]])
        
        mcdat[i]['DEPENDENCE'] = -1 # for unknown
        mcdat[i]['LOCSRC'] = 'AUST'
        mcdat[i]['LOCSRC'] = 'AUST'
        
        # get pref mw
        if not isnan(gadat[idx[0]]['mag_mww']):
            mcdat[i]['PREFMW'] = gadat[idx[0]]['mag_mww']
            mcdat[i]['PREFMWSRC'] = 'AUST'
        elif not isnan(gadat[idx[0]]['mag_mw']):
            mcdat[i]['PREFMW'] = gadat[idx[0]]['mag_mw']
            mcdat[i]['PREFMWSRC'] = 'AUST'
        elif not isnan(gadat[idx[0]]['mag_mwp']):
            mcdat[i]['PREFMW'] = gadat[idx[0]]['mag_mwp']
            mcdat[i]['PREFMWSRC'] = 'AUST'
            
        # set other mags
        if not isnan(gadat[idx[0]]['mag_ms']):
            mcdat[i]['PREFMS'] = gadat[idx[0]]['mag_ms']
            mcdat[i]['PREFMSSRC'] = 'AUST'
        if not isnan(gadat[idx[0]]['mag_mb']):
            mcdat[i]['PREFmb'] = gadat[idx[0]]['mag_mb']
            mcdat[i]['PREFmbSRC'] = 'AUST'
        if not isnan(gadat[idx[0]]['mag_ml']):
            mcdat[i]['PREFML'] = gadat[idx[0]]['mag_ml']
            mcdat[i]['PREFMLSRC'] = 'AUST'
            
        # get ML zone
        ml_zone = get_au_ml_zone([mc['LON']], [mc['LAT']])
        if ml_zone[0].endswith('Australia'):
            mcdat[i]['MLREGION'] = ml_zone[0].strip('ustralia')
        else:
            mcdat[i]['MLREGION'] = 'Other'
        
#####################################################################
# Add missing GCMT data (see ISC cat)
#####################################################################        

        
#####################################################################
# now dump all data to pkl
#####################################################################        

pklfile = open('merged_cat.pkl', 'wb')
pickle.dump(mcdat, pklfile, protocol=-1)
pklfile.close()