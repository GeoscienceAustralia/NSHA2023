# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:36:57 2023

@author: u56903
"""
import pandas as pd
import pickle
from numpy import array, where, isnan, unique, mean
from io_catalogues import parse_ga_event_query
from misc_tools import dictlist2array
from mag_tools import get_au_ml_zone
from obspy import UTCDateTime

###############################################################################
# parse Phil's merged catalogue with revised MLs
###############################################################################

#mcdf = pd.read_csv('Merged_Catalogue_Test.csv')
mcdf = pd.read_csv('Merged_Catalogue.csv')

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
    # add datetime
    mcdat[i]['DATETIME'] = UTCDateTime(mcdat[i]['DATESTR'].replace(' ','T'))
    
    # get ga event index
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
        if ml_zone[0].endswith('Australia') or ml_zone[0].endswith('A'):
            mcdat[i]['MLREGION'] = ml_zone[0].strip('ustralia')
            if mcdat[i]['MLREGION'] == 'SWA':
                mcdat[i]['MLREGION'] = 'WA'
        else:
            mcdat[i]['MLREGION'] = 'Other'
        
###############################################################################
# Append additional NEAC data
###############################################################################

print('append additional neac mags here') # but maybe not, as MLs not great

#####################################################################
# Add missing GCMT data (see ISC cat)
#####################################################################        

# parse mw data
mwfile = 'combined_au_mw.dat'
lines = open(mwfile).readlines()[1:]

mwdat = []
for line in lines:
    dat = line.strip().split(',')
    mdat = {'dt':UTCDateTime(dat[0]), 'lon':float(dat[1]), 'lat':float(dat[2]),
            'mw':float(dat[4]), 'src':dat[5]}
    mwdat.append(mdat)
    
# get unique events and get mean mags
mw_events = unique(dictlist2array(mwdat,'dt'))
for mwe in mw_events:
    mwidx = []
    mwval = []
    mwsrc = []
    for i, mwd in enumerate(mwdat):
        if mwd['dt'] == mwe:
            mwidx.append(i)
            mwval.append(mwd['mw'])
            mwsrc.append(mwd['src'])
    
    # assign mean mag to event (note, will use last src) 
    mean_mw = mean(array(mwval))
    for i in mwidx:
        mwdat[i]['mw'] = mean_mw
        mwdat[i]['src'] = '; '.join(mwsrc)

# now loop through merged cat and add Mws
for i, mc in enumerate(mcdat):
    for mwd in mwdat:
        # add mw data
        if mc['DATETIME'] > mwd['dt']-60. and mc['DATETIME'] < mwd['dt']+60. \
           and mc['MX_ORIGML'] > mwd['mw']-1.0 and mc['MX_ORIGML'] < mwd['mw']+1.0:
            mcdat[i]['PREFMW'] = mwd['mw']
            mcdat[i]['PREFMWSRC'] = mwd['src']
            print(mwd['src'], mwd['mw'])    
        
#####################################################################
# now dump all data to pkl
#####################################################################        

pklfile = open('merged_cat.pkl', 'wb')
pickle.dump(mcdat, pklfile, protocol=-1)
pklfile.close()