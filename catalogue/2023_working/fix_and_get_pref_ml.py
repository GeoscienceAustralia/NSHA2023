# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 15:57:11 2023

@author: u56903
"""
import pickle
from numpy import array, where, isnan, unique, mean
from io_catalogues import parse_ga_event_query
from misc_tools import dictlist2array
from mag_tools import get_au_ml_zone
from obspy import UTCDateTime
from mapping_tools import distance
from data_fmt_tools import get_station_distance_stadat, return_sta_data, \
                           return_all_au_station_data
'''
def get_station_details(sta):
    # get AU station data for corrections
    
    distkm, azim = get_station_distance_stadat(sta, eqlo, eqla)
'''
###############################################################################
# load data
###############################################################################

# get_station_distance_stadat(sta, eqlo, eqla), 
# first, load pkl
mcdat = pickle.load(open('merged_cat.pkl', 'rb'))

# read station data
sta_dat = return_all_au_station_data()

###############################################################################
# loop thru events
###############################################################################

for i, mc in enumerate(mcdat):
    # reset bool variables
    add_W_A_correction = False # for changes in W-A magnification: False assumes 2800; True assumes 2080
    add_GG91_HV_corr = False # for bugs in implementation of GG91 post MGO: add 0.13 if True
    
    ###############################################################################    
    # for western & central Australia - set to Gaull & Gregson (1991)
    ###############################################################################
    
    if mc['MLREGION'] == 'WCA':
        # get from and to eqns based on date and agency
        if mc['PREFMLSRC'] == 'ADE':
            if mc['DATETIME'] < UTCDateTime(1968,1,1):
                legacyML = 'R35'
                targetML = 'GG91'
            elif mc['DATETIME'] > UTCDateTime(2010,1,1):
                legacyML = 'BJ84'
                targetML = 'GG91'
                add_W_A_correction = True
        
        # for MGO/BMR
        else:
            if mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'GG91'
                
        # fix magcalc bugs
        if mc['DATETIME'] >= UTCDateTime(2000,1,1) \
           and mc['DATETIME'] < UTCDateTime(2008,7,1) and mc['PREFMLSRC'] == 'AUST':
               legacyML = 'GG86'
               targetML = 'GG91'
               add_GG91_HV_corr = True
               
        # fix Antelope & SCP bugs
        elif mc['DATETIME'] >= UTCDateTime(2008,7,1) and mc['PREFMLSRC'] == 'AUST':
            add_GG91_HV_corr = True
        
    ###############################################################################
    # for South Australia/Flinders Ranges region - set to Greenhalgh & Singh (1986)
    ###############################################################################
    
    if mc['MLREGION'] == 'SA':
        if mc['PREFMLSRC'] == 'ADE':
            if mc['DATETIME'] < UTCDateTime(1968,1,1):
                legacyML = 'R35'
                targetML = 'GS86'
            elif mc['DATETIME'] > UTCDateTime(2010,1,1):
                legacyML = 'BJ84'
                targetML = 'GS86'
                add_W_A_correction = True
    
        # assume other agencies (including GA/AGSO/BMR) used Richter
        else:
            if mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'GS86'
                
        # fix magcalc bugs
        if mc['DATETIME'] >= UTCDateTime(2000,1,1) \
           and mc['DATETIME'] < UTCDateTime(2008,7,1) and mc['PREFMLSRC'] == 'AUST':
               legacyML = 'MLM92'
               targetML = 'GS86'
            
    ###############################################################################            
    # for Eastern Australia - set to Michael-Lieba & Malafant (1992)
    ###############################################################################
    
    if mc['MLREGION'] == 'EA':
        if mc['PREFMLSRC'] == 'MEL' or mc['PREFMLSRC'] == 'SRC' \
            or mc['PREFMLSRC'] == 'RC' or mc['PREFMLSRC'] == 'GG':
            if mc['DATETIME'] < UTCDateTime(1998,1,1):
                legacyML = 'R35'
                targetML = 'MLM92'
            elif mc['DATETIME'] >= UTCDateTime(1998,1,1) \
                and mc['DATETIME'] < UTCDateTime(2016,1,1):
                legacyML = 'BJ84'
                targetML = 'MLM92'
                add_W_A_correction = True
        
        elif mc['PREFMLSRC'] == 'ADE':
            if mc['DATETIME'] >= UTCDateTime(1998,1,1) \
                and mc['DATETIME'] < UTCDateTime(2016,1,1):
                legacyML = 'BJ84'
                targetML = 'MLM92'
                add_W_A_correction = True
        
        # assume other agencies (including GA/AGSO/BMR) used Richter
        else:
            if mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'MLM92'
                
        # fix magcalc bugs
        if mc['DATETIME'] >= UTCDateTime(2000,1,1) \
           and mc['DATETIME'] < UTCDateTime(2008,7,1) and mc['PREFMLSRC'] == 'AUST':
               legacyML = 'GG91'
               targetML = 'MLM92'
               
    ###############################################################################            
    # get recording stations relative to event
    ###############################################################################
    trng_dist = []
    for sta in sta_dat:           
           
           if mc['DATETIME'] >= UTCDateTime(sta['startdate']) \
               and mc['DATETIME'] <= UTCDateTime(sta['enddate']):
               
               # calc distance
               dist = distance(mc['LAT'], mc['LON'], sta['stla'], sta['stlo'])[0]
               
               # max dist
               if dist <= 1200.:
                   trng_dist.append(dist)
                   #sta_dist.append(sta['sta'])


