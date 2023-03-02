# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 15:57:11 2023

@author: u56903
"""
import pickle
import pandas as pd
from numpy import array, where, nan, isnan, unique, mean, delete, loadtxt
#from io_catalogues import parse_ga_event_query
#from misc_tools import dictlist2array
#from mag_tools import get_au_ml_zone
from obspy import UTCDateTime
from mapping_tools import distance
from data_fmt_tools import return_all_au_station_data


def get_ml_cor(repi, eqdep, legacyML, targetML, ml):
    from calculate_magnitudes import get_ml_corrections
    from numpy import sqrt, nan
    
    logA = 0 # dummy value
    rhyp = sqrt(repi**2 + eqdep**2)
        
    mldict = get_ml_corrections(logA, rhyp, eqdep)
    #print(mldict)
    
    # if legacy is null - set nan
    if legacyML == 'null':
        legacy_corr = nan
        target_corr = nan
        corr_diff = 0.0
        sta_logA = nan
        
    else:
        legacy_corr = mldict[legacyML]
        target_corr = mldict[targetML]
        
        # get difference between ML relationships
        corr_diff = target_corr - legacy_corr
        
        # back-engineer log A value - for testing
        sta_logA = ml - legacy_corr
    
    return {'ml':ml, 'rhyp':rhyp, 'corr_diff':corr_diff, 'legacy_corr':legacy_corr, 
            'target_corr':target_corr, 'sta_logA': sta_logA}
    
###############################################################################
# load data
###############################################################################

# get_station_distance_stadat(sta, eqlo, eqla), 
# first, load pkl
mcdat = pickle.load(open('merged_cat.pkl', 'rb'))

# read station data
sta_dat = return_all_au_station_data()

# load coeffs for ML (pre-1950), MD, or MP correction
r0, r1 = loadtxt('ml_revision_reg.csv', delimiter=',')

###############################################################################
# loop thru events
###############################################################################

for i, mc in enumerate(mcdat):
    #print('\n\n')
    print(mc['DATETIME'])
    
    # reset variables
    add_W_A_correction = False # for changes in W-A magnification: False assumes 2800; True assumes 2080
    add_DC_W_A_correction = False # as above, but where damping = 0.8 for both magnifications
    add_GG91_HV_corr = False # for bugs in implementation of GG91 post MGO: add 0.13 if True
    legacyML = 'null'
    targetML = 'null'
    event_mlcor = nan
    
    # set min distance considered so effects of saturation can be considered - from Allen (2021)
    if mc['PREFML'] >= 4.0 and mc['PREFML'] < 4.5 and mc['DATETIME'] < UTCDateTime(2010,1,1):
        min_considered_dist = 75. # km
    elif mc['PREFML'] >= 4.5 and mc['PREFML'] < 5.0 and mc['DATETIME'] < UTCDateTime(2010,1,1):
        min_considered_dist = 150. # km
    elif mc['PREFML'] >= 5.0 and mc['DATETIME'] < UTCDateTime(2010,1,1):
        min_considered_dist = 250. # km
    else:
        min_considered_dist = 80. # km
        
    
    ###############################################################################    
    # for western & central Australia - set to Gaull & Gregson (1991)
    ###############################################################################
    
    if mc['MLREGION'] == 'WA' or mc['MLREGION'] == 'WCA':
        # get from and to eqns based on date and agency
        if mc['PREFMLSRC'] == 'ADE':
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1968,1,1):
                legacyML = 'R35'
                targetML = 'GG91'
            
            elif mc['DATETIME'] >= UTCDateTime(1968,1,1) \
                and mc['DATETIME'] < UTCDateTime(1998,1,1):
                legacyML = 'GS86'
                targetML = 'GG91'
               
            # assume eqlocl with W-A Gain = 2080 
            elif mc['DATETIME'] >= UTCDateTime(1998,1,1) \
                and mc['DATETIME'] < UTCDateTime(2010,1,1):
                legacyML = 'GS86'
                targetML = 'GG91'
                add_W_A_correction = True
                
            # assume ADE using eqFocus with W-A Gain=2080
            elif mc['DATETIME'] > UTCDateTime(2010,1,1):
                 legacyML = 'BJ84'
                 targetML = 'GG91'
                 add_W_A_correction = True
                 
        
        if mc['PREFMLSRC'] == 'AUST':
            # fix MagCalc bugs where SA ML was used for WA
            if mc['DATETIME'] >= UTCDateTime(2000,1,1) \
               and mc['DATETIME'] < UTCDateTime(2008,7,1):
                   legacyML = 'GS86'
                   targetML = 'GG91'
                   add_GG91_HV_corr = True
                   
            # fix Antelope & SCP bugs
            elif mc['DATETIME'] >= UTCDateTime(2008,7,1):
                legacyML = 'null'
                targetML = 'null'
                add_GG91_HV_corr = True
                add_W_A_correction = True
                
            elif mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'GG91'
                
            # if MGO location, asume MGO ML too
            elif mc['DATETIME'] >= UTCDateTime(1990,1,1) \
                and mc['DATETIME'] < UTCDateTime(2000,1,1):
                legacyML = 'null'
                targetML = 'null'
                        
        # for MGO/BMR
        else:
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'GG91'
            
    ###############################################################################
    # for South Australia/Flinders Ranges region - set to Greenhalgh & Singh (1986)
    ###############################################################################
    
    if mc['MLREGION'] == 'SA':
        if mc['PREFMLSRC'] == 'ADE':
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1968,1,1):
                legacyML = 'R35'
                targetML = 'GS86'
           
            # assume ADE using eqFocus with W-A Gain=2080
            elif mc['DATETIME'] >= UTCDateTime(2010,1,1):
                legacyML = 'BJ84'
                targetML = 'GS86'
                add_W_A_correction = True
            
            # make no change, but assume eqlocl & eqFocus with W-A Gain = 2080 
            elif mc['DATETIME'] >= UTCDateTime(1998,1,1) \
                and mc['DATETIME'] < UTCDateTime(2010,1,1):
                legacyML = 'null'
                targetML = 'null'
                add_W_A_correction = True
            
            # make no change
            elif mc['DATETIME'] >= UTCDateTime(1960,1,1):
                legacyML = 'null'
                targetML = 'null'
        
        elif mc['PREFMLSRC'] == 'AUST':
            # fix MagCalc bugs where SEA ML was used for SA
            if mc['DATETIME'] >= UTCDateTime(2000,1,1) \
               and mc['DATETIME'] < UTCDateTime(2008,7,1): 
                legacyML = 'MLM92'
                targetML = 'GS86'
            
            # add Antelope/SCP W-A correction
            elif mc['DATETIME'] > UTCDateTime(2008,7,1):
                legacyML = 'null'
                targetML = 'null'
                add_DC_W_A_correction = True
                
            elif mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'GS86'
                
            # make no change - assume events recalculated from 2010
            elif  mc['DATETIME'] < UTCDateTime(2010,1,1):
                legacyML = 'null'
                targetML = 'null'
               
        # assume other agencies (including GA/AGSO/BMR) used Richter
        else:
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1990,1,1):
                legacyML = 'R35'
                targetML = 'GS86'
                                        
    ###############################################################################            
    # for Eastern Australia - set to Michael-Lieba & Malafant (1992)
    ###############################################################################
    
    if mc['MLREGION'] == 'SEA' or mc['MLREGION'] == 'EA':
        if mc['PREFMLSRC'] == 'MEL' or mc['PREFMLSRC'] == 'SRC' \
            or mc['PREFMLSRC'] == 'RC' or mc['PREFMLSRC'] == 'GG':
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1998,1,1):
                legacyML = 'R35'
                targetML = 'MLM92'
            
            # assume SRC using eqLocl/eqFocus with W-A Gain=2080
            elif mc['DATETIME'] >= UTCDateTime(1998,1,1) \
                and mc['DATETIME'] < UTCDateTime(2016,1,1):
                legacyML = 'BJ84'
                targetML = 'MLM92'
                add_W_A_correction = True
            
            # assume SRC using eqFocus with W-A Gain=2080
            elif mc['DATETIME'] >= UTCDateTime(2016,1,1):
                 legacyML = 'null'
                 targetML = 'null'
                 add_W_A_correction = True
                
        elif mc['PREFMLSRC'] == 'AUST':
            # fix MagCalc bugs where WA ML was used for SEA
            if mc['DATETIME'] >= UTCDateTime(2000,1,1) \
               and mc['DATETIME'] < UTCDateTime(2008,7,1):
                   legacyML = 'GG91'
                   targetML = 'MLM92'
            
            # add Antelope/SCP W-A correction
            elif mc['DATETIME'] > UTCDateTime(2008,7,1):
                legacyML = 'null'
                targetML = 'null'
                add_W_A_correction = True
                
            elif mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1991,1,1):
                legacyML = 'R35'
                targetML = 'MLM92'
                
        elif mc['PREFMLSRC'] == 'ADE':
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1968,1,1):
                legacyML = 'R35'
                targetML = 'MLM92'
                
            elif mc['DATETIME'] >= UTCDateTime(2010,1,1) \
                and mc['DATETIME'] < UTCDateTime(2018,1,1):
                legacyML = 'BJ84'
                targetML = 'MLM92'
                add_W_A_correction = True
                
            elif mc['DATETIME'] >= UTCDateTime(1960,1,1):
                legacyML = 'GS86'
                targetML = 'MLM92'
            
                
        elif mc['PREFMLSRC'] == 'BRS':
            if mc['DATETIME'] >= UTCDateTime(1960,1,1):
                # assume R35 all the time as per correspondence with Jack Rynn
                legacyML = 'R35'
                targetML = 'MLM92'
            
        # Assumed W-A sensitivity of 2080
        elif mc['PREFMLSRC'] == 'Allen (unpublished)':
            legacyML = 'null'
            targetML = 'null'
            add_DC_W_A_correction = True
        
        # assume other agencies (including GA/AGSO/BMR) used Richter
        else:
            if mc['DATETIME'] >= UTCDateTime(1960,1,1) \
                and mc['DATETIME'] < UTCDateTime(1991,1,1):
                legacyML = 'R35'
                targetML = 'MLM92'
                
    ###############################################################################            
    # get recording stations relative to event
    ###############################################################################
    event_dists = []
    event_stas = []
    for sta in sta_dat:           
           
           if mc['DATETIME'] >= UTCDateTime(sta['startdate']) \
               and mc['DATETIME'] <= UTCDateTime(sta['enddate']):
               
               # calc distance
               dist = distance(mc['LAT'], mc['LON'], sta['stla'], sta['stlo'])[0]
               
               # max dist
               if dist <= 1100.:
                   event_dists.append(dist)
                   event_stas.append(sta['sta'])
    
    # find if any stations between 80-180 km - lower distance changed from 50 km
    event_dists = array(event_dists)
    event_stas = array(event_stas)
    
    sta_mlcor = []
    event_mlcor = nan # set initial value
    ml_corrections = nan
    
    # check depth is not nan
    if isnan(mc['DEP']):
        depth = 0.
    else:
        depth = mc['DEP']
    
    # find duplicate stations - load station sets
    lines = open('station_sets.csv').readlines()
    sta_sets = []
    for line in lines:
        sta_sets.append(set(line.strip().split(',')))
    
    # find duplicate stations - find indexes
    didx = []
    for ss in sta_sets:
        cnt = 0
        didx_temp = []
        for j, sta in enumerate(event_stas):
            if sta in ss:
                cnt += 1
                didx_temp.append(j)
    
        if cnt >= 2:
            didx.append(didx_temp[1])
            
    # remove duplicate stations
    event_dists = delete(event_dists, didx)
    event_stas = delete(event_stas, didx)
    
    #idx = where(event_dists <= 180.)[0]
    #print(event_stas)
    
    ###############################################################################            
    # get corrections
    ###############################################################################
    
    '''
    # if different magnitude types preferred, ignore - still want to keep best ML estimate if not preferred
    if mc['MX_TYPE'] == 'mb' or mc['MX_TYPE'] == 'MS':
        legacyML = 'null'
        targetML = 'null'
        event_mlcor = nan
    '''
    
    # if ML(< 1950), use regression from regress_ml_legacy_target.py
    if mc['DATETIME'] < UTCDateTime(1960,1,1) and isnan(mc['PREFML']) == False:
        event_mlcor = (r0 * mc['PREFML'] + r1) - mc['PREFML']
        legacyML = 'null'
        targetML = 'regression'
        min_repi = nan
    
    # if MD or MP, use regression from regress_ml_legacy_target.py
    elif mc['MX_TYPE'] == 'MP' or mc['MX_TYPE'] == 'MD':
        if isnan(mc['PREFML']) == False:
            event_mlcor = (r0 * mc['PREFML'] + r1) - mc['PREFML']
            legacyML = 'null'
            targetML = 'regression'
            min_repi = nan
                
    else:
        # loop thru stations found for each event    
        if len(event_dists) > 0:
            
            mcd = min_considered_dist
            
            if mcd < 180.:
                idx = where((event_dists >= mcd) & (event_dists <= 180.))[0]
            else:
                idx = []
            
            if len(idx) > 0 and mcd < 180:
                # get corrections for all stations in range
                for repi in event_dists[idx]:
                    ml_corrections = get_ml_cor(repi, depth, legacyML, targetML, mc['PREFML'])
                    min_repi = min(event_dists)
                    
                    sta_mlcor.append(ml_corrections['corr_diff'])
                
                # get final correction
                event_mlcor = mean(array(sta_mlcor))
                #print('Mean ML coor: '+str(event_mlcor))
                #print(sta_mlcor)
                
            # else, get closest distant site     
            else:
                
                idx = where((event_dists > mcd) & (event_dists <= 1100.))[0]
                
                # get minimum distance
                if len(idx) > 0:
                    min_repi = min(event_dists[idx])
                
                    ml_corrections = get_ml_cor(min_repi, depth, legacyML, targetML, mc['PREFML'])
                    
                    # get final correction
                    event_mlcor = ml_corrections['corr_diff']
                    
                else:
                    min_repi = nan
                    event_mlcor = nan
        
        # now add correction for incorrect use of V to H correction for GG91
        #print('\n'+str(event_mlcor))
        if add_GG91_HV_corr == True:
            event_mlcor += 0.13
            
        if add_W_A_correction == True:
            # add approximate mangitude difference after accounting for W-A gain and damping
            '''
            For more details, see: Glanville, H., T. Allen, B. Stepin, and C. Bugden 
            (2020). Seismic monitoring of the NSW CSG areas: monitoring of seismic 
            activity in the CSG production area of Camden and the seismicity of the 
            region, Geoscience Australia Record 2020/20, Canberra, 39 pp, 
            doi: 10.11636/Record.2020.020.
            
            Using ML_2080 dependent relationship as outlined in: The 2023 National 
            Seismic Hazard Assessment for Australia: Earthquake Hypocentre Catalogue
            '''
            # ML2800_corr = m0 * ML2080 + m1
            m0 = -0.00686502963174
            m1 = 0.123141814037
            ML2800_corr = m0 * (mc['PREFML'] + event_mlcor) + m1
            event_mlcor += ML2800_corr
        
        # if W-A damping = 0.8 for both W-A gain settings
        if add_DC_W_A_correction == True:
            event_mlcor += 0.13
        
        # if event_mlcor == 0.0, set to nan so REVML_2023 == nan
        if event_mlcor == 0.0:
            event_mlcor = nan
        
    # add revised ML to mccat
    mcdat[i]['REVML_2023'] = mc['PREFML'] + event_mlcor
    mcdat[i]['legacyML'] = legacyML
    mcdat[i]['targetML'] = targetML
    mcdat[i]['min_repi'] = min_repi
    
    '''
    # for test events
    #print(event_stas)
    print( mc['MLREGION'])
    print(', '.join((legacyML, targetML)))
    #print(event_dists)
    #print(sta_mlcor)
    print(ml_corrections)
    print(event_mlcor)
    print('Rev ML = ' + str(mc['PREFML'] + event_mlcor))
    '''
#####################################################################
# now dump all data to pkl
#####################################################################        

pklfile = open('merged_cat_revised_ml.pkl', 'wb')
pickle.dump(mcdat, pklfile, protocol=-1)
pklfile.close()

#####################################################################
# now convert to pandas and write to csv
#####################################################################        

# convert dicts back to dataframe
mcdf = pd.DataFrame(mcdat)
mcdf.to_csv('merged_cat_revised_ml.csv', sep=',')
