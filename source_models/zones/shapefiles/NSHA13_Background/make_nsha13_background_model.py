import shapefile
from shapely.geometry import Polygon
from numpy import ones_like, array
try:
    from tools.nsha_tools import get_field_data
    from tools.source_shapefile_builder import get_preferred_catalogue, \
                                               get_completeness_model, get_aus_shmax_vectors, \
                                               get_rate_adjust_factor, build_source_shape, \
                                               get_ul_seismo_depths, get_neotectonic_domain_params, \
                                               aggregate_intraslab_sources
except:
    print('Add PYTHONPATH to NSHA23 root directory')

###############################################################################

''' START MAIN CODE HERE '''
   
###############################################################################
# parse NSHA13 shp and prep data
###############################################################################

domshp = 'NSHA13_Background_NSHA18_Merged.shp'

print('Reading source shapefile...')
sf = shapefile.Reader(domshp)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get field data
src_codes = get_field_data(sf, 'code', 'str')
src_names1 = get_field_data(sf, 'Name', 'str')
src_names2 = get_field_data(sf, 'SRC_NAME', 'str')
domains = get_field_data(sf, 'DOMAIN', 'float')
mmax1 = get_field_data(sf, 'max_mag', 'float')
mmax2 = get_field_data(sf, 'MMAX_BEST', 'float')
trt = get_field_data(sf, 'trt', 'str')
usd = get_field_data(sf, 'usd', 'float')
lsd = get_field_data(sf, 'lsd', 'float')
hd = get_field_data(sf, 'hd1', 'float')
stk = get_field_data(sf, 'strike1', 'float')
dip = get_field_data(sf, 'dip1', 'float')
rke = get_field_data(sf, 'rake1', 'float')

# merge source names
src_names = []
mmax = []
for sn1, sn2 in zip(src_names1, src_names2):
    if sn1.strip() != '':
        src_names.append(sn1)
    elif sn2.strip() != '':
        src_names.append(sn2)
    else:
        src_names.append('null')
        
# merge Mmax
for mm1, mm2 in zip(mmax1, mmax2):
    if mm1 == 0.0:
        mmax.append(mm2)
    elif mm2 == 0.0:
        mmax.append(mm1)
    else:
        mmax.append(nan)
        
# set domain for unset domains
trt_new = []
for i in range(0,len(trt)):
    if trt[i] == 'Oceanic':
        domains[i] = 8
    elif trt[i] == 'Active':
        domains[i] = 9
    elif trt[i] == 'Extended' and domains[i] == 0.:
        domains[i] = 7
    elif trt[i] == 'Interface':
        domains[i] = 10
    elif trt[i] == 'Intraslab':
        domains[i] = 11
    
    if trt[i] == 'NCratonic':
        trt_new.append('Non_cratonic')
    else:
        trt_new.append(trt[i])

###############################################################################
# load neotectonic domains parameters
###############################################################################

'''
# set domestic domain numbers based on neotectonic domains
neo_domains, neo_min_rmag, neo_mmax, neo_trt, neo_bval_fix, neo_bval_sig_fix \
    = get_neotectonic_domain_params(sf, trt_new)
'''    
# set domestic domain numbers based on neotectonic domains
refshpfile = '..//reference_shp.txt'
#refshpfile = '..//reference_shp_mx.txt' # for testing only!
neo_domains, neo_min_rmag, neo_mmax, neo_trt, neo_bval_fix, neo_bval_sig_fix \
    = get_neotectonic_domain_params(sf, trt_new, refshpfile)

    
for i in range(0, len(domains)):
    if neo_domains[i] > 0 and neo_domains[i] < 8:
        domains[i] = neo_domains[i]
        mmax[i] = neo_mmax[i]

zone_class = list(domains)[:]
'''
# reset Gawler Craton to Flinders due to b-value similarities
zone_class[0] = 2.

# reset West Coast to extended
zone_class[7] = 6.
domains[7] = 6.

# reset West Coast Passive Margin to extended
zone_class[10] = 7.
domains[10] = 7.

# reset Ottway/Gippsland to extended
zone_class[6] = 5.
domains[6] = 5.


'''
# reset Southern Oceanic buffer
zone_class[51] = 8.
domains[51] = 8.
zone_class[39] = 8.
domains[39] = 8.

# reset EBGZ margin extended
zone_class[54] = 7.
domains[54] = 7.

# reset Northwest buffer
zone_class[52] = 7.
domains[52] = 7.

# reset NWO to Oceanic
zone_class[26] = 8.
domains[26] = 8.

# reset Tasmania to Non-cratonic
zone_class[41] = 4.
domains[41] = 4.

###############################################################################
#  set intraslab aggregation class
###############################################################################
print('\n!!! REMEMBER TO RESET SRM_200_300 SOURCE CODE !!!!\n')
for i, src_code in enumerate(src_codes):
   zone_class[i] =  aggregate_intraslab_sources(src_code, zone_class[i])
      
###############################################################################
#  set pref strike/dip/rake
###############################################################################

pref_stk = []
pref_dip = []
pref_rke = []
for i in range(0,len(stk)):
    if stk[i] == 0.0 and domains[i] <= 8:
        pref_stk.append(-999)
        pref_dip.append(-999)
        pref_rke.append(-999)
    else:
        pref_stk.append(stk[i])
        pref_dip.append(dip[i])
        pref_rke.append(rke[i])
        
# set depth parameters
dep_b = []
dep_u = []
dep_l = []
for i in range(0,len(lsd)):
    if domains[i] < 8:
        lsd[i] = 20.
        if trt[i] == 'Cratonic':
            dep_b.append(5.0)
            dep_u.append(2.5)
            dep_l.append(10.)
        else:
            dep_b.append(10.)
            dep_u.append(5.0)
            dep_l.append(15.)
    
    # don't care about depth logic tree outside Australia
    else:
        dep_b.append(hd[i])
        dep_u.append(-999)
        dep_l.append(-999)

# fix preferred upper/lower seismo depths from Domains
usd, lsd = get_ul_seismo_depths(src_codes, usd, lsd)

###############################################################################
# get preferred catalogues 
###############################################################################

# get preferred catalogues for each zone
prefCat = get_preferred_catalogue(domshp)

# fix catalogue for source zones
prefCat[38] = 'NSHA23CAT_V0.1_hmtk_declustered.csv'
prefCat[39] = 'NSHA23CAT_V0.1_hmtk_declustered.csv'
prefCat[54] = 'NSHA23CAT_V0.1_hmtk_declustered.csv'

###############################################################################
# load 2018 completeness models
###############################################################################
single_mc = 0
ycomp, mcomp, min_rmag_ignore = get_completeness_model(src_codes, shapes, domains, single_mc)

min_rmag = neo_min_rmag

# use manual modification
for i in range(0,len(trt_new)):
    if trt_new[i] == 'Active':
        min_rmag[i] = 5.75
    elif trt_new[i] == 'Intraslab':
        min_rmag[i] = 5.75

min_rmag[26] = 3.5 # NWO
min_rmag[3] = 6.1 # TAFS
min_rmag[11] = 6.0 # NBOT
min_rmag[12] = 6.1 # NBT
min_rmag[54] = 3.3 # EBGZ


'''
min_rmag[45] = 3.5 # NECS
#min_rmag[50] = 3.2 # CARP
min_rmag[19] = 3.5 # SEOB
min_rmag[18] = 3.5 # SWOB
min_rmag[14] = 3.2 # ZN6b
min_rmag[44] = 3.5 # TP
#min_rmag[4]  = 3.1 # ZN1b - Gawler, cases b to skyrocket!
'''
'''
min_rmag[52] = 3.3 # KMBY
min_rmag[55] = 3.3 # NACR
min_rmag[54] = 3.3 # NAOR
min_rmag[49] = 3.3 # PLBR
min_rmag[53] = 3.5 # WAPM
min_rmag[48] = 3.2 # YLGN
min_rmag[66] = 3.5 # NWB1
'''

# SEOB - multi-corner
ycomp[51] = '1980;1964;1900'
mcomp[51] = '3.5;5.0;6.0'
min_rmag[51] = 3.5

###############################################################################
# load Rajabi SHMax vectors 
###############################################################################

shmax_pref, shmax_sig = get_aus_shmax_vectors(src_codes, shapes)

###############################################################################
# get rate adjustment factors 
###############################################################################

origshp = 'NSHA13_BACKGROUND_NSHA18_May2016.shp'
newField = 'code'
origField = 'CODE'
rte_adj_fact = get_rate_adjust_factor(domshp, newField, origshp, origField)
              
###############################################################################
# write initial shapefile
###############################################################################

outshp = 'NSHA13_Background_NSHA18.shp'
bval_fix = -99 * ones_like(rte_adj_fact)
bval_sig_fix = -99 * ones_like(rte_adj_fact)

build_source_shape(outshp, shapes, src_names, src_codes, zone_class, \
                   rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                   min_rmag, mmax, bval_fix, bval_sig_fix, \
                   ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                   shmax_pref, shmax_sig, trt_new, domains, prefCat)


