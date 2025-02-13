# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:07:46 2015

@author: tallen
"""
# https://gis.stackexchange.com/questions/130199/changing-color-of-raster-images-based-on-their-data-values-gdal
#from pyrate.shared import Ifg, DEM
#import subprocess
#print('\nUse: conda activate gdal\n')
import sys
import os
import tempfile
import numpy as np

def gen_color_file(input_file):
    fp, temp_file = tempfile.mkstemp(suffix='.txt')

    dem = DEM(input_file)
    dem.open()
    phase_data = dem.height_band.ReadAsArray()

    max_ph = np.nanmax(phase_data)
    min_ph = np.nanmin(phase_data)
    range_ph = max_ph-min_ph
    colors = ['black', 'blue', 'yellow', 'orange', 'red', 'white']
    with open(temp_file, 'w') as f:
        for i, c in enumerate(colors[:-1]):
            f.write(str(int(min_ph + (i + 1)*range_ph/len(colors))) +
                    ' ' + c + '\n')
        f.write(str(int(max_ph - range_ph/len(colors))) +
                ' ' + colors[-1] + '\n')
    os.close(fp)
    return temp_file
'''
cmd = "gdaldem color-relief " + input_file \
      + ' ' + color_file + ' ' + output_file
      subprocess.check_call(cmd, shell=True)
'''    
'''
# read sol file
'''
from sys import argv
from os import sep, path, mkdir, system, getcwd
from numpy import array, mgrid, nan, shape, hstack, isinf, log, exp, interp
from scipy.interpolate import griddata
from osgeo import osr, gdal
import shapefile	
from shapely.geometry import Point, Polygon
from tools.oq_tools import return_annualised_haz_curves
from misc_tools import cmap2rgb, remove_last_cmap_colour
from gmt_tools import cpt2colormap
import warnings
warnings.filterwarnings("ignore")

#from misc_tools import dictlist2array

maxlat = -90
minlat = 90
maxlon = -180
minlon = 180

'''
# THIS IS HOW TO RUN ME IN IPYTHON

run csv2geoTIFF.py ../../source_models/complete_model/final/results_maps_SA10/hazard_curve-mean-SA(1.0)_1.csv
'''
##############################################################################
# parse hazard grid
##############################################################################

hazCurveGridFile = argv[1]

# check to see if geotiff folder exists
if path.isdir('geotiff') == False:
    mkdir('geotiff')
    
# make out geoTIFF    
period = hazCurveGridFile.split(sep)[-2].split('_')[-3]

siteClass = 'SC_' + hazCurveGridFile.split(sep)[-2].split('_')[-1]

# parse grid file
gridDict, imls, investigation_time = return_annualised_haz_curves(hazCurveGridFile)

##############################################################################
# interpolate hazard curve and fill dictionary
##############################################################################

# set grid return periods for hazard curve
probs = array([0.02,0.01375,0.01,0.00445,0.002,0.0021,0.001,0.0005,0.000671,0.000404,0.0002,0.0001])

grddict = []
alon = []
alat = []

for site in gridDict:
    poe ='poe' # short term fix
    #interpHaz = exp(interp(log(probs[::-1]), log(site[period+'_probs_annual'][::-1]), log(imls[::-1])))[::-1]
    interpHaz = exp(interp(log(probs[::-1]), log(site[poe+'_probs_annual'][::-1]), log(imls[::-1])))[::-1]
    
    # fill a temp dictionary
    tmpdict = {'lon':site['lon'], 'lat':site['lat']}
    
    for p, ih in zip(probs, interpHaz):
        if ih < 1E-20:
            tmpdict['P'+str(p)] = 1E-20
        else:
            tmpdict['P'+str(p)] = ih
    
    grddict.append(tmpdict)
    alon.append(tmpdict['lon'])
    alat.append(tmpdict['lat'])
    
    # get bbox 
    if tmpdict['lon'] > maxlon:
        maxlon = tmpdict['lon']
    if tmpdict['lon'] < minlon:
        minlon = tmpdict['lon']
        
    if tmpdict['lat'] > maxlat:
        maxlat = tmpdict['lat']
    if tmpdict['lat'] < minlat:
        minlat = tmpdict['lat']

print('BBOX', '/'.join((str(minlon), str(maxlon), str(minlat), str(maxlat))))
##############################################################################
# make mesh
##############################################################################

resolution = 0.05 # degrees
invres = int(1./resolution)
keys = ['P0.0021', 'P0.000671', 'P0.000404'] # probabilities
pc50 = ['0.1', '0.033', '0.02']
for key, p50 in zip(keys, pc50):
    print('Making', key, 'grid mesh...')
    """
    # first make z data array
    ahaz = []
    for grdval in grddict:
        ahaz.append(grdval[key])
    
    # space at define resolution    
    nlon = int(round((maxlon-minlon)*invres))
    nlat = int(round((maxlat-minlat)*invres))
    grid_x, grid_y = mgrid[minlon:maxlon:complex(nlon), minlat:maxlat:complex(nlat)]  # lon: 1000 @ 0.1 deg; lat: 500 @ 0.1 deg
    
    # interpolate z data array
    grid_z = griddata((array(alon), array(alat)), array(ahaz), (grid_x, grid_y), \
                      method='cubic', fill_value=nan) # scipy.interpolate (linear, nearest, cubic)
    #grid_z = griddata(array(alon), array(alat), array(ahaz), grid_x, grid_y, interp='linear') # matplotlib
    
    
    # mask grid points outside defined grid to avoid extrapolation
    print('Masking', key, 'grid...')
    if getcwd().startswith('/nas'):
        inshape = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/shapefiles//au_maritime_boundary_digitised.shp'
    else:
        inshape = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/postprocessing/maps/shapefiles/au_maritime_boundary_digitised.shp'
    
    sf = shapefile.Reader(inshape)
    sf = sf.shapes()
    poly = Polygon(sf[0].points)
    
    flat_x = grid_x.flatten()
    flat_y = grid_y.flatten()
    flat_z = grid_z.flatten()
    
    # now loop through points
    print('Cropping points...')
    for i in range(0, len(flat_x)):
        point = Point(flat_x[i], flat_y[i])
        if point.within(poly) == False:
            flat_z[i] = nan
    
    # reshape xyz
    grid_z = flat_z.reshape(shape(grid_z))
    grid_z = grid_z[:,::-1]
    
    '''
    write mesh
    '''
    
    print('Writing', key, 'geoTIFF...')
    #output_file = path.join('geotiff', '_'.join(('nsha18',period, p50+'.tiff')))
    output_file = path.join('geotiff', '_'.join(('nsha23',period, siteClass, p50+'.tiff')))
    
    # Create gtif
    nbands = 1
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_file, nlon, nlat, nbands, gdal.GDT_Float32)
    
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    #dst_ds.SetGeoTransform([ minlon, 0.1, 0, minlat, 0, 0.1 ] ) - this works, but map is upsidedown
    dst_ds.SetGeoTransform([ minlon, resolution, 0, maxlat, 0, -resolution])
      
    # set the reference info 
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS("WGS84") # alt NAD27, NAD83
    #srs.SetLCC() # to set Lambert Conformal Conic
    dst_ds.SetProjection(srs.ExportToWkt())
    
    # write the band
    dst_ds.GetRasterBand(1).WriteArray(grid_z.T)
    dst_ds = None # to close file

    # testing
    #src_ds = gdal.Open(path.join('geotiff', '_'.join(('nsha18',period, p50+'.tiff'))))
    """
    ##############################################################################
    # make gdal cpt file and recolour
    ##############################################################################
    
    # set bounds for colours
    if p50 == '0.1' or p50 == '0.033':

        if period == 'PGA':
            if p50 == '0.033':
                bounds = array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.14, 0.20, 0.26, 0.32, 0.4])
            else:
                bounds = array([0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24])
        elif period == 'SA005' or period == 'SA01'  \
           or period == 'SA03' or period == 'SA05':
            bounds = array([0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24])
        elif  period == 'SA07'  or period == 'SA10':
            bounds = array([0, 0.002, 0.004, 0.007, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.08])
        elif period == 'SA02':
            if p50 == '0.033':
                bounds = array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.26, 0.36, 0.48, 0.6])
            else:
                bounds = array([0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24])
        else:
            bounds = array([0, 0.001, 0.002, 0.004, 0.007, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.045, 0.06])
    else:
        if period == 'PGA' or period == 'SA005' or period == 'SA01'  \
           or period == 'SA03' or period == 'SA05':
            bounds = array([0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.12, 0.16, 0.22, 0.30, 0.40, 0.5])
        elif period == 'SA02':
            bounds = array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.26, 0.36, 0.48, 0.6])
        elif period == 'SA07':
            bounds = array([0, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24, 0.36])
        elif period == 'SA15' or period == 'SA10':
            bounds = array([0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.2, 0.3])
        else:
            bounds = array([0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.015, 0.03, 0.04, 0.05, 0.06, 0.1, 0.16])
    ncolours = 13
    
       
    if getcwd().startswith('/nas'):
        cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/cw1-013_mod.cpt'
    else:
        cptfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/postprocessing/maps/cw1-013_mod.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncolours, rev=True)
    rgbTable = cmap2rgb(cmap, ncolours)[0] * 255
    
    # make gdal cpt file
    cpttxt = ''
    i = 0
    for bound, rgb in zip(bounds[:-1], rgbTable[:-1]):
        cpttxt += '\t'.join((str(bound), str(rgb[0]), str(rgb[1]), str(rgb[2]))) + '\n'
        cpttxt += '\t'.join((str(bounds[i+1]-0.0000001), str(rgb[0]), str(rgb[1]), str(rgb[2]))) + '\n'
        i += 1
    cpttxt += '\t'.join((str(bounds[-1]), str(rgbTable[-1][0]), str(rgbTable[-1][1]), str(rgbTable[-1][2]))) + '\n'
    cpttxt += '\t'.join(('nv','0','0','0','0'))
    
        
    f = open('gdal_cpt.dat', 'w')
    f.write(cpttxt)
    f.close()
    
    ##############################################################################
    # write geotiff
    ##############################################################################
    
    # recolour geoTiff
    intiff = path.join('geotiff', '_'.join(('nsha23',period, siteClass, p50+'.tiff')))
    outtiff = intiff[0:-5]+'_colour.tiff'
    		
    '''
    # set nan values to no data
    gdaledit_cmd = ' '.join(('/Users/trev/opt/miniconda3/envs/py37/bin/gdal_edit.py -a_nodata -nan',intiff))
    system(gdaledit_cmd)
    
    gdaldem_cmd = ' '.join(('gdaldem color-relief', intiff, 'gdal_cpt.dat', outtiff))
    #gdaldem color-relief jotunheimen.tif color_relief.txt jotunheimen_colour_relief.tif
    system(gdaldem_cmd)
    '''
    
    ##############################################################################
    # write qml
    ##############################################################################
    
    qmlfile = 'template.qml'
    qlines = open(qmlfile).readlines()
    
    newqml = ''
    for ql in qlines:
       if ql.strip().startswith('classificationMin='):
           newline = '"'.join(('<rasterrenderer classificationMin=', str(bounds[1]),' classificationMax=',str(bounds[-1]),' alphaBand="-1" nodataColor="" band="1" opacity="1" type="singlebandpseudocolor">'))+'\n'
           newqml += newline
       elif ql.strip().startswith('<colorrampshader minimumValue='):
           newline = '"'.join(('<colorrampshader minimumValue=', str(bounds[1]+1E-19),' labelPrecision="4" maximumValue=',str(bounds[-1]-1E-19),' classificationMode="1" clip="0" colorRampType="INTERPOLATED">'))+'\n'
           newqml += newline
       elif ql.strip().startswith('</colorramp>'):
           newqml += '</colorramp>\n'
           newqml += '"'.join(('<item color="#a9c8e4" alpha="255" value=', str(bounds[1]), ' label=', str('%0.4f' % bounds[1]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#8daac6" alpha="255" value=', str(bounds[2]), ' label=', str('%0.4f' % bounds[2]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#728caa" alpha="255" value=', str(bounds[3]), ' label=', str('%0.4f' % bounds[3]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#4f6581" alpha="255" value=', str(bounds[4]), ' label=', str('%0.4f' % bounds[4]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#8d9971" alpha="255" value=', str(bounds[5]), ' label=', str('%0.4f' % bounds[5]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#c5c356" alpha="255" value=', str(bounds[6]), ' label=', str('%0.4f' % bounds[6]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#fbee3b" alpha="255" value=', str(bounds[7]), ' label=', str('%0.4f' % bounds[7]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#f1ba3b" alpha="255" value=', str(bounds[8]), ' label=', str('%0.4f' % bounds[8]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#e8853c" alpha="255" value=', str(bounds[9]), ' label=', str('%0.4f' % bounds[9]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#de513c" alpha="255" value=', str(bounds[10]), ' label=', str('%0.4f' % bounds[10]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#b44b3f" alpha="255" value=', str(bounds[11]), ' label=', str('%0.4f' % bounds[11]), '/>'))+'\n'
           newqml += '"'.join(('<item color="#8a4343" alpha="255" value=', str(bounds[12]), ' label=', str('%0.4f' % bounds[12]), '/>'))+'\n'
       else:
           newqml += ql
           
    # now write to file
    newqmlfile = intiff[:-4] + 'qml'
    f = open(newqmlfile, 'w')
    f.write(newqml)
    f.close()

