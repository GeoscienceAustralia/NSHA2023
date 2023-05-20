from shapely.geometry import Point, Polygon
import shapefile
from numpy import array
from mapping_tools import distance, reckon, get_line_parallels

# set outer bbox
minlon = 109
maxlon = 158
minlat = -47
maxlat = -7.5

shpfile = '../source_models/zones/shapefiles/Other/hazard_crop.shp'
faultshp = '../source_models/faults/FSM/FSD_simple_faults_plus_vert2.shp'
#csvfile = '../source_models/zones/2018_mw/Domains_multi_mc/results_maps_PGA/hazard_map-mean_mapping.csv'

# build grid
step = 12.5 # km
lons = []
lats = []

lat = minlat
while lat < maxlat:
    lon = minlon
    lons.append(lon)
    lats.append(lat)
    
    while lon < maxlon:
        
        # iterate lon
        lon, lat_ignore = reckon(lat, lon, step, 90.)
        
        # add arrays
        lons.append(lon)
        lats.append(lat)
        
    # iterate lat
    lon, lat = reckon(lat, minlon, step, 0.)        
        
# add fault sources and buffers
print('Reading fault shapefile...')
sf = shapefile.Reader(faultshp)
shapes = sf.shapes()
for shape in shapes:
   for point in shape.points:
       lons.append(point[0])
       lats.append(point[1])
       
   # add buffer - 2 km
   posazpts, negazpts = get_line_parallels(shape.points, 2.0)
   for pap, nap in zip(posazpts, negazpts):
       lons.append(pap[0])
       lats.append(pap[1])
       lons.append(nap[0])
       lats.append(nap[1])
       
   # add buffer - 5 km
   posazpts, negazpts = get_line_parallels(shape.points, 5.0)
   for pap, nap in zip(posazpts, negazpts):
       lons.append(pap[0])
       lats.append(pap[1])
       lons.append(nap[0])
       lats.append(nap[1])

lons = array(lons)
lats = array(lats)

# down sample lo/la for outsite AU polygon
dlon = lons[range(0, len(lons), 21)]
dlat = lats[range(0, len(lats), 21)]
    
#parese shapefile
print('Reading source shapefile...')
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
poly = Polygon(shapes[0].points)

print('Looping thru sites...')
# loop through points
inlo = []
inla = []
outtxt = ''

for lo, la in zip(lons, lats):
    pt = Point(lo, la)
    
    # if in poly, or every tenth point
    if pt.within(poly):
       
       outtxt += ','.join((str('%0.4f' % lo), str('%0.4f' % la))) + '\n'
       
# now add downsampled outside
for lo, la in zip(dlon, dlat):
    pt = Point(lo, la)
    
    # if in poly, or every tenth point
    if pt.within(poly) == False:
       
       outtxt += ','.join((str('%0.4f' % lo), str('%0.4f' % la))) + '\n'


# write to file
f = open('nsha23_map_sites.csv', 'w')
f.write(outtxt)
f.close()