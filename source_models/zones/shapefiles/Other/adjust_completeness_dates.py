import shapefile
import pprint as p
from mapping_tools import get_field_data, get_field_index


shpfile = 'gridded_polygons_3d_completeness.shp'

print('Reading source shapefile...')
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
recs = sf.records()
fields = sf.fields[1:]

# get completeness years
ycomp = get_field_data(sf, 'YCOMP', 'str')
yidx = get_field_index(sf, 'YCOMP')

mcomp = get_field_data(sf, 'MCOMP', 'str')
midx = get_field_index(sf, 'MCOMP')

# loop thru recs
for i, rec in enumerate(recs):

    if rec[yidx] == '1985;1970;1960;1905;1880' and rec[midx] == '3.05;4.05;4.55;6.05;6.45':
        recs[i][yidx] = '1985;1970;1960;1930;1905;1880'
        recs[i][midx] = '3.05;4.05;4.55;5.55;6.05;6.45'

    if rec[yidx] == '2000;1980;1964;1900' and rec[midx] == '3.4;3.5;5.0;6.0':
        recs[i][yidx] = '2000;1980;1964;1930;1900'
        recs[i][midx] = '3.45;4.55;5.05;5.55;6.05'
        
    if rec[yidx] == '1990;1980;1964;1900' and rec[midx] == '3.0;3.5;5.0;6.0':
        recs[i][yidx] = '2000;1980;1964;1930;1900'
        recs[i][midx] = '3.45;4.55;5.05;5.55;6.05'
        
    if rec[yidx] == '1995;1970;1960;1905;1880' and rec[midx] == '3.05;3.95;4.55;6.05;6.55':
        recs[i][yidx] = '1995;1970;1960;1930;1905;1880'
        recs[i][midx] = '3.05;3.95;4.55;5.55;6.05;6.55'
        
    if rec[yidx] == '1980;1970;1960;1930;1900;1880' and rec[midx] == '3.05;3.95;4.45;5.45;6.05;6.45':
        recs[i][yidx] = '1980;1970;1960;1930;1900;1880'
        recs[i][midx] = '3.05;3.95;4.45;5.45;6.05;6.45'
        
    if rec[yidx] == '1985;1975;1960;1905;1880' and rec[midx] == '3.05;3.55;4.55;6.05;6.45':
        recs[i][yidx] = '1985;1975;1960;1930;1905;1880'
        recs[i][midx] = '3.05;3.55;4.55;5.55;6.05;6.45'
        
    if rec[yidx] == '1980;1970;1960;1950;1885' and rec[midx] == '2.95;3.45;4.05;4.55;5.65':
        recs[i][yidx] = '1980;1970;1960;1950;1910;1885'
        recs[i][midx] = '3.05;3.45;4.05;4.55;5.65;6.45'


# re-write shape
outshp = shpfile[:-4]+'_adj.shp'
w = shapefile.Writer(outshp, shapeType=5)

# set fields
for field in fields:
    w.field(field[0], field[1], field[2], field[3])

for i, shape in enumerate(shapes):

    # set shape polygon
    w.poly([shape.points])
    
    # set record
    w.record(recs[i][0], recs[i][1], recs[i][2], recs[i][3], recs[i][4], recs[i][5], recs[i][6], recs[i][7], recs[i][8], recs[i][9], recs[i][10], \
             recs[i][11], recs[i][12], recs[i][13], recs[i][14], recs[i][15], recs[i][16], recs[i][17], recs[i][18], recs[i][19], recs[i][20], \
             recs[i][21], recs[i][22], recs[i][23], recs[i][24], recs[i][25], recs[i][26], recs[i][27], recs[i][28], recs[i][29], recs[i][30], \
             recs[i][31], recs[i][32], recs[i][33], recs[i][34], recs[i][35], recs[i][36])
    
w.close()

# write projection file
print(outshp)
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'w')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()








