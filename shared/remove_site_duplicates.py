from numpy import loadtxt

dat = loadtxt('nsha23_map_sites.csv', delimiter=',')

newtxt = ','.join((str(dat[0][0]), str(dat[0][1]))) + '\n'

for i, d in enumerate(dat[1:]):
	# here i is previous index
	if d[0] == dat[i][0] and d[1] == dat[i][1]:
		print('skipping', d)
	else:
		newtxt += ','.join((str(d[0]), str(d[1]))) + '\n'
		
f = open('nsha23_map_sites.hold', 'w')
f.write(newtxt)
f.close()