# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:14:25 2023

@author: u56903
"""
import pickle
from numpy import arange, array, delete, isnan, where, loadtxt, concatenate, interp
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import scipy.odr.odrpack as odrpack
#from scipy.odr import Model, Data, ODR
from misc_tools import dictlist2array
import matplotlib as mpl
mpl.style.use('classic')

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x >= hx
    xmod[idx] = 1
    return xmod

def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_free(c, x):
    from numpy import zeros_like
    hx = c[3] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)

    ans1 = modx_lo * (c[0] * x + c[1])
    yhinge = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yhinge)

    return ans1 + ans2

def bilinear_reg_fix(c, x):
    from numpy import zeros_like
    hx = 4.5 #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)

    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yarea)

    return ans1 + ans2

###############################################################################
# load data
###############################################################################

# get_station_distance_stadat(sta, eqlo, eqla), 
# first, load pkl
mcdat = pickle.load(open('merged_cat_pref_mags.pkl', 'rb'))

###############################################################################
# get ML-MW arrays
###############################################################################

ms = dictlist2array(mcdat, 'PREFMS')
mw = dictlist2array(mcdat, 'PREFMW')
evdt = dictlist2array(mcdat, 'DATETIME')
msref = dictlist2array(mcdat, 'PREFMSSRC')

# remove data < 1963
idx = evdt < UTCDateTime(1950,1,1)
ms = delete(ms, idx)
mw = delete(mw, idx)
evdt = delete(evdt, idx)
msref = delete(msref, idx)

# ignore AUST data between 2010 and 2018.5
idx = where((evdt >= UTCDateTime(2010,1,1)) & (evdt <= UTCDateTime(2018,6,1)) \
            & (msref == 'AUST'))[0]
    
ms = delete(ms, idx)
mw = delete(mw, idx)
evdt = delete(evdt, idx)
msref = delete(msref, idx)

###############################################################################
# regress all corrections
###############################################################################
def ortho_lin_reg(c, x):
    return c[0] * x + c[1]

def ortho_quad_reg(c, x):
    #return c[0] * x**2 + c[1] * x + c[2]
    #return c[0] * (x-7.0)**2 + c[1] * (x-7.0) + 7.0
    return c[0] * x**2 + c[1]

# exlude nan data
idx = where((isnan(ms) == False) & (isnan(mw) == False))[0] 
data = odrpack.RealData(ms[idx], mw[idx])
'''
sx = interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
data = odrpack.RealData(ms[idx], mw[idx], sx=sx[idx], sy=sy[idx])
'''
quad_reg = odrpack.Model(ortho_quad_reg)
odr = odrpack.ODR(data, quad_reg, beta0=[0.075, 3.])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
c0 = out.beta[0]
c1 = out.beta[1]
out.pprint()

# write regression coeffs
outtxt = ','.join((str(c0), str(c1)))
f = open('mw-ms_empirical_coeffs.csv', 'w')
f.write(outtxt)
f.close()

# export ms-mw data
txt = 'DATETIME,ms,MW\n'
for i in range(0, len(idx)):
    txt += ','.join((str(evdt[idx[i]]), str(ms[idx[i]]), str(mw[idx[i]]))) + '\n'
                     
f = open('ms_mw_data_list.csv', 'w')
f.write(txt)
f.close()

"""
###############################################################################
# bi-linear regression
###############################################################################

############### bilinear auto
data = odrpack.RealData(ms[idx], mw[idx])
'''
sx = interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
'''
data = odrpack.RealData(ms[idx], mw[idx]) # , sx=sx[idx], sy=sy[idx]

bilin_reg = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data, bilin_reg, beta0=[1.0, -0.5, 1.5, 5.])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
out.pprint()

b0 = out.beta[0]
b1 = out.beta[1]
b2 = out.beta[2]
hx = out.beta[3] # x hinge point
"""
###############################################################################
# plot data
###############################################################################

fig = plt.figure(1, figsize=(8,8))
ax = plt.subplot(111)
"""
idx1 = where(mwref=='Allen et al. (2006)')[0]
idx2 = where(mwref=='Allen (2012)')[0]
idx3 = where(mwref=='Ghasemi et al (2016)')[0]
idx4 = where(mwref=='AUST')[0]
idx5 = []
for i, mr in enumerate(mwref):
    if mr != 'Ghasemi et al (2016)' and mr != 'Allen (2012)' \
       and mr != 'Allen et al. (2006)' and mr != 'AUST':
        idx5.append(i)

d1 = plt.plot(ms[idx1], mw[idx1], 'o', ms=6, c='#6e7645', alpha=1.0) #, label='Allen et al. (2006)')
d2 = plt.plot(ms[idx2], mw[idx2], '^', ms=6, c='#ca7700', alpha=1.0) #, label='Allen (2012)')
d3 = plt.plot(ms[idx3], mw[idx3], 'v', ms=6, c='#988642', alpha=1.0) #, label='Ghasemi et al. (2017)')
d4 = plt.plot(ms[idx4], mw[idx4], 'd', ms=6, c='#a5d867', alpha=1.0) #, label='AUST')
d5 = plt.plot(ms[idx5], mw[idx5], 'h', ms=6, c='#003145', alpha=1.0) #, label='Other')

leg1 = plt.legend([d1[0], d2[0], d3[0], d4[0], d5[0]], \
                  ['Allen et al. (2006)', 'Allen (2012)', 'Ghasemi et al. (2017)', 'AUST', 'Other'], \
                  loc=4, numpoints=3, fontsize=14)
"""
plt.plot(ms[idx], mw[idx], 'o', ms=6, c='#606f74', alpha=0.5, label='Data')
#plt.plot(ml_legacy_2018, ml_rev_2018, 'x', ms=6, c='0.75', label='2018 Rev')

###############################################################################
# plot regression
###############################################################################
plt.plot([2,7],[2,7],'k--', label='1:1')

# plt 2023 relationship
xplt = arange(2, 7, 0.01)
#yplt = c0 * (xplt-7.0)**2 + c1 * (xplt-7.0) + 7.0
yplt = c0 * xplt**2 + c1
plt.plot(xplt, yplt, '-', lw=2, c='#773775', label='2023 Quadratic')

'''
yplt = b1 + b0 * xplt
yhinge = b1 + b0 * hx
idx = xplt > hx
yplt[idx] = b2 * (xplt[idx]-hx) + yhinge
plt.plot(xplt, yplt, '-', lw=2, c='#773775', label='2023 Bi-linear')
'''
###############################################################################
# plot other models
###############################################################################

# ################## TA-model
mw_ta = 0.060977193461713*(xplt-7.)**2 + 0.988674360853538*(xplt-7.) + 7.
# ################## Johnston
log_mw_j = 24.66 - 1.083*xplt + 0.192*xplt**2
mw_j = 2*log_mw_j / 3 - 10.7
# ################## Scordilis (2006)
mw_s = concatenate((0.67 * xplt[xplt<=6.1]+2.07,
                        0.99 * xplt[xplt>6.1] + 0.08))
# ################## Das etal 2010
mw_das = concatenate((0.67 * xplt[xplt<=6.1]+2.12,
                        1.06 * xplt[xplt>6.1] - 0.38))
#mw_das = 1.538*xplt - 2.538
# ################## Youngs (2012)
mw_yon = 2.654 + 0.334*xplt + 0.040*xplt**2

# ################## Di Giacomo (2015)
mw_dig = concatenate((0.67 * xplt[xplt<=6.47]+2.13,
                        1.10 * xplt[xplt>6.47] - 0.67))

'''
# plt 2023 simulated (2800) relationship
s0, s1, s2 = loadtxt('mw-ms_coeffs_2800.csv', delimiter=',', skiprows=1)
yplt = s0 * xplt**2 + s1 * xplt + s2
plt.plot(xplt, yplt, '-', lw=2, c='#0b5e4a', label='NSHA23')
'''
# plt 2018 relationship
c0 = 0.075
c1 = 3.357
yplt = c0 * xplt**2 + c1
plt.plot(xplt, yplt, '-', lw=2, c='#637c6b', label='NSHA18')

plt.plot(xplt,mw_j,'-', c='#72c7e7', lw=2,label='Johnston (1996)')
plt.plot(xplt,mw_yon,'-', c='#0b5e4a', lw=2,label='Youngs (2012)')
#plt.plot(xplt,mw_ta,'-', c='#cb6c37',lw=2,label='Allen (2012)')
plt.plot(xplt,mw_dig,'-', c='#b43b3b',lw=2,label='Di Giacomo et al. (2015)')

plt.grid(which='both')
plt.xlabel(r'$\mathregular{M_{S}}$', fontsize=20)
plt.ylabel(r'$\mathregular{M_{W}}$', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.xlim([3,7])
plt.ylim([3,7])
plt.legend(loc=2, numpoints=3, fontsize=12)
#plt.gca().add_artist(leg1)

plt.savefig('mw-ms.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
