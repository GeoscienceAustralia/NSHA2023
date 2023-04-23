# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:14:25 2023

@author: u56903
"""
import pickle
from numpy import arange, array, delete, isnan, where, loadtxt, zeros_like, hstack
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import scipy.odr.odrpack as odrpack
#from scipy.odr import Model, Data, ODR
from misc_tools import dictlist2array, get_mpl2_colourlist
import matplotlib as mpl
mpl.style.use('classic')
cols = get_mpl2_colourlist()

###############################################################################
# load data
###############################################################################

# get_station_distance_stadat(sta, eqlo, eqla), 
# first, load pkl
mcdat = pickle.load(open('merged_cat_pref_mags.pkl', 'rb'))

###############################################################################
# get ML-MW arrays
###############################################################################

ml_2800 = dictlist2array(mcdat, 'PREFML_2023')
mw = dictlist2array(mcdat, 'PREFMW')
evdt = dictlist2array(mcdat, 'DATETIME')
mwref = dictlist2array(mcdat, 'PREFMWSRC')

# remove data < 1963
idx = evdt < UTCDateTime(1990,1,1)
ml_2800 = delete(ml_2800, idx)
mw = delete(mw, idx)
evdt = delete(evdt, idx)
mwref = delete(mwref, idx)

###############################################################################
# add data from recent events
###############################################################################
# events: 2021 Marble Bar, 2022 Arthur River
newml = array([5.55+0.13, 4.753+0.13, 5.0+0.13]) # last one just a guess
newmw = array([5.323, 4.5, 4.69268914])
newdt = array([UTCDateTime('2021-11-13T13:05:52.663'), UTCDateTime('2022-01-24T21:24:47.666'), \
              UTCDateTime('2023-04-17T09:02:57.386')])
newref = array(['AUST', 'AUST', 'AUST'])

# load W-A correction coeffs
fn = loadtxt('wa_sensitivity_coeffs.csv', delimiter=',', skiprows=1)

# conv 2080 to 2800
norm_mean = fn[5]
norm_std  = fn[6]
x_norm = (newml - norm_mean)/norm_std
ML2800_corr = fn[0] * x_norm**4 + fn[1] * x_norm**3 + fn[2] * x_norm**2 + fn[3] * x_norm + fn[4]
newml += ML2800_corr

# concat
ml_2800 = hstack((ml_2800, newml))
mw = hstack((mw, newmw))
evdt = hstack((evdt, newdt))
mwref = hstack((mwref, newref))

###############################################################################
# regress all corrections
###############################################################################
def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x > hx
    xmod[idx] = 1
    return xmod

def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x <= hx
    xmod[idx] = 1
    return xmod

def ortho_lin_reg(c, x):
    return c[0] * x + c[1]

def ortho_quad_reg(c, x):
    return c[0] * x**2 + c[1] * x + c[2]

def lin2quad_reg(c, x):
   from numpy import zeros_like
   hx = c[4] # hinge magnitude
   ans2 = zeros_like(x)
   ans1 = zeros_like(x)

   modx_lo = lowside(x, hx)
   modx_hi = highside(x, hx)
   
   ans1 = modx_lo * (c[0] * x + c[1])
   yhinge = c[0] * hx + c[1]
   ans2 = modx_hi * (c[2]*(x-hx) + c[3]*(x-hx)**2 + yhinge)

   return ans1 + ans2

# exlude nan data
idx = where((isnan(ml_2800) == False) & (isnan(mw) == False) \
             & (mwref != 'Ghasemi et al (2016)'))[0] 
data = odrpack.RealData(ml_2800[idx], mw[idx])
'''
sx = interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
data = odrpack.RealData(mb[mb>=3.5], mw[mb>=3.5], sx=sx[mb>=3.5], sy=sy[mb>=3.5])
'''
# do straight quadratic
quad_reg = odrpack.Model(ortho_quad_reg)
odr = odrpack.ODR(data, quad_reg, beta0=[0.09, 0.1, 2.0])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
c0 = out.beta[0]
c1 = out.beta[1]
c2 = out.beta[2]

# write regression coeffs
outtxt = ','.join((str(c0), str(c1), str(c2)))
f = open('mw-ml_empirical_coeffs.csv', 'w')
f.write(outtxt)
f.close()

#################
# do lin2quad

l2q_reg = odrpack.Model(lin2quad_reg)
odr = odrpack.ODR(data, l2q_reg, beta0=[0.6, 1.0, 0.7, 0.1, 4.5])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
lqc0 = out.beta[0]
lqc1 = out.beta[1]
lqc2 = out.beta[2]
lqc3 = out.beta[3]
lqc4 = out.beta[4]

'''
# first get linear section
idx = where((isnan(ml_2800) == False) & (isnan(mw) == False) \
             & (ml_2800 < 4.0) & (mwref != 'Ghasemi et al (2016)'))[0] 
l2q_reg = odrpack.Model(ortho_lin_reg)
odr = odrpack.ODR(data, l2q_reg, beta0=[0.67, 1.0])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
lqc0 = out.beta[0]
lqc1 = out.beta[1]
'''

# write regression coeffs
outtxt = ','.join((str(lqc0), str(lqc1), str(lqc2), str(lqc3), str(lqc4)))
f = open('mw-ml_lin2quad_coeffs.csv', 'w')
f.write(outtxt)
f.close()

# export ml-mw data
txt = 'DATETIME,ML,MW\n'
for i in range(0, len(idx)):
    txt += ','.join((str(evdt[idx[i]]), str(ml_2800[idx[i]]), str(mw[idx[i]]))) + '\n'
                     
f = open('ml_mw_data_list.csv', 'w')
f.write(txt)
f.close()
###############################################################################
# plot data
###############################################################################

fig = plt.figure(1, figsize=(8,8))
ax = plt.subplot(111)

idx1 = where(mwref=='Allen et al. (2006)')[0]
idx2 = where(mwref=='Allen (2012)')[0]
idx3 = where(mwref=='Ghasemi et al (2016)')[0]
idx4 = where(mwref=='AUST')[0]
idx5 = where(mwref=='GCMT')[0]
idx6 = []
for i, mr in enumerate(mwref):
    if mr != 'Ghasemi et al (2016)' and mr != 'Allen (2012)' \
       and mr != 'Allen et al. (2006)' and mr != 'AUST' and mr != 'GCMT':
        idx6.append(i)

d1 = plt.plot(ml_2800[idx1], mw[idx1], 'o', ms=6, c=cols[0], alpha=1.0) #, label='Allen et al. (2006)')
d2 = plt.plot(ml_2800[idx2], mw[idx2], '^', ms=6, c=cols[1], alpha=1.0) #, label='Allen (2012)')
d3 = plt.plot(ml_2800[idx3], mw[idx3], 'v', ms=6, c=cols[2], alpha=1.0) #, label='Ghasemi et al. (2017)')
d4 = plt.plot(ml_2800[idx4], mw[idx4], 'D', ms=6, c=cols[3], alpha=1.0) #, label='AUST')
d5 = plt.plot(ml_2800[idx5], mw[idx5], 's', ms=6, c=cols[4], alpha=1.0) #, label='Other')
d6 = plt.plot(ml_2800[idx6], mw[idx6], 'h', ms=6, c=cols[5], alpha=1.0) #, label='Other')

leg1 = plt.legend([d1[0], d2[0], d3[0], d4[0], d5[0], d6[0]], \
                  ['Allen et al. (2006)', 'Allen (2012)', 'Ghasemi et al. (2017)', 'AUST', 'GCMT', 'Other'], \
                  loc=4, numpoints=3, fontsize=14)

#plt.plot(ml_2800[idx], mw[idx], 'o', ms=6, c='#606f74', alpha=0.5, label='Data')
#plt.plot(ml_legacy_2018, ml_rev_2018, 'x', ms=6, c='0.75', label='2018 Rev')

###############################################################################
# plot models
###############################################################################
plt.plot([2,7],[2,7],'k--', label='1:1')

# plt 2023 quadratic relationship
xplt = arange(2, 7, 0.01)
yplt = c0 * xplt**2 + c1 * xplt + c2
plt.plot(xplt, yplt, '-', lw=2, c='#773775', label='2023 Empirical')

# plt 2023 lin2quad relationship
'''
ans1 = modx_lo * (c[0] * x + c[1])
yhinge = c[0] * hx + c[1]
ans2 = modx_hi * (c[2] * (x-hx) + c[3] * (x-hx)**2 + yhinge)
'''
yplt = zeros_like(xplt)
idx = xplt <= lqc4
yplt[idx] = lqc0 * xplt[idx] + lqc1
idx = xplt > lqc4
yplt[idx] = lqc0 * lqc4 + lqc1 \
            + lqc2 * (xplt[idx]-lqc4) + lqc3 * (xplt[idx]-lqc4)**2
plt.plot(xplt, yplt, '-', lw=2, c='dodgerblue', label='2023 Linear-Quadratic')

# plt 2023 simulated (2800) relationship
s0, s1, s2 = loadtxt('mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)
yplt = s0 * xplt**2 + s1 * xplt + s2
plt.plot(xplt, yplt, '-', lw=2, c='#0b5e4a', label='2023 Simulated (W-A 2800)')

# plt 2018 relationship
c0 = 0.042
c1 = 0.481
c2 = 1.395
yplt = c0 * xplt**2 + c1 * xplt + c2
plt.plot(xplt, yplt, '-', lw=2, c='#cb6c37', label='NSHA18 (W-A 2080)')

plt.grid(which='both')
plt.xlabel(r'$\mathregular{M_{L(2800)}}$', fontsize=20)
plt.ylabel(r'$\mathregular{M_{W}}$', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.xlim([2,7])
plt.ylim([2,7])
plt.legend(loc=2, numpoints=3, fontsize=14)
plt.gca().add_artist(leg1)

plt.savefig('mw-ml.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
