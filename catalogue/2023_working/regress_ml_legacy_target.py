# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:14:25 2023

@author: u56903
"""
import pickle
from numpy import arange, array, delete, isnan, where
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import scipy.odr.odrpack as odrpack
#from scipy.odr import Model, Data, ODR
from misc_tools import dictlist2array
import matplotlib as mpl
mpl.style.use('classic')

###############################################################################
# load data
###############################################################################

# get_station_distance_stadat(sta, eqlo, eqla), 
# first, load pkl
mcdat = pickle.load(open('merged_cat_revised_ml.pkl', 'rb'))

###############################################################################
# get data arrays
###############################################################################

ml_legacy = dictlist2array(mcdat, 'PREFML')
ml_rev = dictlist2array(mcdat, 'REVML_2023')
ml_legacy_2018 = ml_legacy
ml_rev_2018 = dictlist2array(mcdat, 'REVML')
ml_target = dictlist2array(mcdat, 'targetML')

# remove null target eqn
idx = ml_target == 'null'
ml_legacy = delete(ml_legacy, idx)
ml_rev = delete(ml_rev, idx)

###############################################################################
# regress
###############################################################################
def ortho_lin_reg(c, x):
    return c[0] * x + c[1]

# exlude nan data
idx = where((isnan(ml_legacy) == False) & (isnan(ml_rev) == False))[0] 
data = odrpack.RealData(ml_legacy[idx], ml_rev[idx])
'''
sx = np.interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
data = odrpack.RealData(mb[mb>=3.5], mw[mb>=3.5], sx=sx[mb>=3.5], sy=sy[mb>=3.5])
'''
lin_reg = odrpack.Model(ortho_lin_reg)
odr = odrpack.ODR(data, lin_reg, beta0=[1.0, 0.0])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
c0 = out.beta[0]
c1 = out.beta[1]

# write regression coeffs
outtxt = ','.join((str(c0), str(c1)))
f = open('ml_revision_reg.csv', 'w')
f.write(outtxt)
f.close()

###############################################################################
# plot
###############################################################################

fig = plt.figure(1, figsize=(8,8))
ax = plt.subplot(111)

plt.plot([2,7],[2,7],'k--')
plt.plot(ml_legacy, ml_rev, '+', ms=6, c='#606f74', label='2023 Revised')
#plt.plot(ml_legacy_2018, ml_rev_2018, 'x', ms=6, c='0.75', label='2018 Rev')

# plt 2023 relationship
xplt = array([2,7])
yplt = c0 * xplt + c1
plt.plot(xplt, yplt, '-', lw=2, c='#773775', label='2023 Regression')

# plt 2018 relationship
c0 = 0.84
c1 = 0.35
'''
# 2021 JoS coeffs
c0 = 0.90
c1 = 0.09
'''
yplt = c0 * xplt + c1
plt.plot(xplt, yplt, '-', lw=2, c='#0b5e4a', label='2018 Regression')

plt.grid(which='both')
plt.xlabel('Legacy ML', fontsize=20)
plt.ylabel('Revised ML', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.xlim([2,7])
plt.ylim([2,7])
plt.legend(loc=2, numpoints=1, fontsize=16)

plt.savefig('ml_revision_reg.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
