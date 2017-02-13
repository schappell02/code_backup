from gcwork import objects
import sc_starset as starset
#import sc_accel_class as acc
#from gcwork import util
#from gcwork import orbits
import sc_young as young
#from pysqlite2 import dbapi2 as sqlite
from sqlite3 import dbapi2 as sqlite
import scipy
import pyfits
from scipy import stats
from scipy import special
from scipy import integrate
from gcwork import starTables
import pickle
import nmpfit_sy
import asciidata, os, sys, pickle
import nmpfit_sy2 as nmpfit_sy
#from pylab import *
import numpy as np
#import pylab as py
import math
#import histNofill
#import matplotlib.axes3d as p3
import pdb
import time
import scipy.optimize as opter
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter 
import pymultinest
import datetime
import time
import threading
from astropy.stats import LombScargle
# Several functions were taken from:
# /ghezgroup/code/python/gcwork/polyfit/accel.py,
# velocity.py and residuals.py, written by JLu.

root = '/g/ghez/align/'
home='/u/pi/'

pi = math.pi

# Mass and Ro from S0-2 - Ghez et al. 2008, Anna Boehle up to 2013
mass = 4.07e6 #3.93e6
masse = 0.6e6 #0.19e6
dist = 7960.0 #7.84e3
#diste = 0.14e3
G = 6.6726e-8
msun = 1.99e33
GM = G * mass * msun
mass_g = mass * msun
masse_g = masse * msun
sec_in_yr = 3.1557e7
cm_in_au = 1.496e13
cm_in_pc = 3.086e18
km_in_pc = 3.086e13
au_in_pc = 206265.0
asy_to_kms = dist * cm_in_au / (1e5 * sec_in_yr)
as_to_km = dist * cm_in_au / (1e5)
density0 = 3.2e4
density0_g = 3.2e4 * msun
density0e = 1.3e4
density0e_g = density0e * msun
GM_as_yr = GM * sec_in_yr**2 * 1e-15 / as_to_km**3


def outLombScargle(label='',min_p=1.0e-3,max_p=1.0e1):
    #min and max period in YEARS
    min_f = 2.0*pi / max_p
    max_f = 2.0*pi / min_p

    RV = np.loadtxt('/u/schappell/'+label+'date_vz_err.dat')
    date = RV[:,0]
    vz = RV[:,1]
    vzerr = RV[:,2]
    
    frequency,power = LombScargle(date,vz,vzerr).autopower(minimum_frequency=min_f,maximum_frequency=max_f)

    period = 2.0*pi/frequency

    np.savetxt('/u/schappell/'+label+'LS_period_power.dat',np.transpose([period,power]),delimiter=' ')



def chi2_sin(param,*needed):
    #pdb.set_trace()
    date_obs, RV_obs, err_obs = needed
    #chi2 = np.array([])
    #for val in param:
        #pdb.set_trace()
    amp = param[0]
    period = param[1]
    offset = param[2]
    chi2 = (RV_obs - (amp*np.sin(2.0*pi*date_obs/period + offset)))**2 / err_obs**2
        #chi2 = np.append(chi2,tmp)

    return np.sum(chi2)


def chi2_sin_curve(date_obs,amp,period,offset):
    #pdb.set_trace()
    #date_obs, RV_obs, err_obs = needed
    #chi2 = np.array([])
    #for val in param:
        #pdb.set_trace()
    #amp = param[0]
    #period = param[1]
    #offset = param[2]
    #chi2 = (RV_obs - (amp*np.sin(2.0*pi*date_obs/period + offset)))**2 / err_obs**2
        #chi2 = np.append(chi2,tmp)

    return amp*np.sin(2.0*pi*date_obs/period + offset)




def mock_fitSin(guess=[50.0,2.0/365.0,1.5]):
    mockData = np.loadtxt('/u/schappell/mockRV.dat')
    use_date = mockData[:,0]
    mockRV = mockData[:,1]
    use_err = mockData[:,2]

    results = minimize(chi2_sin,guess,method='nelder-mead',args=(use_date,mockRV,use_err))
    param = results.x
    
    print ''
    print 'Results:'
    print 'Amp = '+str(param[0])+' km/s'
    print 'Period = '+str(param[1])+' years'
    print 'Offset = '+str(param[2])+' radians'
    print ''
