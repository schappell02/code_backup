from gcwork import objects
import sc_starset as starset
import sc_accel_class as acc
from gcwork import util
from gcwork import orbits
import sc_young as young
from pysqlite2 import dbapi2 as sqlite
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
from pylab import *
import numpy as np
import pylab as py
import math
import histNofill
#import matplotlib.axes3d as p3
import pdb
import time
import scipy.optimize as opter
from scipy.optimize import fsolve
#from scipy.optimize import minimize
from matplotlib.ticker import ScalarFormatter 
#import pymultinest
import datetime
import time
import threading
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


def loadPop(alnDir='13_08_21/',
            align = 'align/align_d_rms_1000_abs_t', root_tmp='/g/ghez/align/',
            poly='polyfit_c/fit',points='points_c/', starlist='yngstars',to_debug=False):
    """Load a specific population of stars.

    Loads align data, omitting stars (according to chosen starlist) 
    Set starlist='all' to keep young stars in list, while still
    making LGS cuts. Set starlist='oldstars' to only include old stars
    (actually, these are non-known-young stars). Set
    starlist='yngstars' to only include known young stars.
    """

    starlists = ['all','oldstars','yngstars']

    if not(starlists.__contains__(starlist)):
        print "Star list must be 'all', 'oldstars', or 'yngstars'"
        print "Using default starlist='oldstars'"
        starlist='yngstars'

    s = starset.StarSet(root_tmp + alnDir + align)
    s.loadPolyfit(root_tmp + alnDir + poly, accel=0, arcsec=1)
    s.loadPolyfit(root_tmp + alnDir + poly, accel=1, arcsec=1)
    s.loadPoints(root_tmp + alnDir + points)

    # Get the most up-to-date young stars from Tuan's sql data base.
    yng = young.youngStarNames()

    # Pull out old stars
    yngstars = []
    names = [star.name for star in s.stars]


#    print 'All the stars'
#    print names
    for name in yng:
	try:
	    idx = names.index(name)
	    star = s.stars[idx]
	    yngstars.append(star)
	except ValueError:
	    continue

    oldstars = []

    if (starlist=='all' or starlist=='oldstars'):
        # Include all stars
        for star in s.stars:
	    oldstars.append(star)

    if starlist=='oldstars':
        # Remove young stars
        for star in yngstars:
            try:
   	        oldstars.remove(star)
            except ValueError:
                continue

    if starlist=='yngstars':
        for star in yngstars:
            r = np.sqrt(star.x**2+star.y**2)
            # Don't include the central arcsecond young stars:
            #if r > 0.8:
	    #    oldstars.append(star)
	    oldstars.append(star)
                # This is called oldstars but it's now the young stars
         

    # Remove Sgr A from list
    try:
        idx = names.index('SgrA')
        star = s.stars[idx]
        oldstars.remove(star)
    except ValueError:
        pass

#    # Remove any stars with < 3 points, as accelerations
#    # cannot be fit to these anyway
#    for star in oldstars:
#        if star.pointsCnt < 3:
#            oldstars.remove(star)

    # Build oldstars 
    for star in oldstars:
	#star.x = star.fitXv.p
	#star.y = star.fitYv.p
	#star.vx = star.fitXv.v * asy_to_kms
	#star.vy = star.fitYv.v * asy_to_kms
	#star.xerr = star.fitXv.perr
	#star.yerr = star.fitYv.perr
	#star.vxerr = star.fitXv.verr * asy_to_kms
	#star.vyerr = star.fitYv.verr * asy_to_kms
	#star.xchi2 = star.fitXv.chi2
	#star.ychi2 = star.fitYv.chi2
	#star.xchi2r = star.fitXv.chi2red
	#star.ychi2r = star.fitYv.chi2red

        # CAUTION! A star may have a velocity fit but NOT
        # have an accel fit (in the case of having only 3 points)
        # Similarly, a star may be in the align files, but not in
        # polyfit files (fit.linearFormal or fit.accelFormal)
        try:
            # In arcsec, mas/yr, mas/yr^2
            star.x0 = star.fitXa.p
            star.y0 = star.fitYa.p
            star.x0e = star.fitXa.perr
            star.y0e = star.fitYa.perr
            star.vx = star.fitXa.v
            star.vy = star.fitYa.v
            star.vxe = star.fitXa.verr
            star.vye = star.fitYa.verr
            star.ax = star.fitXa.a
            star.ay = star.fitYa.a
            star.axe = star.fitXa.aerr
            star.aye = star.fitYa.aerr
            star.cnt = star.velCnt
	    star.xchi2 = star.fitXa.chi2
	    star.ychi2 = star.fitYa.chi2
	    star.xchi2r = star.fitXa.chi2red
	    star.ychi2r = star.fitYa.chi2red
            star.t0x = star.fitXa.t0
            star.t0y = star.fitYa.t0
        except AttributeError, e:
            print '%s does not have acceleration fit' % star.name
            oldstars.remove(star)
    
    s.stars = oldstars

    if starlist=='all':
        tag = 'stars'
    if starlist=='oldstars':
        tag = 'old stars'
    if starlist=='yngstars':
        tag = 'young stars'

    print 'Found %d %s' % (len(oldstars),tag)

    return s


def nzErr(xerr, yerr, vxerr, vyerr, year_x, year_y, mag, alnDir = '13_08_21/', chainsDir = 'efit/chains_S0-2_newRV2/'):
    """
    Need to run this function to update the errors in X and Y of a given star. There is an
    error in position and velocity of sgr* from S0-2's efit. This uncertainty affects
    errors in positions of the stars, but not their velocity or acceleration errors in
    X and Y.

    xerr, yerr - Initial set of x and y coords in "

    year - epoch of these positions in years

    alnDir - align directory

    chainsDir - chains directory inside alnDir

    output - Final set of x and y errors in "
    """

    #Read in values for error in position and velocity of sgr*
    origin_val = asciidata.open('/g/ghez/align/' + alnDir + chainsDir + 'efit_summary.txt')
    ori_x0e = origin_val[18][0]
    ori_y0e = origin_val[19][0]
    ori_vxe = origin_val[20][0]
    ori_vye = origin_val[21][0]
    t_0 = 2000.0 #hard coded t_0 of sgr*

 #   magBins=np.array([9,11,12,13,14,15,16,17,18,19,20,21])
 #   deltaArr=np.array([3.5,71.0,58.0,210.0,300.0,650.0,700.0,1100.0,1900.0,2200.0,3000.0])*1e-6

#    delta = mag*0.0
#    for i in range(len(mag)):
#        for j in range(len(deltaArr)):
#            if ((mag[i] > magBins[j]) & (mag[i] <= magBins[j+1])):
#                delta[i]=deltaArr[j]

    #Update errors
    xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((year_x - t_0)*ori_vxe)**2)
    yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((year_y - t_0)*ori_vye)**2)
    vxerr = np.sqrt(vxerr**2 + ori_vxe**2)
    vyerr = np.sqrt(vyerr**2 + ori_vye**2)

    return xerr, yerr, vxerr, vyerr


def accelInfo(alnDir='14_06_18/', root_tmp='/g/ghez/align/',
              align = 'align/align_d_rms_1000_abs_t', updateErr = True,Rcut = 1.7,
              poly='polyfit_nz/fit', points='points_nz/', polyj='polyfit_nzj/fit',
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',
              starlist='all', magCut=22, nEpochs=14):

    """
    Gets radial acceleration and other info on stars that make mag, cnt, and 2D radius cuts
    """

    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)

    names = s.getArray('name')
#    print "stars that are in hist Accel:"
#    print names

    # In arcsec, mas/yr, mas/yr^2
    x0 = s.getArray('x0')
    y0 = s.getArray('y0')
    x0e = s.getArray('x0e')
    y0e = s.getArray('y0e')
    vx = s.getArray('vx')
    vy = s.getArray('vy')
    vxe = s.getArray('vxe')
    vye = s.getArray('vye')
    ax = s.getArray('ax')
    ay = s.getArray('ay')
    axe = s.getArray('axe')
    aye = s.getArray('aye')
    t0x = s.getArray('t0x')
    t0y = s.getArray('t0y')
    #xchi2r = s.getArray('xchi2r')
    #ychi2r = s.getArray('ychi2r')
    chi2xv = s.getArray('fitXv.chi2')
    chi2yv = s.getArray('fitYv.chi2')
    chi2xa = s.getArray('fitXa.chi2')
    chi2ya = s.getArray('fitYa.chi2')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')


    #Load Jerk values
    fitRoot = root_tmp + alnDir + polyj
    fitFile = fitRoot + '.accelFormal'
    t0File = fitRoot + '.t0'

    _fit = asciidata.open(fitFile)
    _t0 = asciidata.open(t0File)

    jt0x = _t0[1].tonumpy()
    jt0y = _t0[2].tonumpy()
    jnames = _fit[0].tonumpy()
    jx0 = _fit[1].tonumpy()
    jvx = _fit[2].tonumpy()
    jax = _fit[3].tonumpy()
    jx = _fit[4].tonumpy()
    jx0e = _fit[5].tonumpy()
    jvxe = _fit[6].tonumpy()
    jaxe = _fit[7].tonumpy()
    jxe = _fit[8].tonumpy()
    chi2xj = _fit[9].tonumpy()
#    jxq = _fit[10].tonumpy()

    jy0 = _fit[11].tonumpy()
    jvy = _fit[12].tonumpy()
    jay = _fit[13].tonumpy()
    jy = _fit[14].tonumpy()
    jy0e = _fit[15].tonumpy()
    jvye = _fit[16].tonumpy()
    jaye = _fit[17].tonumpy()
    jye = _fit[18].tonumpy()
    chi2yj = _fit[19].tonumpy()
#    jyq = _fit[20].tonumpy()

    r = np.sqrt(x0**2 + y0**2)
    idx = np.where((cnt > nEpochs) & (mag < magCut) & (r < Rcut))[0]
    x0 = x0[idx]
    y0 = y0[idx]
    x0e = x0e[idx]
    y0e = y0e[idx]
    vx = vx[idx]
    vy = vy[idx]
    vxe = vxe[idx]
    vye = vye[idx]
    ax = ax[idx]
    ay = ay[idx]
    axe = axe[idx]
    aye = aye[idx]
    cnt = cnt[idx]
    mag = mag[idx]
    #xchi2r = xchi2r[idx]
    #ychi2r = ychi2r[idx]
    t0x = t0x[idx]
    t0y = t0y[idx]
    chi2xv = chi2xv[idx]
    chi2yv = chi2yv[idx]
    chi2xa = chi2xa[idx]
    chi2ya = chi2ya[idx]
    chi2xj = chi2xj[idx]
    chi2yj = chi2yj[idx]
    names = [names[nn] for nn in idx]

    #Do F tests in this function, as opposed to making it own function
    #otherwise have to important data points twice
    #F tests from SYelda and accel_class, which calls on plotting functions

    signif = scipy.special.erfc(pvalue/math.sqrt(2.0))

    xFval = np.zeros(len(chi2xa))
    yFval = np.zeros(len(chi2xa))
    xFprob = np.zeros(len(chi2xa))
    yFprob = np.zeros(len(chi2xa))
    dof = np.zeros(len(chi2xa))

    good = np.where(chi2xa > 0)[0]
    dof[good] = cnt[good] - 3.0 #dof for acceleration

    xFval[good] = ((chi2xv[good] - chi2xa[good])/1.0)/(chi2xa[good]/dof[good])
    yFval[good] = ((chi2yv[good] - chi2ya[good])/1.0)/(chi2ya[good]/dof[good])
    
    xFprob[good] = stats.f.sf(xFval[good], 1, dof[good])
    yFprob[good] = stats.f.sf(yFval[good], 1, dof[good])

    # Round to nearest decimal
    xFprob = np.around(xFprob,decimals=4)
    yFprob = np.around(yFprob,decimals=4)

    pass_f = np.zeros(len(names))

    for ss in range(len(names)):
        xFp = xFprob[ss]
        yFp = yFprob[ss]

        if ((xFp < signif) | (yFp < signif)):
            pass_f[ss] = 1 # Accel fit
        else:
            pass_f[ss] = 0 # Velocity fit

    #Now F test between accel and jerk

    xFval = np.zeros(len(chi2xa))
    yFval = np.zeros(len(chi2xa))
    xFprob = np.zeros(len(chi2xa))
    yFprob = np.zeros(len(chi2xa))
    dof = np.zeros(len(chi2xa))

    good = np.where(chi2xj > 0)[0]
    dof[good] = cnt[good] - 4.0 #dof for jerk

    xFval[good] = ((chi2xa[good] - chi2xj[good])/1.0)/(chi2xj[good]/dof[good])
    yFval[good] = ((chi2ya[good] - chi2yj[good])/1.0)/(chi2yj[good]/dof[good])
    
    xFprob[good] = stats.f.sf(xFval[good], 1, dof[good])
    yFprob[good] = stats.f.sf(yFval[good], 1, dof[good])

    # Round to nearest decimal
    xFprob = np.around(xFprob,decimals=4)
    yFprob = np.around(yFprob,decimals=4)

    pass_fj = np.zeros(len(names))

    for ss in range(len(names)):
        xFp = xFprob[ss]
        yFp = yFprob[ss]

        if ((xFp < signif) | (yFp < signif)):
            pass_fj[ss] = 1 # Jerk fit
        else:
            pass_fj[ss] = 0 # Accel fit

    for i in range(len(names)):
        if (str(names[i]) == 'S0-16'):
            pass_f[i] == 1

    # S0-16 has such a high inclination that it is being seen as a velocity, force the f test vel/acc to be acc

    idex = np.where((pass_fj ==1) & (pass_f == 1))[0]

    for i in idex:
        jdex = np.where(jnames == names[i])
        x0[i] = jx0[jdex]
        y0[i] = jy0[jdex]
        x0e[i] = jx0e[jdex]
        y0e[i] = jy0e[jdex]
        vx[i] = jvx[jdex]
        vy[i] = jvy[jdex]
        vxe[i] = jvxe[jdex]
        vye[i] = jvye[jdex]
        ax[i] = jax[jdex]
        ay[i] = jay[jdex]
        axe[i] = jaxe[jdex]
        aye[i] = jaye[jdex]
        t0x[i] = jt0x[jdex]
        t0y[i] = jt0y[jdex]

    r = np.hypot(x0, y0)
    ar = ((ax*x0) + (ay*y0)) / r
    at = ((ax*y0) - (ay*x0)) / r
    are =  (axe*x0/r)**2 + (aye*y0/r)**2
    are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    are =  np.sqrt(are)
    ate =  (axe*y0/r)**2 + (aye*x0/r)**2
    ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
    ate =  np.sqrt(ate)

    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, alnDir=alnDir, chainsDir = chainsDir)

    if updateErr:
        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17.1])
        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5
        delta = mag*0.0
        for i in range(len(mag)):
            for j in range(len(atBins)-1):
                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
                    delta[i] = deltaArr[j]

        ate = np.sqrt(ate**2 + delta**2)
        are = np.sqrt(are**2 + delta**2)


    #Calculating z0 aand a_z and their errors
    x0_cm = x0 * cm_in_au * dist
    y0_cm = y0 * cm_in_au * dist
    x0e_cm = x0e * cm_in_au * dist
    y0e_cm = y0e * cm_in_au * dist
    r_cm = r * dist * cm_in_au
    re_cm = np.sqrt((x0_cm*x0e_cm/r_cm)**2 + (y0_cm*y0e_cm/r_cm)**2)
    ar_cmss = ar * asy_to_kms * 1e5 / sec_in_yr
    are_cmss = are * asy_to_kms * 1e5 / sec_in_yr
#    z0_cm = np.sqrt(abs(((GM * r_cm / ar_cmss)**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
#    r3d_cm = np.sqrt(x0_cm**2 + y0_cm**2 + z0_cm**2)
#    r3de_cm = np.sqrt(((x0_cm*x0e_cm)**2 + (y0_cm*y0e_cm)**2 + (z0_cm*z0e_cm)**2)/r3d_cm**2)

    return names, r_cm, ar_cmss, are_cmss





def gamma_C_accel(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                  poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
                  Rcut=1.7,nEpochs=14.,magCut=15.5,grange=[-5.0,1.9],maxBreak=5.0,alrange=[-10.0,10.0],
                  maxDelta=10.0,maxZ=1.0e20,globalt0=2006.0,maxR=2.0,flag='brokenPower',polyj='polyfit_nzj/fit',
                  pvalue=4.0):
    #Uses MultiNest in c to determine gamma
    #maxBreak in pc
    #maxZ in arcsec


    maxR *= dist * cm_in_au
    maxZ *= dist * cm_in_au

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    names, r2d, ar, are = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                    Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,
                                    chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs)

    #min_accel = -1.0*GM / r2d**2
    #max_accel = -1.0*GM*r2d / (np.sqrt(r2d**2 + maxZ**2))**3
    #maxr = np.sqrt(maxR**2 + maxZ**2)
    #maxz_star = np.sqrt(maxr**2 - r2d**2)

    dbfile = '/g/ghez/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
    oldProb = np.zeros(len(names))
    #Determine a star's prob of being old

    for i in range(len(names)):
        tmpName = str(names[i])
        try:
            cur.execute('SELECT name,young,old FROM stars WHERE name=?', [tmpName])
            for row in cur:
                if (row[1] == 'F') | (row[2] == 'T'):
                    oldProb[i] = 1.0
                    
        except:
            oldProb[i] = 0.0

    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
    for i in range(len(names)):
        tmpName = str(names[i])
        try:
            cur.execute('SELECT name,probOldSimPrior FROM unknownSims WHERE name=?', [tmpName])
            for row in cur:
                oldProb[i] = float(row[1])
        except:
            continue
        if ((tmpName == 'S0-38') | (tmpName == 'S0-49')):
            oldProb[i] = 1.0

    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/' + alnDir + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    #pdex = np.where(oldProb == 1.0)[0]
    num_stars = len(pdex)
    tmprint = 0
    print 'Stars used:'
    for p in pdex:
        print tmprint+5,names[p],oldProb[p]
        tmprint += 1


    #Save info to file for c to read
    savenames = [names[pp] for pp in pdex]
    r2d = r2d[pdex]
    ar = ar[pdex]
    are = are[pdex]
    oldProb = oldProb[pdex]
    np.savetxt('/u/schappell/code/c/stars_mn.dat',np.transpose([r2d,ar,are,oldProb]),delimiter=' ')

    os.system('g++ sc_mn.cpp')
