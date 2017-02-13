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
from matplotlib.ticker import ScalarFormatter 
#import pymultinest
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
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',nonRadial=0,
              starlist='all', magCut=22, lowMag=-10.0, nEpochs=14,ak_correct=True):

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

    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

    if (ak_correct == True):
        for mm in range(len(names)):
            try:
                cur.execute('SELECT Ak_sch FROM stars WHERE name=?', [str(names[mm])])
                for row in cur:
                    mag[mm] += 2.7 - row[0]
            except:
                pdb.set_trace()

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
    if (nonRadial==0):
        idx = np.where((cnt > nEpochs) & (mag < magCut) & (mag > lowMag) & (r < Rcut))[0]
    else:
        in_gcows = np.zeros(len(mag))
        plate_scale = 0.00995
        gcows_zero = np.array([1500.0,1500.0])
        gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')
        gcows_dex = np.where(gcows > 0)
        np.savetxt('/u/schappell/code/c/gcows_field.dat',np.transpose([gcows_dex[1],gcows_dex[0]]),delimiter=' ')
        #saves x and y positions, in pixels
        for ii in range(len(mag)):
            #pdb.set_trace()
            try:
                x_pixels = int(round((-1.0*x0[ii] / plate_scale) + gcows_zero[0]))
                y_pixels = int(round((y0[ii] / plate_scale) + gcows_zero[1]))
                in_gcows[ii] = gcows[y_pixels,x_pixels]
            except:
                continue
            if ((x_pixels < 0.0) | (y_pixels < 0.0)):
                in_gcows[ii] = 0
        idx = np.where((cnt > nEpochs) & (mag < magCut) & (mag > lowMag) & (in_gcows==1) & (r < Rcut))[0]
        #pdb.set_trace()
        
    #for ii in range(len(names)):
     #   if (names[ii] == 'S0-17'):
     #       pdb.set_trace()
     #       idx = np.append(idx,ii)
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
    #chi2xj = chi2xj[idx]
    #chi2yj = chi2yj[idx]
    names = [names[nn] for nn in idx]

    chi2xj_tmp = np.array([])
    chi2yj_tmp = np.array([])
    for tname in names:
        tdex = np.where(jnames == tname)[0]
        chi2xj_tmp = np.append(chi2xj_tmp,chi2xj[tdex])
        chi2yj_tmp = np.append(chi2yj_tmp,chi2yj[tdex])


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

    good = np.where(chi2xj_tmp > 0)[0]
    dof[good] = cnt[good] - 4.0 #dof for jerk

    xFval[good] = ((chi2xa[good] - chi2xj_tmp[good])/1.0)/(chi2xj_tmp[good]/dof[good])
    yFval[good] = ((chi2ya[good] - chi2yj_tmp[good])/1.0)/(chi2yj_tmp[good]/dof[good])
    
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
        if ((str(names[i]) == 'S0-16') | (str(names[i]) == 'S0-17')):
            pass_f[i] = 1
            pass_fj[i] = 1
            print 'S0-16 or S0-17 updated accel and jerk'

    # S0-16 and S0-17 have such high inclinations that they are being seen as velocity, force the f test vel/acc to be acc

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

#    pdb.set_trace()

    return names, r_cm, ar_cmss, are_cmss




def PMNGamma_NOPE(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                  poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
                  Rcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxZ=100.0):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    pdb.set_trace()
    names = [names[nn] for nn in cut]


    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
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
        try:
            cur.execute('SELECT name,probOldSimPrior FROM unknownSims WHERE name=?', [tmpName])
            for row in cur:
                oldProb[i] = row[1]
        except:
            continue
        if (tmpName == 'S0-32'):
            oldProb[i] = 0.0


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #define prior
    #single power law
    #uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        cube[1] = (cube[1] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
        cube[2] = (cube[2] * (2.0 * maxR) - maxR) * dist * cm_in_au
        cube[3] = cube[3] * maxZ * dist * cm_in_au #z0

    #define model
    #position part, orbit fit up to jerk
    def pos_lnL(gmod,xmod,ymod,zmod,xdata,ydata,tdata,xe_data,ye_data,FinvX,FinvY):
        amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
        xtilda = (xdata - xmod - xmod*amod*(tdata - globalt0)**2.0)
        ytilda = (ydata - ymod - ymod*amod*(tdata - globalt0)**2.0)

        axarray = np.zeros(2)
        ayarray = np.zeros(2)
        axarray[0] = np.sum(xtilda*(tdata - globalt0)/xe_data**2.0)
        axarray[1] = np.sum(xtilda*(tdata - globalt0)**3.0/xe_data**2.0)
        ayarray[0] = np.sum(ytilda*(tdata - globalt0)/ye_data**2.0)
        ayarray[1] = np.sum(ytilda*(tdata - globalt0)**3.0/ye_data**2.0)

        x_af = np.dot(axarray,FinvX)
        y_af = np.dot(ayarray,FinvY)

        return (np.sum(xtilda**2.0/xe_data**2.0) - np.dot(x_af,axarray) +
                np.sum(ytilda**2.0/ye_data**2.0) - np.dot(y_af,ayarray))/-2.0

    #rho part, integrating over 3D cube in space
    def rho_lnL(gmod,xmod,ymod,zmod):
        norm=integrate.tplquad(lambda zprime,yprime,xprime:(xprime**2.0 + yprime**2.0 + zprime**2.0)**(gmod/-2.0),
                               -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, 
                               lambda xprime, yprime: 0.0, lambda xprime, yprime: maxZ)
        norm = norm[0]
        return math.log((xmod**2.0 + ymod**2.0 + zmod**2.0)**(gmod/-2.0) / norm)

    #entire ln L

    #place holders for now to define function
    gmod,xmod,ymod,zmod = None, None, None, None
    xp,yp,time,xerr,yerr,inv_xf,inv_yf = None, None, None, None, None, None, None

    def SPLloglike(cube,ndim,nparams):
        gmod,xmod,ymod,zmod = cube[0],cube[1],cube[2],cube[3]
        lnL_1 = pos_lnL(gmod,xmod,ymod,zmod,xp,yp,time,xerr,yerr,inv_xf,inv_yf)
        lnL_2 = rho_lnL(gmod,xmod,ymod,zmod)

        return lnL_1 + lnL_2


    parameters = ['Gamma','x0','y0','z0']
    n_params = len(parameters)

    print ''
    print ''
    print ''
    print 'Stars and prob of being old:'

    for n in range(len(names)):
        if (oldProb[n] > 0.0):
            tmpName = names[n]
            print tmpName,oldProb[n]

            pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

            time = pointsTab[0].tonumpy()
            xp = pointsTab[1].tonumpy() * dist * cm_in_au
            yp = pointsTab[2].tonumpy() * dist * cm_in_au
            xerr = pointsTab[3].tonumpy()
            yerr = pointsTab[4].tonumpy()

            if (updateErr == True):
                xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
                yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


            xerr = xerr * dist * cm_in_au
            yerr = yerr * dist * cm_in_au
            time = time * sec_in_yr
            t_0 = t_0 * sec_in_yr
            globalt0 = globalt0 * sec_in_yr

            p_error2 = 1.0e300
            xfarray = np.zeros([2,2])
            xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
            xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
            xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
            xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

            inv_xf = np.linalg.inv(xfarray)

            yfarray = np.zeros([2,2])
            yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
            yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
            yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
            yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2

            inv_yf = np.linalg.inv(yfarray)


            #run MultiNest for this star
            pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+tmpName+'_pmn',
                            resume = False,verbose = True)
            
            pdb.set_trace()


def testPMN(alnDir='14_06_18/', root_tmp='/g/ghez/align/',align = 'align/align_d_rms_1000_abs_t', 
            poly='polyfit_nz/fit', points='points_nz/', polyj='polyfit_nzj/fit',star='S1-13',
            chainsDir = 'efit/chains_S0-2_newRV2/', starlist='all', verbose = True):

    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)

    names = s.getArray('name')
#    print "stars that are in hist Accel:"
#    print names

    # In arcsec, mas/yr, mas/yr^2
    x0 = s.getArray('x0')
    y0 = s.getArray('y0')
    x0e = s.getArray('x0e')
    y0e = s.getArray('y0e')
    #vx = s.getArray('vx')
    #vy = s.getArray('vy')
    #vxe = s.getArray('vxe')
    #vye = s.getArray('vye')
    ax = s.getArray('ax')
    ay = s.getArray('ay')
    axe = s.getArray('axe')
    aye = s.getArray('aye')
    t0x = s.getArray('t0x')
    t0y = s.getArray('t0y')
    xchi2r = s.getArray('xchi2r')
    ychi2r = s.getArray('ychi2r')
        
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

    r = np.hypot(x0, y0)
    ar = ((ax*x0) + (ay*y0)) / r
    at = ((ax*y0) - (ay*x0)) / r
    are =  (axe*x0/r)**2 + (aye*y0/r)**2
    are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    are =  np.sqrt(are)
    #ate =  (axe*y0/r)**2 + (aye*x0/r)**2
    #ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
    #ate =  np.sqrt(ate)

    def prior(cube,ndim,nparams):
        #cube[0] = cube[0] * (1.9 + 2.0) - 2.0 #gamma
        cube[0] = cube[0] * 1.0 * dist * cm_in_au #z0

    def loglike(cube,ndim,nparams):
        gmod,zmod = cube[0],cube[1]
        a_p = -1.0*GM*rtmp/(math.sqrt(rtmp**2.0 + zmod**2.0))**3.0
        prob_ap = math.exp(-1.0*(atmp - a_p)**2.0/(2.0*aetmp**2.0))
        #norm1 = abs(aetmp*math.sqrt(pi/2.0)*(special.erf((atmp - np.min(a_p))/(math.sqrt(2.0)*aetmp)) - 
        #                                     special.erf((atmp - np.max(a_p))/(math.sqrt(2.0)*aetmp))))
        prob_a = prob_ap
        #rho_tmp = (math.sqrt(rtmp**2.0 + zmod**2.0))**(-1.0*gmod)
        #Rcut = 3.0 * dist * cm_in_au
        #rcut_tmp = math.sqrt(Rcut**2.0 + (100.0*dist*cm_in_au)**2.0)
        #norm2a = integrate.quad(lambda x: x**(2.0-gmod), 0.0, rcut_tmp)
        #norm2b = integrate.quad(lambda x: x**2.0 * (rtmp**2.0 + x**2)**(gmod/-2.0), 0.0, 
        #                        math.sqrt(rcut_tmp**2.0 - Rcut**2.0))
        #norm2 = norm2a[0] - norm2b[0]
        #prob_z = rho_tmp/norm2

        return math.log(prob_a)

    parameters = ['z0']
    n_params = len(parameters)

    for i in range(len(names)):
        tmpName = names[i]
        if (tmpName == str(star)):
            rtmp = r[i] * dist * cm_in_au
            atmp = ar[i] * asy_to_kms * 1e5 / sec_in_yr
            aetmp = are[i] * asy_to_kms * 1e5 / sec_in_yr

            #pdb.set_trace()

            pymultinest.run(loglike,prior,n_params,outputfiles_basename='pmnOld/test_pmn',
                            resume = False,verbose = True)
            
            pdb.set_trace()




def grrPMN():

    #let's try fitting a line with PMN

    xdata = np.array([0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0])
    ydata = np.array([0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0])

    def prior(cube,ndim,nparams):
        cube[0] = cube[0] * 10.0 - 5.0
        cube[1] = cube[1] * 10.0 - 5.0

    def loglike(cube,ndim,nparams):
        bint,slope = cube[0],cube[1]
        like = np.sum((ydata - (bint + xdata*slope))**2.0/-2.0)

        return like

    parameters = ['intercept','slope']
    n_params = len(parameters)

    pymultinest.run(loglike,prior,n_params,outputfiles_basename='pmnOld/slope',
                    resume = False,verbose = True)






def PMNGamma_xyz(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                 poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
                 Rcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxZ=10.0,
                 flag='fullTEST'):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    #pdb.set_trace()
    names = [names[nn] for nn in cut]
    xprint = xfcut[cut]
    yprint = yfcut[cut]


    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
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
        try:
            cur.execute('SELECT name,probOldSimPrior FROM unknownSims WHERE name=?', [tmpName])
            for row in cur:
                oldProb[i] = row[1]
        except:
            continue
        if (tmpName == 'S0-32'):
            oldProb[i] = 0.0


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*
        pdb.set_trace()

    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    num_stars = len(pdex)

    for p in pdex:
        print names[p],xprint[p]*dist*cm_in_au,yprint[p]*dist*cm_in_au

    pdb.set_trace()


    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        for i in range(len(pdex)):
            cube[i*3+1] = (cube[i*3+1] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            cube[i*3+2] = (cube[i*3+2] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            cube[i*3+3] = cube[i*3+3] * maxZ * dist * cm_in_au #z0
            

    def SPLloglike(cube,ndim,nparams):
        gmod = cube[0]
        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            tmpName = names[p]
            lnL_1 = pos_lnL(tmpName,cube[i*3+1],cube[i*3+2],cube[i*3+3])
            lnL_2 = rho_lnL(gmod,cube[i*3+1],cube[i*3+2],cube[i*3+3])
            total_lnL += oldProb[p]*(lnL_1 + lnL_2)

        return total_lnL

    n_params = 3*len(pdex)+1


    #position part, orbit fit up to jerk
    def pos_lnL(starName,xmod,ymod,zmod,t_0=t_0,globalt0=globalt0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)

        amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
        xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
        ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

        axarray = np.zeros(2)
        ayarray = np.zeros(2)
        axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
        axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
        ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
        ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

        x_af = np.dot(axarray,inv_xf)
        y_af = np.dot(ayarray,inv_yf)

        return (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0

    #rho part, integrating over 3D cube in space
    def rho_lnL(gmod,xmod,ymod,zmod):
        norm=integrate.tplquad(lambda zprime,yprime,xprime:(xprime**2.0 + yprime**2.0 + zprime**2.0)**(gmod/-2.0),
                               -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, 
                               lambda xprime, yprime: 0.0, lambda xprime, yprime: maxZ)
        norm = norm[0]
        return math.log((xmod**2.0 + ymod**2.0 + zmod**2.0)**(gmod/-2.0) / norm)


    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')




def PMNGamma_singlePower(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                         poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
                         Rcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxZ=10.0,
                         flag='fullTEST'):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    #pdb.set_trace()
    names = [names[nn] for nn in cut]
    #xprint = xfcut[cut]
    #yprint = yfcut[cut]

    maxZ *= dist * cm_in_au
    maxR *= dist * cm_in_au


    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
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
        try:
            cur.execute('SELECT name,probOldSimPrior FROM unknownSims WHERE name=?', [tmpName])
            for row in cur:
                oldProb[i] = row[1]
        except:
            continue
        if (tmpName == 'S0-32'):
            oldProb[i] = 0.0



    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    num_stars = len(pdex)
    print 'Stars used:'
    for p in pdex:
        print names[p],oldProb[p]


    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        for i in range(len(pdex)):
            #cube[i*3+1] = (cube[i*3+1] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            #cube[i*3+2] = (cube[i*3+2] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            cube[i+1] = cube[i+1] * maxZ  #z0
            

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        gmod = cube[0]
        #norm=integrate.tplquad(lambda zprime,yprime,xprime:(xprime**2.0 + yprime**2.0 + zprime**2.0)**(gmod/-2.0),
        #                       -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, 
        #                       lambda xprime, yprime: 0.0, lambda xprime, yprime: maxZ)
        norm=integrate.dblquad(lambda xprime,yprime: maxZ*(maxZ**2.0+xprime**2.0+yprime**2.0)**(gmod/-2.0)*
                               (maxZ**2.0/(yprime**2.0+xprime**2.0)+1.0)**(gmod/2.0)*
                               scipy.special.hyp2f1(0.5,gmod/2.0,1.5,-1.0*maxZ**2.0/(xprime**2.0+yprime**2.0)),
                               0.0,maxR, lambda xprime: 0.0, lambda xprime: maxR)
        norm = 4.0*norm[0]
        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            tmpName = names[p]
            tmplnL = star_lnL(tmpName,gmod,cube[i+1],norm)
            total_lnL += tmplnL*oldProb[p]

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL


    def star_lnL(starName,gmod,zmod,norm,t_0=t_0,globalt0=globalt0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)


        def like_toSum(xmod,ymod,zmod=zmod,gmod=gmod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0,inv_xf=inv_xf,inv_yf=inv_yf,norm=norm):
            amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
            xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
            ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

            axarray = np.zeros(2)
            ayarray = np.zeros(2)
            axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
            axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
            ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
            ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

            x_af = np.dot(axarray,inv_xf)
            y_af = np.dot(ayarray,inv_yf)

            lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                       np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0
            return math.exp(lnL_pos)*(xmod**2 + ymod**2 + zmod**2)**(gmod/-2.0) / norm

        def int_y(maxR = maxR):
            return [-1.0*maxR, maxR]

        def int_x(ymod,maxR=maxR):
            return [-1.0*maxR, maxR]
        
        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        if (starlike == 0.0):
            starlike = 1e-323
        return math.log(starlike)


    n_params = len(pdex)+1

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def PMNGamma(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
	     poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
	     Rcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxZ=10.0,
             maxBreak=10.0,alrange=[-10.0,10.0],maxDelta=10.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    names = [names[nn] for nn in cut]
    #xprint = xfcut[cut]
    #yprint = yfcut[cut]
    #pdb.set_trace()
    maxZ *= dist * cm_in_au
    maxR *= dist * cm_in_au


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
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    num_stars = len(pdex)
    tmprint = 0
    print 'Stars used:'
    for p in pdex:
        print tmprint+5,names[p],oldProb[p]
        tmprint += 1


    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        cube[1] = cube[1] * (alrange[1] - alrange[0]) + alrange[0] #alpha
        cube[2] = cube[2] * maxDelta #delta
        cube[3] = cube[3] * maxBreak * dist * cm_in_au #r break
        for i in range(len(pdex)):
            #cube[i*3+1] = (cube[i*3+1] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            #cube[i*3+2] = (cube[i*3+2] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            cube[i+4] = cube[i+4] * maxZ  #z0
            

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        gmod = cube[0]
        almod = cube[1]
        demod = cube[2]
        brmod = cube[3]
        #norm=integrate.tplquad(lambda zprime,yprime,xprime:(xprime**2.0 + yprime**2.0 + zprime**2.0)**(gmod/-2.0),
        #                       -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, 
        #                       lambda xprime, yprime: 0.0, lambda xprime, yprime: maxZ)
        #norm=integrate.dblquad(lambda xprime,yprime: maxZ*(maxZ**2.0+xprime**2.0+yprime**2.0)**(gmod/-2.0)*
        #                       (maxZ**2.0/(yprime**2.0+xprime**2.0)+1.0)**(gmod/2.0)*
        #                       scipy.special.hyp2f1(0.5,gmod/2.0,1.5,-1.0*maxZ**2.0/(xprime**2.0+yprime**2.0)),
        #                       0.0,maxR, lambda xprime: 0.0, lambda xprime: maxR)
        norm = integrate.dblquad(lambda Rprime,zprime: 2.0*pi*Rprime*(Rprime**2+zprime**2)**(gmod/-2.0)*
                                 (1.0+(math.sqrt(Rprime**2+zprime**2)/brmod)**demod)**((gmod-almod)/demod),
                                 0.0,maxR, lambda Rprime: 0.0, lambda Rprime: maxZ)
        norm = norm[0]
        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            tmpName = names[p]
            tmplnL = star_lnL(tmpName,gmod,almod,demod,brmod,cube[i+4],norm)
            total_lnL += tmplnL*oldProb[p]

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL


    def star_lnL(starName,gmod,almod,demod,brmod,zmod,norm,t_0=t_0,globalt0=globalt0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)


        def like_toSum(xmod,ymod,zmod=zmod,gmod=gmod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0,inv_xf=inv_xf,inv_yf=inv_yf,norm=norm):
            amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
            xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
            ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

            axarray = np.zeros(2)
            ayarray = np.zeros(2)
            axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
            axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
            ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
            ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

            x_af = np.dot(axarray,inv_xf)
            y_af = np.dot(ayarray,inv_yf)

            lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                       np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0
            like_den = ((xmod**2+ymod**2+zmod**2)**(gmod/-2.0)*
                        (1.0+(math.sqrt(xmod**2+ymod**2+zmod**2)/brmod)**demod)**((gmod-almod)/demod) / norm)
            return math.exp(lnL_pos)*like_den

        def int_y(maxR = maxR):
            return [-1.0*maxR, maxR]

        def int_x(ymod,maxR=maxR):
            xlim = math.sqrt(maxR**2-ymod**2)
            return [-1.0*xlim, xlim]
        
        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        if (starlike == 0.0):
            starlike = 1e-323
        return math.log(starlike)


    n_params = len(pdex)+4

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')



def PMNGamma_noJerk(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                    poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
                    Rcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxZ=10.0,
                    maxBreak=10.0,alrange=[-10.0,10.0],maxDelta=10.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    #pdb.set_trace()
    names = [names[nn] for nn in cut]
    #xprint = xfcut[cut]
    #yprint = yfcut[cut]

    maxZ *= dist * cm_in_au
    maxR *= dist * cm_in_au

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
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    num_stars = len(pdex)
    tmprint = 0
    print 'Stars used:'
    for p in pdex:
        print tmprint+5,names[p],oldProb[p]
        tmprint += 1


    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        cube[1] = cube[1] * (alrange[1] - alrange[0]) + alrange[0] #alpha
        cube[2] = cube[2] * maxDelta #delta
        cube[3] = cube[3] * maxBreak * dist * cm_in_au #r break
        for i in range(len(pdex)):
            #cube[i*3+1] = (cube[i*3+1] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            #cube[i*3+2] = (cube[i*3+2] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            cube[i+4] = cube[i+4] * maxZ  #z0
            

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        gmod = cube[0]
        almod = cube[1]
        demod = cube[2]
        brmod = cube[3]
        #norm=integrate.tplquad(lambda zprime,yprime,xprime:(xprime**2.0 + yprime**2.0 + zprime**2.0)**(gmod/-2.0),
        #                       -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, 
        #                       lambda xprime, yprime: 0.0, lambda xprime, yprime: maxZ)
        #norm=integrate.dblquad(lambda xprime,yprime: maxZ*(maxZ**2.0+xprime**2.0+yprime**2.0)**(gmod/-2.0)*
        #                       (maxZ**2.0/(yprime**2.0+xprime**2.0)+1.0)**(gmod/2.0)*
        #                       scipy.special.hyp2f1(0.5,gmod/2.0,1.5,-1.0*maxZ**2.0/(xprime**2.0+yprime**2.0)),
        #                       0.0,maxR, lambda xprime: 0.0, lambda xprime: maxR)
        norm = integrate.dblquad(lambda Rprime,zprime: 2.0*pi*Rprime*(Rprime**2+zprime**2)**(gmod/-2.0)*
                                 (1.0+(math.sqrt(Rprime**2+zprime**2)/brmod)**demod)**((gmod-almod)/demod),
                                 0.0,maxR, lambda Rprime: 0.0, lambda Rprime: maxZ)
        norm = norm[0]
        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            tmpName = names[p]
            tmplnL = star_lnL(tmpName,gmod,almod,demod,brmod,cube[i+4],norm)
            total_lnL += tmplnL*oldProb[p]

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL


    def star_lnL(starName,gmod,almod,demod,brmod,zmod,norm,t_0=t_0,globalt0=globalt0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        #p_error2 = 1.0e300
        #xfarray = np.zeros([2,2])
        #xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        #xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        #xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        #xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        #inv_xf = np.linalg.inv(xfarray)
        
        #yfarray = np.zeros([2,2])
        #yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        #yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        #yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        #yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        #inv_yf = np.linalg.inv(yfarray)


        def like_toSum(xmod,ymod,zmod=zmod,gmod=gmod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0,norm=norm):
            amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
            xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
            ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

            #axarray = np.zeros(2)
            #ayarray = np.zeros(2)
            #axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
            #axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
            #ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
            #ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

            #x_af = np.dot(axarray,inv_xf)
            #y_af = np.dot(ayarray,inv_yf)

            #lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
            #           np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0

            lnL_pos = (-1.0*np.sum(xtilda**2/(2.0*xerr**2))+(np.sum(xtilda*(time-globalt0)/(2.0*xerr**2)))**2/
                        np.sum((time-globalt0)**2/(2.0*xerr**2))-
                        1.0*np.sum(ytilda**2/(2.0*yerr**2))+(np.sum(ytilda*(time-globalt0)/(2.0*yerr**2)))**2/
                        np.sum((time-globalt0)**2/(2.0*yerr**2)))
            like_den = ((xmod**2+ymod**2+zmod**2)**(gmod/-2.0)*
                        (1.0+(math.sqrt(xmod**2+ymod**2+zmod**2)/brmod)**demod)**((gmod-almod)/demod) / norm)
            return math.exp(lnL_pos)*like_den

        def int_y(maxR = maxR):
            return [-1.0*maxR, maxR]

        def int_x(ymod,maxR=maxR):
            xlim = math.sqrt(maxR**2-ymod**2)
            return [-1.0*xlim, xlim]
        
        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        if (starlike == 0.0):
            starlike = 1e-323
        return math.log(starlike)


    n_params = len(pdex)+4

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')




def PMNGamma_noXY(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                  poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
                  Rcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxZ=10.0,
                  maxBreak=10.0,alrange=[-10.0,10.0],maxDelta=10.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    names = [names[nn] for nn in cut]
    #xprint = xfcut[cut]
    #yprint = yfcut[cut]
    #pdb.set_trace()
    maxZ *= dist * cm_in_au
    maxR *= dist * cm_in_au


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
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    num_stars = len(pdex)
    tmprint = 0
    print 'Stars used:'
    for p in pdex:
        print tmprint+5,names[p],oldProb[p]
        tmprint += 1


    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        cube[1] = cube[1] * (alrange[1] - alrange[0]) + alrange[0] #alpha
        cube[2] = cube[2] * maxDelta #delta
        cube[3] = cube[3] * maxBreak * dist * cm_in_au #r break
        for i in range(len(pdex)):
            #cube[i*3+1] = (cube[i*3+1] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            #cube[i*3+2] = (cube[i*3+2] * (2.0 * maxR) - maxR) * dist * cm_in_au #x0 and y0
            cube[i+4] = cube[i+4] * maxZ  #z0
            

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        gmod = cube[0]
        almod = cube[1]
        demod = cube[2]
        brmod = cube[3]
        #norm=integrate.tplquad(lambda zprime,yprime,xprime:(xprime**2.0 + yprime**2.0 + zprime**2.0)**(gmod/-2.0),
        #                       -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, 
        #                       lambda xprime, yprime: 0.0, lambda xprime, yprime: maxZ)
        #norm=integrate.dblquad(lambda xprime,yprime: maxZ*(maxZ**2.0+xprime**2.0+yprime**2.0)**(gmod/-2.0)*
        #                       (maxZ**2.0/(yprime**2.0+xprime**2.0)+1.0)**(gmod/2.0)*
        #                       scipy.special.hyp2f1(0.5,gmod/2.0,1.5,-1.0*maxZ**2.0/(xprime**2.0+yprime**2.0)),
        #                       0.0,maxR, lambda xprime: 0.0, lambda xprime: maxR)
        norm = integrate.dblquad(lambda Rprime,zprime: 2.0*pi*Rprime*(Rprime**2+zprime**2)**(gmod/-2.0)*
                                 (1.0+(math.sqrt(Rprime**2+zprime**2)/brmod)**demod)**((gmod-almod)/demod),
                                 0.0,maxR, lambda Rprime: 0.0, lambda Rprime: maxZ)
        norm = norm[0]
        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            tmpName = names[p]
            tmplnL = star_lnL(tmpName,gmod,almod,demod,brmod,cube[i+4],norm)
            total_lnL += tmplnL*oldProb[p]

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL


    def star_lnL(starName,gmod,almod,demod,brmod,zmod,norm,t_0=t_0,globalt0=globalt0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)


        def like_toSum(zmod=zmod,gmod=gmod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0,inv_xf=inv_xf,inv_yf=inv_yf,norm=norm):

            findmin = np.argmax(xp**2+yp**2)
            xmod = xp[findmin]
            ymod = yp[findmin]
            amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
            xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
            ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

            axarray = np.zeros(2)
            ayarray = np.zeros(2)
            axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
            axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
            ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
            ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

            x_af = np.dot(axarray,inv_xf)
            y_af = np.dot(ayarray,inv_yf)

            lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                       np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0
            like_den = ((xmod**2+ymod**2+zmod**2)**(gmod/-2.0)*
                        (1.0+(math.sqrt(xmod**2+ymod**2+zmod**2)/brmod)**demod)**((gmod-almod)/demod) / norm)
            return math.exp(lnL_pos)*like_den

        #def int_y(maxR = maxR):
         #   return [-1.0*maxR, maxR]

        #def int_x(ymod,maxR=maxR):
         #   xlim = math.sqrt(maxR**2-ymod**2)
          #  return [-1.0*xlim, xlim]
        
        #starlike = integrate.nquad(like_toSum,[int_x,int_y])
        #starlike = starlike[0]
        starlike = like_toSum()
        if (starlike == 0.0):
            starlike = 1e-323
        return math.log(starlike)


    n_params = len(pdex)+4

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def PMN_z(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
          poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
          globalt0 = 2013.318,maxR=2.0,maxA=1.2,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    maxR *= dist * cm_in_au


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * maxA
            

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        amod = cube[0]

        total_lnL = star_lnL(star,amod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,amod,t_0=t_0,globalt0=globalt0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)


        def like_toSum(xmod,ymod,amod=amod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0,inv_xf=inv_xf,inv_yf=inv_yf):
            #amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
            xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
            ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

            axarray = np.zeros(2)
            ayarray = np.zeros(2)
            axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
            axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
            ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
            ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

            x_af = np.dot(axarray,inv_xf)
            y_af = np.dot(ayarray,inv_yf)

            lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                       np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0
            return math.exp(lnL_pos)

        def int_y(maxR = maxR):
            return [-1.0*maxR, maxR]

        def int_x(ymod,maxR=maxR):
            xlim = math.sqrt(maxR**2-ymod**2)
            return [-1.0*xlim, xlim]
        
        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        if (starlike == 0.0):
            starlike = 1e-323
        return math.log(starlike)


    n_params = 1

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def PMN_z_noJerk(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                 poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
                 globalt0 = 2013.318,grange=[-2.0,1.9],maxR=2.0,maxA=1.2,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma

    #maxZ *= dist * cm_in_au
    maxR *= dist * cm_in_au


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*


    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * maxA  #z0
            

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        amod = cube[0]
        
        total_lnL = star_lnL(star,amod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL


    tmpName = star
    def star_lnL(starName,amod,t_0=t_0,globalt0=globalt0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        #p_error2 = 1.0e300
        #xfarray = np.zeros([2,2])
        #xfarray[0,0] = np.sum((time - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        #xfarray[0,1] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        #xfarray[1,0] = np.sum((time - globalt0)**4.0/(xerr**2.0))
        #xfarray[1,1] = np.sum((time - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        #inv_xf = np.linalg.inv(xfarray)
        
        #yfarray = np.zeros([2,2])
        #yfarray[0,0] = np.sum((time - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        #yfarray[0,1] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        #yfarray[1,0] = np.sum((time - globalt0)**4.0/(yerr**2.0))
        #yfarray[1,1] = np.sum((time - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        #inv_yf = np.linalg.inv(yfarray)


        def like_toSum(xmod,ymod,amod=amod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0):
            #amod = -0.5*GM / (xmod**2.0 + ymod**2.0 + zmod**2.0)**(3.0/2.0)
            xtilda = (xp - xmod - xmod*amod*(time - globalt0)**2.0)
            ytilda = (yp - ymod - ymod*amod*(time - globalt0)**2.0)

            #axarray = np.zeros(2)
            #ayarray = np.zeros(2)
            #axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
            #axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
            #ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
            #ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

            #x_af = np.dot(axarray,inv_xf)
            #y_af = np.dot(ayarray,inv_yf)

            #lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
            #           np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0

            lnL_pos = (-1.0*np.sum(xtilda**2/(2.0*xerr**2))+(np.sum(xtilda*(time-globalt0)/(2.0*xerr**2)))**2/
                        np.sum((time-globalt0)**2/(2.0*xerr**2))-
                        1.0*np.sum(ytilda**2/(2.0*yerr**2))+(np.sum(ytilda*(time-globalt0)/(2.0*yerr**2)))**2/
                        np.sum((time-globalt0)**2/(2.0*yerr**2)))

            return math.exp(lnL_pos)

        def int_y(maxR = maxR):
            return [-1.0*maxR, maxR]

        def int_x(ymod,maxR=maxR):
            xlim = math.sqrt(maxR**2-ymod**2)
            return [-1.0*xlim, xlim]
        
        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        if (starlike == 0.0):
            starlike = 1e-323
        return math.log(starlike)


    n_params = 1

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def testAccel(r2d,ar,are,maxZ=0.05):
    r2d = r2d * dist * cm_in_au
    maxZ = maxZ * cm_in_pc

    ar = ar * asy_to_kms * 0.1 / sec_in_yr
    are = are * asy_to_kms * 0.1 / sec_in_yr

    #zzs = np.array([i*maxZ/100.0 for i in range(101)])
    #atmp = -1.0*GM*r2d / (r2d**2 + zzs**2)**(3.0/2.0)
    #toplot = -0.5*(ar - atmp)**2 / are**2.0

    #pdb.set_trace()

    def testprior(cube,ndim,nparams):
        cube[0] = cube[0] * maxZ  #z0


    def testloglike(cube,ndim,nparams):
        zmod = cube[0]

        amod = -1.0*GM*r2d / (r2d**2 + zmod**2)**(3.0/2.0)

        lnL = -0.5*(ar - amod)**2 / are**2.0
        return lnL

    n_params = 1
    pymultinest.run(testloglike,testprior,n_params,outputfiles_basename='pmnOld/testz_accel_')




def PMN_xva(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
            poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
            globalt0 = 2013.318,maxR=2.0,maxZ=100.0,maxV=1.0e4,maxA=1.0e4,maxJ=0.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxR in "
    #maxZ in mpc
    #maxV in km/s
    #maxA in micro-as/yr^2
    #maxA in micro-as/yr^3

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_pc / 1000.
    #maxV *= 1.0e5
    #maxA *= asy_to_kms * 0.1 / sec_in_yr


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * 2.0 * maxR - maxR #x0 in cm
        cube[1] = cube[1] * 2.0 * maxR - maxR #y0 in cm
        #cube[2] = cube[2] * maxZ #z0 in cm
        cube[2] = cube[2] * 2.0 * maxV - maxV #vx in cm/s
        cube[3] = cube[3] * 2.0 * maxV - maxV #vy in cm/s
        cube[4] = cube[4] * 2.0 * maxA - maxA #ax in cm/s^2
        cube[5] = cube[5] * 2.0 * maxA - maxA #ay in cm/s^2
        if (maxJ > 0.0):
            cube[6] = cube[6] * 2.0 * maxJ - maxJ #x and y jerk terms
            cube[7] = cube[7] * 2.0 * maxJ - maxJ
        

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        x0mod = cube[0] * dist * cm_in_au
        y0mod = cube[1] * dist * cm_in_au
        vxmod = cube[2] * 1.0e5
        vymod = cube[3] * 1.0e5
        axmod = cube[4] * asy_to_kms * 0.1 / sec_in_yr
        aymod = cube[5] * asy_to_kms * 0.1 / sec_in_yr
        if (maxJ > 0.0):
            jxmod = cube[6] * asy_to_kms * 0.1 / sec_in_yr**2
            jymod = cube[7] * asy_to_kms * 0.1 / sec_in_yr**2
            
            total_lnL = star_lnL(star,x0mod,y0mod,vxmod,vymod,axmod,aymod,jxmod=jxmod,jymod=jymod)
        else:
            total_lnL = star_lnL(star,x0mod,y0mod,vxmod,vymod,axmod,aymod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,x0mod,y0mod,vxmod,vymod,axmod,aymod,t_0=t_0,globalt0=globalt0,maxR=maxR,jxmod=0.0,jymod=0.0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        xmod = x0mod + vxmod*(time - globalt0) + 0.5 * axmod * (time - globalt0)**2
        ymod = y0mod + vymod*(time - globalt0) + 0.5 * aymod * (time - globalt0)**2
        if (maxJ > 0.0):
            xmod += jxmod * (time - globalt0)**3 / 6.0
            ymod += jymod * (time - globalt0)**3 / 6.0

        lnX = -0.5 * (xp - xmod)**2 / xerr**2
        lnY = -0.5 * (yp - ymod)**2 / yerr**2

        return np.sum(lnX) + np.sum(lnY)


    n_params = 6
    if (maxJ > 0.0):
        n_params = 8

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def PMN_xvz(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
            poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
            globalt0 = 2013.318,maxR=2.0,maxZ=100.0,maxV=1.0e4,maxJ=0.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxR in "
    #maxZ in mpc
    #maxV in km/s
    #maxJ in micro-as/yr^3

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_pc / 1000.
    #maxV *= 1.0e5
    #maxA *= asy_to_kms * 0.1 / sec_in_yr


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * 2.0 * maxR - maxR #x0 in cm
        cube[1] = cube[1] * 2.0 * maxR - maxR #y0 in cm
        cube[2] = cube[2] * maxZ #z0 in cm
        cube[3] = cube[3] * 2.0 * maxV - maxV #vx in cm/s
        cube[4] = cube[4] * 2.0 * maxV - maxV #vy in cm/s
        if (maxJ > 0.0):
            cube[5] = cube[5] * 2.0 * maxJ - maxJ #x and y jerk terms
            cube[6] = cube[6] * 2.0 * maxJ - maxJ
        

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        x0mod = cube[0] * dist * cm_in_au
        y0mod = cube[1] * dist * cm_in_au
        z0mod = cube[2] * cm_in_pc / 1000.0
        vxmod = cube[3] * 1.0e5
        vymod = cube[4] * 1.0e5
        if (maxJ > 0.0):
            jxmod = cube[5] * asy_to_kms * 0.1 / sec_in_yr**2
            jymod = cube[6] * asy_to_kms * 0.1 / sec_in_yr**2
            
            total_lnL = star_lnL(star,x0mod,y0mod,z0mod,vxmod,vymod,jxmod=jxmod,jymod=jymod)
        else:
            total_lnL = star_lnL(star,x0mod,y0mod,z0mod,vxmod,vymod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,x0mod,y0mod,z0mod,vxmod,vymod,t_0=t_0,globalt0=globalt0,jxmod=0.0,jymod=0.0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        amod = -1.0 * GM  / (math.sqrt(x0mod**2 + y0mod**2 + z0mod**2))**3

        xmod = x0mod + vxmod*(time - globalt0) + 0.5 * amod * xp * (time - globalt0)**2
        ymod = y0mod + vymod*(time - globalt0) + 0.5 * amod * yp * (time - globalt0)**2
        if (maxJ > 0.0):
            xmod += jxmod * (time - globalt0)**3 / 6.0
            ymod += jymod * (time - globalt0)**3 / 6.0

        lnX = -0.5 * (xp - xmod)**2 / xerr**2
        lnY = -0.5 * (yp - ymod)**2 / yerr**2

        return np.sum(lnX) + np.sum(lnY)


    n_params = 5
    if (maxJ > 0.0):
        n_params = 7

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')







def PMN_vz_forcexy(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                   poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
                   maxZ=100.0,maxV=1.0e4,maxJ=0.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxR in "
    #maxZ in mpc
    #maxV in km/s
    #maxJ in micro-as/yr^3

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_pc / 1000.
    #maxV *= 1.0e5
    #maxA *= asy_to_kms * 0.1 / sec_in_yr


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * maxZ #z0 in cm
        cube[1] = cube[1] * 2.0 * maxV - maxV #vx in cm/s
        cube[2] = cube[2] * 2.0 * maxV - maxV #vy in cm/s
        if (maxJ > 0.0):
            cube[3] = cube[3] * 2.0 * maxJ - maxJ #x and y jerk terms
            cube[4] = cube[4] * 2.0 * maxJ - maxJ
        

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        z0mod = cube[0] * cm_in_pc / 1000.0
        vxmod = cube[1] * 1.0e5
        vymod = cube[2] * 1.0e5
        if (maxJ > 0.0):
            jxmod = cube[3] * asy_to_kms * 0.1 / sec_in_yr**2
            jymod = cube[4] * asy_to_kms * 0.1 / sec_in_yr**2
            
            total_lnL = star_lnL(star,z0mod,vxmod,vymod,jxmod=jxmod,jymod=jymod)
        else:
            total_lnL = star_lnL(star,z0mod,vxmod,vymod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,z0mod,vxmod,vymod,t_0=t_0,jxmod=0.0,jymod=0.0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)

        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        tmpdex = np.where(time > 2007.0)[0]
        forceDex = np.median(tmpdex)
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = time[forceDex]
        x0mod = xp[forceDex]
        y0mod = yp[forceDex]

        amod = -1.0 * GM  / (math.sqrt(x0mod**2 + y0mod**2 + z0mod**2))**3

        xmod = x0mod + vxmod*(time - globalt0) + 0.5 * amod * xp * (time - globalt0)**2
        ymod = y0mod + vymod*(time - globalt0) + 0.5 * amod * yp * (time - globalt0)**2
        if (maxJ > 0.0):
            xmod += jxmod * (time - globalt0)**3 / 6.0
            ymod += jymod * (time - globalt0)**3 / 6.0

        lnX = -0.5 * (xp - xmod)**2 / xerr**2
        lnY = -0.5 * (yp - ymod)**2 / yerr**2

        return (np.sum(lnX) + np.sum(lnY))


    n_params = 3
    if (maxJ > 0.0):
        n_params = 5

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')




def PMN_vel_z(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
              poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
              globalt0 = 2013.318,maxR=2.0,maxZ=100.0,maxV=1.0e4,maxJ=0.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxR in "
    #maxZ in mpc
    #maxV in km/s
    #maxJ in micro-as/yr^3

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_pc / 1000.
    #maxV *= 1.0e5
    #maxA *= asy_to_kms * 0.1 / sec_in_yr


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * maxZ #z0 in cm
        cube[1] = cube[1] * 2.0 * maxV - maxV #vx in cm/s
        cube[2] = cube[2] * 2.0 * maxV - maxV #vy in cm/s
        if (maxJ > 0.0):
            cube[3] = cube[3] * 2.0 * maxJ - maxJ #x and y jerk terms
            cube[4] = cube[4] * 2.0 * maxJ - maxJ
        

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        z0mod = cube[0] * cm_in_pc / 1000.0
        vxmod = cube[1] * 1.0e5
        vymod = cube[2] * 1.0e5
        if (maxJ > 0.0):
            jxmod = cube[3] * asy_to_kms * 0.1 / sec_in_yr**2
            jymod = cube[4] * asy_to_kms * 0.1 / sec_in_yr**2
            
            total_lnL = star_lnL(star,z0mod,vxmod,vymod,jxmod=jxmod,jymod=jymod)
        else:
            total_lnL = star_lnL(star,z0mod,vxmod,vymod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,z0mod,vxmod,vymod,t_0=t_0,globalt0=globalt0,jxmod=0.0,jymod=0.0,maxR=maxR):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)

        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr

        def like_toSum(x0mod,y0mod,z0mod=z0mod,vxmod=vxmod,vymod=vymod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,
                       globalt0=globalt0,jxmod=jxmod,jymod=jymod):

            amod = -1.0 * GM  / (math.sqrt(x0mod**2 + y0mod**2 + z0mod**2))**3

            xmod = x0mod + vxmod*(time - globalt0) + 0.5 * amod * xp * (time - globalt0)**2
            ymod = y0mod + vymod*(time - globalt0) + 0.5 * amod * yp * (time - globalt0)**2
            if (maxJ > 0.0):
                xmod += jxmod * (time - globalt0)**3 / 6.0
                ymod += jymod * (time - globalt0)**3 / 6.0

            lnX = -0.5 * (xp - xmod)**2 / xerr**2
            lnY = -0.5 * (yp - ymod)**2 / yerr**2

            return math.exp(np.sum(lnX) + np.sum(lnY))

        def int_y(maxR=maxR):
            return [-1.0*maxR, maxR]

        def int_x(y0mod,maxR=maxR):
            return [-1.0*maxR,maxR]

        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        return math.log(starlike)


    n_params = 3
    if (maxJ > 0.0):
        n_params = 5

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')




def PMN_intv(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
             poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
             globalt0 = 2013.318,maxR=2.0,maxZ=100.0,maxV=1.0e4,maxJ=0.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxR in "
    #maxZ in mpc
    #maxV in km/s
    #maxJ in micro-as/yr^3

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_pc / 1000.
    #maxV *= 1.0e5
    #maxA *= asy_to_kms * 0.1 / sec_in_yr


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * 2.0 * maxR - maxR #x0 in cm
        cube[1] = cube[1] * 2.0 * maxR - maxR #y0 in cm
        cube[2] = cube[2] * maxZ #z0 in cm
        #cube[3] = cube[3] * 2.0 * maxV - maxV #vx in cm/s
        #cube[4] = cube[4] * 2.0 * maxV - maxV #vy in cm/s
        if (maxJ > 0.0):
            cube[3] = cube[3] * 2.0 * maxJ - maxJ #x and y jerk terms
            cube[4] = cube[4] * 2.0 * maxJ - maxJ
        

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        x0mod = cube[0] * dist * cm_in_au
        y0mod = cube[1] * dist * cm_in_au
        z0mod = cube[2] * cm_in_pc / 1000.0
        #vxmod = cube[3] * 1.0e5
        #vymod = cube[4] * 1.0e5
        if (maxJ > 0.0):
            jxmod = cube[3] * asy_to_kms * 0.1 / sec_in_yr**2
            jymod = cube[4] * asy_to_kms * 0.1 / sec_in_yr**2
            
            total_lnL = star_lnL(star,x0mod,y0mod,z0mod,jxmod=jxmod,jymod=jymod)
        else:
            total_lnL = star_lnL(star,x0mod,y0mod,z0mod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,x0mod,y0mod,z0mod,t_0=t_0,globalt0=globalt0,jxmod=0.0,jymod=0.0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        globalt0 = globalt0 * sec_in_yr
        amod = -1.0 * GM  / (math.sqrt(x0mod**2 + y0mod**2 + z0mod**2))**3

        

        xtilda = xp - (x0mod + 0.5 * amod * xp * (time - globalt0)**2)
        ytilda = yp - (y0mod + 0.5 * amod * yp * (time - globalt0)**2)
        if (maxJ > 0.0):
            xtilda -= jxmod * (time - globalt0)**3 / 6.0
            ytilda -= jymod * (time - globalt0)**3 / 6.0

        lnX = -0.5 * (np.sum(xtilda**2 / xerr**2) - (np.sum(xtilda*(time - globalt0)/xerr**2))**2 * np.sum(xerr**2/(time - globalt0)**2))
        lnY = -0.5 * (np.sum(ytilda**2 / yerr**2) - (np.sum(ytilda*(time - globalt0)/yerr**2))**2 * np.sum(yerr**2/(time - globalt0)**2))

        lnL_pos = (-1.0*np.sum(xtilda**2/(2.0*xerr**2))+(np.sum(xtilda*(time-globalt0)/(2.0*xerr**2)))**2/
                    np.sum((time-globalt0)**2/(2.0*xerr**2))-
                    1.0*np.sum(ytilda**2/(2.0*yerr**2))+(np.sum(ytilda*(time-globalt0)/(2.0*yerr**2)))**2/
                    np.sum((time-globalt0)**2/(2.0*yerr**2)))

        return lnL_pos


    n_params = 3
    if (maxJ > 0.0):
        n_params = 5

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')




def PMN_intvxy(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
               poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
               maxZ=100.0,maxJ=0.0,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxR in "
    #maxZ in mpc
    #maxV in km/s
    #maxJ in micro-as/yr^3

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_pc / 1000.
    #maxV *= 1.0e5
    #maxA *= asy_to_kms * 0.1 / sec_in_yr


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        #cube[0] = cube[0] * 2.0 * maxR - maxR #x0 in cm
        #cube[1] = cube[1] * 2.0 * maxR - maxR #y0 in cm
        cube[0] = cube[0] * maxZ #z0 in cm
        #cube[3] = cube[3] * 2.0 * maxV - maxV #vx in cm/s
        #cube[4] = cube[4] * 2.0 * maxV - maxV #vy in cm/s
        if (maxJ > 0.0):
            cube[1] = cube[1] * 2.0 * maxJ - maxJ #x and y jerk terms
            cube[2] = cube[2] * 2.0 * maxJ - maxJ
        

    #pdb.set_trace()
    #was getting a division by zero with dblquad integrating over all x and y's tested
    #integrate over posative quad and multiply result by 4
    def SPLloglike(cube,ndim,nparams):
        #x0mod = cube[0] * dist * cm_in_au
        #y0mod = cube[1] * dist * cm_in_au
        z0mod = cube[0] * cm_in_pc / 1000.0
        #vxmod = cube[3] * 1.0e5
        #vymod = cube[4] * 1.0e5
        if (maxJ > 0.0):
            jxmod = cube[1] * asy_to_kms * 0.1 / sec_in_yr**2
            jymod = cube[2] * asy_to_kms * 0.1 / sec_in_yr**2
            
            total_lnL = star_lnL(star,z0mod,jxmod=jxmod,jymod=jymod)
        else:
            total_lnL = star_lnL(star,z0mod)

#        if (total_lnL == 0.0):
#            total_lnL = 1e-323
        return total_lnL

    tmpName = star

    def star_lnL(starName,z0mod,t_0=t_0,jxmod=0.0,jymod=0.0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        time = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((time - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((time - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        time = time * sec_in_yr
        t_0 = t_0 * sec_in_yr
        forceDex = len(xp)-1
        globalt0 = time[forceDex]
        minX = xp[forceDex] - 3.0*xerr[forceDex]
        maxX = xp[forceDex] + 3.0*xerr[forceDex]
        minY = yp[forceDex] - 3.0*yerr[forceDex]
        maxY = yp[forceDex] + 3.0*yerr[forceDex]


        def like_toSum(x0mod,y0mod,z0mod=z0mod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,time=time,globalt0=globalt0,jxmod=jxmod,jymod=jymod):
            amod = -1.0 * GM  / (math.sqrt(x0mod**2 + y0mod**2 + z0mod**2))**3
            xtilda = xp - (x0mod + 0.5 * amod * xp * (time - globalt0)**2)
            ytilda = yp - (y0mod + 0.5 * amod * yp * (time - globalt0)**2)
            if ((jxmod != 0.0) | (jymod != 0.0)):
                xtilda -= jxmod * (time - globalt0)**3 / 6.0
                ytilda -= jymod * (time - globalt0)**3 / 6.0

            lnL_pos = (-1.0*np.sum(xtilda**2/(2.0*xerr**2))+(np.sum(xtilda*(time-globalt0)/(2.0*xerr**2)))**2/
                        np.sum((time-globalt0)**2/(2.0*xerr**2))-
                        1.0*np.sum(ytilda**2/(2.0*yerr**2))+(np.sum(ytilda*(time-globalt0)/(2.0*yerr**2)))**2/
                        np.sum((time-globalt0)**2/(2.0*yerr**2)))

            return math.exp(lnL_pos)

        def int_y(minY=minY,maxY=maxY):
            return [minY,maxY]

        def int_x(y0mod,minX=minX,maxX=maxX):
            return [minX,maxX]

        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]
        if (starlike == 0.0):
            starlike = 1.0e-323
        return math.log(starlike)


    n_params = 1
    if (maxJ > 0.0):
        n_params = 3

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def PMN_intvjxyMIN(star,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                   poly='polyfit_nz/fit',points='points_nz/',chainsDir='efit/chains_S0-2_newRV2/',
                   maxZ=100.0,globalt0=2013.318,flag='brokenPower'):
    #Uses PyMultiNest to determine gamma
    #maxA in micro-as/yr^2
    #maxZ in mpc


    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    #Define priors
    #single power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * maxZ
        

    def SPLloglike(cube,ndim,nparams):
        z0mod = cube[0] * 1.0e-3 * cm_in_pc

        total_lnL = star_lnL(star,z0mod)

        return total_lnL

    tmpName = star

    def star_lnL(starName,z0mod,t_0=t_0,globalt0=globalt0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + tmpName + '.points')

        timep = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((timep - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((timep - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        timep *= sec_in_yr
        t_0 *= sec_in_yr
        globalt0 *= sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((timep - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((timep - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((timep - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((timep - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((timep - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((timep - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((timep - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((timep - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)

        def like_toSum(x0mod,y0mod,z0mod=z0mod,xp=xp,yp=yp,
                       xerr=xerr,yerr=yerr,timep=timep,globalt0=globalt0,
                       inv_xf=inv_xf,inv_yf=inv_yf,toMinimize=False):

            amod = -1.0 * GM / (math.sqrt(z0mod**2 + y0mod**2 +x0mod**2))**3
            xtilda = xp - (x0mod + 0.5 * amod * x0mod * (timep - globalt0)**2)
            ytilda = yp - (y0mod + 0.5 * amod * y0mod * (timep - globalt0)**2)

            axarray = np.zeros(2)
            ayarray = np.zeros(2)
            axarray[0] = np.sum(xtilda*(timep - globalt0)/xerr**2.0)
            axarray[1] = np.sum(xtilda*(timep - globalt0)**3.0/xerr**2.0)
            ayarray[0] = np.sum(ytilda*(timep - globalt0)/yerr**2.0)
            ayarray[1] = np.sum(ytilda*(timep - globalt0)**3.0/yerr**2.0)

            x_af = np.dot(axarray,inv_xf)
            y_af = np.dot(ayarray,inv_yf)

            lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                       np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0


            if (toMinimize==True):
                return lnL_pos*-1.0
            else:
                return math.exp(lnL_pos)

        def lnL_toMin(xy):

            x0mod = xy[0]
            y0mod = xy[1]
            neg_lnL = like_toSum(x0mod,y0mod,toMinimize=True)

            return neg_lnL


        res = minimize(lnL_toMin,[xp[-1],yp[-1]],method='nelder-mead')
        xy_maxL = res.x

        lnLmin = lnL_toMin(xy_maxL)

        minX = xy_maxL[0] * 1.0
        maxX = xy_maxL[0] * 1.0
        minY = xy_maxL[1] * 1.0
        maxY = xy_maxL[1] * 1.0
        delta_xy = 1.0e-4 * dist * cm_in_au


        tmp_lnL = lnLmin * 1.0
        while (tmp_lnL <= (math.log(1.0e4) + lnLmin)):
            minX -= delta_xy
            tmp_lnL = lnL_toMin([minX,xy_maxL[1]])

        tmp_lnL = lnLmin * 1.0
        while (tmp_lnL <= (math.log(1.0e4) + lnLmin)):
            minY -= delta_xy
            tmp_lnL = lnL_toMin([xy_maxL[0],minY])

        tmp_lnL = lnLmin * 1.0
        while (tmp_lnL <= (math.log(1.0e4) + lnLmin)):
            maxY += delta_xy
            tmp_lnL = lnL_toMin([xy_maxL[0],maxY])

        tmp_lnL = lnLmin * 1.0
        while (tmp_lnL <= (math.log(1.0e4) + lnLmin)):
            maxX += delta_xy
            tmp_lnL = lnL_toMin([maxX,xy_maxL[1]])


        def int_y(minY=minY,maxY=maxY):
            return [minY,maxY]

        def int_x(y0mod,minX=minX,maxX=maxX):
            return [minX,maxX]


        starlike = integrate.nquad(like_toSum,[int_x,int_y])
        starlike = starlike[0]


        if (starlike == 0.0):
            starlike = 1.0e-323
        return math.log(starlike)


    n_params = 1

    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='pmnOld/'+flag+'_')





def gamma_broken(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                 poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains_S0-2_newRV2/',
                 Rcut=1.7,nEpochs=14.,magCut=15.5,grange=[-5.0,1.9],maxBreak=5.0,alrange=[-10.0,10.0],
                 maxDelta=10.0,maxZ=1.0e20,globalt0=2006.0,maxR=2.0,flag='brokenPower',makeMap=False):
    #Uses PyMultiNest to determine gamma
    #maxBreak in pc
    #maxZ in arcsec


    maxR *= dist * cm_in_au
    maxZ *= dist * cm_in_au

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')
    xfcut = s.getArray('x0')
    yfcut = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut = np.where((mag < magCut) & (cnt > nEpochs) & ((xfcut**2.0 + yfcut**2.0) <= Rcut**2.0))[0]
    names = [names[nn] for nn in cut]

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
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
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


    #Define priors
    #broken power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        cube[1] = cube[1] * (alrange[1] - alrange[0]) + alrange[0] #alpha
        cube[2] = cube[2] * maxDelta #delta
        cube[3] = cube[3] * maxBreak #r break
#        for i in range(len(pdex)):
#            cube[i+4] = cube[i+4] * maxZ  #z0
        

    def SPLloglike(cube,ndim,nparams,maxR=maxR,maxZ=maxZ):
        print 'Begin loglike func: '+str(datetime.datetime.now())
        gmod = cube[0]
        almod = cube[1]
        demod = cube[2]
        brmod = cube[3] * cm_in_pc
        #z0mod = cube[0] * 1.0e-3 * cm_in_pc
        print gmod,almod,demod,brmod

        def density_norm(Rprime,zprime,gmod=gmod,almod=almod,demod=demod,brmod=brmod):
            return (Rprime*(Rprime**2+zprime**2)**(gmod/-2.0)*
                    (1.0+(math.sqrt(Rprime**2+zprime**2)/brmod)**demod)**((gmod-almod)/demod))

        norm = integrate.nquad(density_norm,[[0.0,maxR],[0.0,maxZ]])
        norm = norm[0]
        if (norm <= 0.0):
            maxZ = 5.0 * brmod * (abs(gmod/almod))**(1.0/demod)
            norm = integrate.nquad(density_norm,[[0.0,maxR],[0.0,maxZ]])
            norm = norm[0]

            if (norm <= 0.0):
                print 'Density normalization is still zero'
                pdb.set_trace()

        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            #print names[p]
            tmplnL = star_lnL(names[p],gmod,almod,demod,brmod,norm)
            total_lnL += tmplnL*oldProb[p]

        #print 'End loglike func'
        #pdb.set_trace()
        return total_lnL


#    def star_lnL(starName,gmod,almod,demod,brmod,z0mod,norm,t_0=t_0,globalt0=globalt0):
    def star_lnL(starName,gmod,almod,demod,brmod,norm,t_0=t_0,globalt0=globalt0):
        pointsTab = asciidata.open(root_tmp + alnDir + points + starName + '.points')

        #pdb.set_trace()
        print starName


        timep = pointsTab[0].tonumpy()
        xp = pointsTab[1].tonumpy() * dist * cm_in_au
        yp = pointsTab[2].tonumpy() * dist * cm_in_au
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        #print math.sqrt(np.mean(xp)**2 + np.mean(yp)**2) / dist / cm_in_au
        #print str(datetime.datetime.now())

        if (updateErr == True):
            xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((timep - t_0)*ori_vxe)**2)
            yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((timep - t_0)*ori_vye)**2)


        xerr = xerr * dist * cm_in_au
        yerr = yerr * dist * cm_in_au
        timep *= sec_in_yr
        t_0 *= sec_in_yr
        globalt0 *= sec_in_yr

        p_error2 = 1.0e300
        xfarray = np.zeros([2,2])
        xfarray[0,0] = np.sum((timep - globalt0)**2.0/(xerr**2.0)) + 1.0/p_error2
        xfarray[0,1] = np.sum((timep - globalt0)**4.0/(xerr**2.0))
        xfarray[1,0] = np.sum((timep - globalt0)**4.0/(xerr**2.0))
        xfarray[1,1] = np.sum((timep - globalt0)**6.0/(xerr**2.0)) + 1.0/p_error2

        inv_xf = np.linalg.inv(xfarray)
        
        yfarray = np.zeros([2,2])
        yfarray[0,0] = np.sum((timep - globalt0)**2.0/(yerr**2.0)) + 1.0/p_error2
        yfarray[0,1] = np.sum((timep - globalt0)**4.0/(yerr**2.0))
        yfarray[1,0] = np.sum((timep - globalt0)**4.0/(yerr**2.0))
        yfarray[1,1] = np.sum((timep - globalt0)**6.0/(yerr**2.0)) + 1.0/p_error2
        
        inv_yf = np.linalg.inv(yfarray)

#        def like_toSum(x0mod,y0mod,z0mod=z0mod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,
#                       timep=timep,globalt0=globalt0,inv_xf=inv_xf,inv_yf=inv_yf,
#                       gmod=gmod,almod=almod,demod=demod,brmod=brmod,norm=norm,
#                       toMinimize=False):
        def like_toSum(x0mod,y0mod,z0mod,xp=xp,yp=yp,xerr=xerr,yerr=yerr,
                       timep=timep,globalt0=globalt0,inv_xf=inv_xf,inv_yf=inv_yf,
                       gmod=gmod,almod=almod,demod=demod,brmod=brmod,norm=norm,
                       toMinimize=False):

            amod = -1.0 * GM / (math.sqrt(z0mod**2 + y0mod**2 +x0mod**2))**3
            xtilda = xp - (x0mod + 0.5 * amod * x0mod * (timep - globalt0)**2)
            ytilda = yp - (y0mod + 0.5 * amod * y0mod * (timep - globalt0)**2)

            axarray = np.zeros(2)
            ayarray = np.zeros(2)
            axarray[0] = np.sum(xtilda*(timep - globalt0)/xerr**2.0)
            axarray[1] = np.sum(xtilda*(timep - globalt0)**3.0/xerr**2.0)
            ayarray[0] = np.sum(ytilda*(timep - globalt0)/yerr**2.0)
            ayarray[1] = np.sum(ytilda*(timep - globalt0)**3.0/yerr**2.0)

            x_af = np.dot(axarray,inv_xf)
            y_af = np.dot(ayarray,inv_yf)

            lnL_pos = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af,axarray) +
                       np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af,ayarray))/-2.0

            like_den = ((x0mod**2+y0mod**2+z0mod**2)**(gmod/-2.0)*
                        (1.0+(math.sqrt(x0mod**2+y0mod**2+z0mod**2)/brmod)**demod)**((gmod-almod)/demod) / norm)

            if (toMinimize==True):
                return -1.0*(lnL_pos + math.log(like_den))
            else:
                like_pos = math.exp(lnL_pos)

                like_return = like_pos*like_den
                if ((like_return == 0.0) | (like_return != like_return)):
                    like_return = 1e-323
                if (like_return == float('Inf')):
                    print 'Got Inf value in the summation'
                    pdb.set_trace()

                return like_return

        #def lnL_toMin(xy):
        def lnL_toMin(xyz):

            #x0mod = xy[0]
            #y0mod = xy[1]
            #neg_lnL = like_toSum(x0mod,y0mod,toMinimize=True)
            x0mod = xyz[0]
            y0mod = xyz[1]
            z0mod = xyz[2]
            neg_lnL = like_toSum(x0mod,y0mod,z0mod,toMinimize=True)

            return neg_lnL

        #print 'Minimize'

        #res = minimize(lnL_toMin,[xp[-1],yp[-1]],method='nelder-mead')
        res = minimize(lnL_toMin,[xp[-1],yp[-1],1.0e16],method='nelder-mead')
        #xy_maxL = res.x
        xyz_maxL = res.x

#        lnLmin = lnL_toMin(xy_maxL)
        lnLmin = lnL_toMin(xyz_maxL)

        if (xyz_maxL[2] >= maxZ):
            print 'Max Z is not high enough'
            pdb.set_trace()

        #print 'Find bounds'

#        minX = xy_maxL[0] * 1.0
#        maxX = xy_maxL[0] * 1.0
#        minY = xy_maxL[1] * 1.0
#        maxY = xy_maxL[1] * 1.0
        minX = xyz_maxL[0] * 1.0
        maxX = xyz_maxL[0] * 1.0
        minY = xyz_maxL[1] * 1.0
        maxY = xyz_maxL[1] * 1.0
#        minZs = xyz_maxL[2] * 1.0
#        maxZs = xyz_maxL[2] * 1.0
        delta_xy = 1.0e-3 * dist * cm_in_au
#        delta_z = xyz_maxL[2] * 0.3

#        pdb.set_trace()

        tmp_lnL = lnLmin * 1.0
        while ((tmp_lnL <= (math.log(1.0e3) + lnLmin)) & (minX > (xyz_maxL[0] - 0.1*dist*cm_in_au))):
            minX -= delta_xy
            #tmp_lnL = lnL_toMin([minX,xy_maxL[1]])
            tmp_lnL = lnL_toMin([minX,xyz_maxL[1],xyz_maxL[2]])

        tmp_lnL = lnLmin * 1.0
        while ((tmp_lnL <= (math.log(1.0e3) + lnLmin)) & (minY > (xyz_maxL[1] - 0.1*dist*cm_in_au))):
            minY -= delta_xy
            #tmp_lnL = lnL_toMin([xy_maxL[0],minY])
            tmp_lnL = lnL_toMin([xyz_maxL[0],minY,xyz_maxL[2]])

        tmp_lnL = lnLmin * 1.0
        while ((tmp_lnL <= (math.log(1.0e3) + lnLmin)) & (maxY < (xyz_maxL[1] + 0.1*dist*cm_in_au))):
            maxY += delta_xy
            #tmp_lnL = lnL_toMin([xy_maxL[0],maxY])
            tmp_lnL = lnL_toMin([xyz_maxL[0],maxY,xyz_maxL[2]])

        tmp_lnL = lnLmin * 1.0
        while ((tmp_lnL <= (math.log(1.0e3) + lnLmin)) & (maxX < (xyz_maxL[0] + 0.1*dist*cm_in_au))):
            maxX += delta_xy
            #tmp_lnL = lnL_toMin([maxX,xy_maxL[1]])
            tmp_lnL = lnL_toMin([maxX,xyz_maxL[1],xyz_maxL[2]])

#        tmp_lnL = lnLmin * 1.0
#        while (tmp_lnL <= (math.log(1.0e3) + lnLmin)):
#            minZs -= delta_z
#            tmp_lnL = lnL_toMin([xyz_maxL[0],xyz_maxL[1],minZs])

#        tmp_lnL = lnLmin * 1.0
#        while (tmp_lnL <= (math.log(1.0e3) + lnLmin)):
#            maxZs += delta_z
#            tmp_lnL = lnL_toMin([xyz_maxL[0],xyz_maxL[1],maxZs])

#        def int_z(maxZ=maxZ):
#            return [0.0,maxZ]

        #def int_y(minY=minY,maxY=maxY):
#        def int_y(z0mod,minY=minY,maxY=maxY):
#            return [minY,maxY]

        #def int_x(y0mod,minX=minX,maxX=maxX):
#        def int_x(z0mod,y0mod,minX=minX,maxX=maxX):
#            return [minX,maxX]

        #print 'Integrate'

       # def fast_int_y(minY=minY,maxY=maxY):
       #     return [minY,maxY]

        #def fast_int_x(y0mod,minX=minX,maxX=maxX):
        #    return [minX,maxX]

        #pdb.set_trace()

        #starlike = integrate.nquad(like_toSum,[int_x,int_y])
#        starlike = integrate.nquad(like_toSum,[int_x,int_y,int_z])
#        starlike = starlike[0]

        #pdb.set_trace()

        def like_toSumX(xtries,minY=minY,maxY=maxY,maxZ=maxZ):
            #def fast_int_z(maxZ=maxZ):
            #    return [0.0,maxZ]
            xtmp=float(0.0) #placeholder

            def like_toSumY(ytries,xtmp=xtmp,maxZ=maxZ):

                ytmp=float(0.0) #placeholder

                def like_toSumZ(ztries,xtmp=xtmp,ytmp=ytmp,maxZ=maxZ):
                    like_forz = np.array([])
                    for ztmp in ztries:
                        tmplike = like_toSum(xtmp,ytmp,ztmp)
                        like_forz = np.append(like_forz,tmplike)
                    return like_forz
                #L_over_z = integrate.quadrature(like_toSumZ,0.0,maxZ,maxiter=40,tol=1.0e-5)
                #return L_over_z[0]

                like_fory = np.array([])
                for ytmp in ytries:
                    tmplikeY = integrate.quadrature(like_toSumZ,0.0,maxZ,maxiter=40,tol=1.0e-5)
                    like_fory = np.append(like_fory,tmplikeY[0])
                return like_fory

           # L_over_yz = integrate.quadrature(like_toSumY,minY,maxY,maxiter=40,tol=1.0e-5)
            #return L_over_yz[0]

            like_forx = np.array([])
            for xtmp in xtries:
                tmplikeX = integrate.quadrature(like_toSumY,minY,maxY,maxiter=40,tol=1.0e-5)
                like_forx = np.append(like_forx,tmplikeX[0])
            return like_forx


        starlike = integrate.quadrature(like_toSumX,minX,maxX,maxiter=40,tol=1.0e-5)
        starlike = starlike[0]

        #pdb.set_trace()

        if (starlike == 0.0):
            starlike = 1.0e-323
        if (starlike == float('Inf')):
            print 'Got Inf values in a star likelihood'
            pdb.set_trace()

        return math.log(starlike)


    n_params = 4 #len(pdex)+4


    progress = pymultinest.ProgressPrinter(n_params=n_params,outputfiles_basename='/u/schappell/pmnOld/'+flag+'_',interval_ms=100000)
    progress.start()
    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='/u/schappell/pmnOld/'+flag+'_',verbose=True)






def gamma_Accel(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
                Rcut=1.7,nEpochs=14.,magCut=15.5,grange=[-5.0,1.9],maxBreak=5.0,alrange=[0.0,10.0],
                maxDelta=10.0,maxZ=10.0,globalt0=2006.0,maxR=2.0,flag='brokenPower',polyj='polyfit_nzj/fit',
                pvalue=4.0):
    #Uses PyMultiNest to determine gamma
    #maxBreak in pc
    #maxZ in pc


    maxR *= dist * cm_in_au
    maxZ *= cm_in_pc

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    names, r2d, ar, are = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                    Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,
                                    chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs)

    min_accel = -1.0*GM / r2d**2
    max_accel = -1.0*GM*r2d / (np.sqrt(r2d**2 + maxZ**2))**3
    maxr = np.sqrt(maxR**2 + maxZ**2)
    maxz_star = np.sqrt(maxr**2 - r2d**2)

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


    #Define priors
    #broken power law, uniform priors for all
    def SPLprior(cube,ndim,nparams):
        cube[0] = cube[0] * (grange[1] - grange[0]) + grange[0] #gamma
        cube[1] = cube[1] * (alrange[1] - alrange[0]) + alrange[0] #alpha
        cube[2] = cube[2] * maxDelta #delta
        cube[3] = cube[3] * maxBreak #r break
#        for i in range(len(pdex)):
#            cube[i+4] = cube[i+4] * maxZ  #z0
        

    def SPLloglike(cube,ndim,nparams,maxR=maxR,maxZ=maxZ):
        print 'Begin loglike func: '+str(datetime.datetime.now())
        gmod = cube[0]
        almod = cube[1]
        demod = cube[2]
        brmod = cube[3] * cm_in_pc
        #z0mod = cube[0] * 1.0e-3 * cm_in_pc
        print gmod,almod,demod,brmod
        pdb.set_trace()


        def density_norm(Rprime,zprime,gmod=gmod,almod=almod,demod=demod,brmod=brmod):
            return (Rprime*(Rprime**2+zprime**2)**(gmod/-2.0)*
                    (1.0+(math.sqrt(Rprime**2+zprime**2)/brmod)**demod)**((gmod-almod)/demod))

        norm = integrate.nquad(density_norm,[[0.0,maxR],[0.0,maxZ]])
        norm = norm[0]

        #def density_SumR(Rtries,maxZ=maxZ):
        #    Rtmp=float(0.0)#placeholder

        #    def density_SumZ(ztries,Rtmp=Rtmp,maxZ=maxZ):
        #        density_sumz = np.array([])
        #        for ztmp in ztries:
        #            tmpsumz = density_norm(Rtmp,ztmp)
        #            density_sumz = np.append(density_sumz,tmpsumz)
        #        return density_sumz

        #    density_sum = np.array([])
        #    for Rtmp in Rtries:
        #        tmpsum = integrate.quadrature(density_SumZ,0.0,maxZ)
        #        density_sum = np.append(density_sum,tmpsum[0])
        #    return density_sum

        #norm = integrate.quadrature(density_SumR,0.0,maxR)
        #norm = norm[0]

        if (norm <= 0.0):
            maxZ = 5.0 * brmod * (abs(gmod/almod))**(1.0/demod)
            norm = integrate.nquad(density_norm,[[0.0,maxR],[0.0,maxZ]])
            norm = norm[0]
            print 'norm was 0 or less originally'

            if (norm <= 0.0):
                print 'Density normalization is still zero'
                pdb.set_trace()

        #pdb.set_trace()
        total_lnL = 0.0
        for i in range(len(pdex)):
            p = pdex[i]
            print names[p]
            tmplnL = star_lnL(p,gmod,almod,demod,brmod,norm)
            total_lnL += tmplnL*oldProb[p]

        #print 'End loglike func'
        pdb.set_trace()
        return total_lnL


#    def star_lnL(starName,gmod,almod,demod,brmod,z0mod,norm,t_0=t_0,globalt0=globalt0):
    def star_lnL(pdex,gmod,almod,demod,brmod,norm,t_0=t_0,globalt0=globalt0):

       # print str(datetime.datetime.now())
       
        def like_toSum(z0mod,pdex=pdex,gmod=gmod,almod=almod,demod=demod,brmod=brmod,norm=norm):

            #pdb.set_trace()
            if (ar[pdex] < 0.0):
                amod = -1.0 * GM * r2d[pdex] / (np.sqrt(r2d[pdex]**2 +z0mod**2))**3
                like_pos = np.exp(-1.0*(ar[pdex] - amod)**2/(2.0*are[pdex]**2))
                norm_pos = abs(are[pdex]*math.sqrt(pi/2.0)*(special.erf((ar[pdex] - min_accel[pdex])/(math.sqrt(2)*are[pdex]))
                                                          - special.erf((ar[pdex] - max_accel[pdex])/(math.sqrt(2)*are[pdex]))))
                if (norm_pos == 0.0):
                    norm_pos = integrate.quad(lambda x: math.exp(-1.0*(ar[pdex] - x)**2/(2.0*are[pdex]**2)), min_accel[pdex], max_accel[pdex])
                    norm_pos = abs(norm_pos[0])

            else:
                like_pos = (r2d[pdex]**2 + z0mod**2)**(2.5)/(3.0*z0mod*GM*r2d[pdex])
                norm_pos = maxz_star[pdex]

            like_den = ((r2d[pdex]**2+z0mod**2)**(gmod/-2.0)*
                        (1.0+(np.sqrt(r2d[pdex]**2+z0mod**2)/brmod)**demod)**((gmod-almod)/demod))


            like_return = like_pos*like_den/ (norm_pos * norm)

            if ((like_return==0.0) | (like_return != like_return)):
                like_return = 1e-323
            #for ii in range(len(like_return)):
            #    if ((like_return[ii] == 0.0) | (like_return[ii] != like_return[ii])):
            #        like_return[ii] = 1e-323
            #    if (like_return[ii] == float('Inf')):
            #        print 'Got Inf value in summation'
            #        pdb.set_trace()

                #if (like_return[ii] < 0.0):
            #pdb.set_trace()
            return like_return


        #starlike = integrate.quadrature(like_toSum,0.0,maxZ,maxiter=40)
        starlike = integrate.quad(like_toSum,0.0,maxZ)
        starlike = starlike[0]

        #pdb.set_trace()

        if (starlike == 0.0):
            starlike = 1.0e-323
        if (starlike == float('Inf')):
            print 'Got Inf value in star likelihood'
            pdb.set_trace()

        pdb.set_trace()
        return math.log(starlike)

    n_params = 4 #len(pdex)+4


    progress = pymultinest.ProgressPrinter(n_params=n_params,outputfiles_basename='/u/schappell/pmnOld/'+flag+'_',interval_ms=30000)
    progress.start()
    pymultinest.run(SPLloglike,SPLprior,n_params,outputfiles_basename='/u/schappell/pmnOld/'+flag+'_',verbose=True)





def gamma_C_accel(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                  poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',innerCut=5.0,
                  Rcut=1.7,nEpochs=14.,magCut=15.5,lowMag=-10.0,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,outerCut=15.0,
                  mnAlpha=6.0,mnDelta=4.0,mnBreak=0.5,max_r=5.0,schodel=False,maserDate='13jullgs1',situation=1,
                  label='',nonRadial=0,ak_correct=True,magCutMaser=17.75,lowMagMaser=-10.0,single=False,testing=False):

    #situation = 1, alpha, delta, and r_break free, stars in maser mosaic not included
               # 2, only gamma free, stars in maser mosaic not included
               # 3, alpha, delta, and r_break free, stars in maser mosaic included
               # 4, only gamma free, stars in maser mosaic included

    if (testing ==False):
        if ((situation != 1) & (situation != 2) & (situation != 3) & (situation != 4)):
            sys.exit("Did not pick available situation, try 1, 2, 3, or 4")
            
        if (schodel==True):
            ak_correct=True
        else:
            ak_correct=False #if using schodel's data, which is ext corrected, have to correct our data
        #otherwise, using stars from maser field, don't ext correct anything

    #Uses MultiNest in c to determine gamma
    #maxBreak in pc
    #maxZ in arcsec


    #maxR *= dist * cm_in_au
    #maxZ *= dist * cm_in_au

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
        names, r2d, ar, are = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                        Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                        chainsDir=chainsDir,starlist='all',magCut=magCut,lowMag=lowMag,nEpochs=nEpochs,nonRadial=nonRadial)

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
            if ((tmpName == 'S0-38') | (tmpName == 'S0-49') | (tmpName == 'S0-35') | (tmpName == 'S1-32')):
                oldProb[i] = 1.0
                if (tmpName == 'S0-61'):
                    oldProb[i] = 0.0

        for i in range(len(names)):
            tmpName = str(names[i])
            try:
                cur.execute('SELECT name,young,old FROM stars WHERE name=?', [tmpName])
                for row in cur:
                    if (row[1] == 'F') | (row[2] == 'T'):
                        oldProb[i] = 1.0
                    elif ((row[1] == 'T') | (row[2] == 'F')):
                        oldProb[i] = 0.0
                    
            except:
                continue

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
        np.savetxt('/u/schappell/code/c/stars_mn'+label+'.dat',np.transpose([r2d,ar,are,oldProb]),delimiter=' ')

    #Rcut *= dist / au_in_pc #Rcut now in pc
    #pdb.set_trace()

    mnRoot = '/u/schappell/pmnOld/data'+label
    if (single==False):
        if ((situation == 2) | (situation == 4)):
            mnRoot += '_fixed_a'+str(mnAlpha)+'_d'+str(mnDelta)+'_br'+str(mnBreak)+'pc'
        else:
            mnRoot += '_free'
        if (situation > 2):
            mnRoot += '_maser_cuts_'+str(innerCut)+'_'+str(outerCut)+'arcsec'
        else:
            mnRoot += '_cut_'+str(Rcut)+'arcsec'
    else:
        mnRoot += '_SPL_cut_'+str(Rcut)+'arcsec'

    mnRoot += '_'

    Rcut *= dist / au_in_pc
    if ((situation > 2) & (single==False)):
        maserStars(maserDate=maserDate,innerCut=innerCut,outerCut=outerCut,magCut=magCutMaser,lmagCut=lowMagMaser,schodel=schodel,label=label)
        outerCut *= dist / au_in_pc
        innerCut *= dist / au_in_pc

    if (single==False):
        os.system('g++ sc_mn_mthread.cpp gauss_legendre.c -o sc_mn -std=c++11 -lpthread -I/u/schappell/code/c -I/u/schappell/code/c/boost/config -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest')
        #os.system('g++ sc_mn_mthread.cpp gauss_legendre.c -o sc_mn -std=c++11 -lpthread -I/u/schappell/code/c -I/u/schappell/code/c/boost/config -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -L/u/schappell/code/c -lmultinest')

        outCommand = './sc_mn '+mnRoot+' '+str(mnAlpha)+' '+str(mnDelta)+' '+str(mnBreak)+' '+str(max_r)+' '+str(Rcut)
        outCommand += ' '+str(situation)+' '+str(innerCut)+' '+str(outerCut)+' '+str(nonRadial)+' '+label

        pdb.set_trace()

    else:
       os.system('g++ sc_mn_single.cpp gauss_legendre.c -o sc_mn -I/u/schappell/code/c -I/u/schappell/code/c/boost/config -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest') 

       outCommand = './sc_mn '+mnRoot+' '+str(max_r)+' '+str(Rcut)+' '+str(nonRadial)+' '+label

    os.system(outCommand)



def testCMN_old(numStars=27,min_r=0.0,max_r=1.0,gamma=-1.0,alpha=4.0,delta=4.0,r_break=0.5,perturb=False,R2d_cut=1.7,
                label='test_r2donly',offset_al=0.0,offset_de=0.0,offset_br=0.0,solver=True,situation=3,innerCut=0.0,outerCut=15.0,
                numStarsM=2127,nonRadial=0):

    #max_r in pc
    #R2d_cut in arcsec
    #r_break in pc
    #scaling_a in pc
    #offset parameters ADD that offset to the true value (alpha, delta, etc) and gives that new value to MN
    #perturb adds random perturbation (+ or -) within one sigma to all accelerations
    #rhoModel sets which model for density used, broken=borken power law (gamma, alpha, delta, and r_break)
    #         dehnen=dehnen (gamma and scaling radius)

    #situation = 1, alpha, delta, and r_break free, stars in maser mosaic not included
               # 2, only gamma free, stars in maser mosaic not included
               # 3, alpha, delta, and r_break free, stars in maser mosaic included
               # 4, only gamma free, stars in maser mosaic included

    if ((situation != 1) & (situation != 2) & (situation != 3) & (situation != 4)):
        sys.exit("Did not pick available situation, try 1, 2, 3, or 4")
        

    if (solver==False):
    #if not using solver
        numt = 100000

    R2d_cut *= dist / au_in_pc #R2d_cut now in pc
    innerCut *= dist / au_in_pc
    outerCut *= dist / au_in_pc

    #read in actual stars
    actual = np.loadtxt('/u/schappell/code/c/stars_gcows_for_mock.dat')
    are_act = actual[:,2]
    pOld_act = actual[:,3]

    if (nonRadial>0):
        plate_scale=0.00995
        gcows_zero = np.array([1500.0,1500.0])
        gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')

    #density profile
    def density(r,gamma=gamma,alpha=alpha,delta=delta,r_break=r_break):
        try:
            (r/r_break)**(-1.0*gamma)*(1.0 + (r/r_break)**delta)**((gamma-alpha)/delta)*r**2
        except:
            pdb.set_trace()
        return (r/r_break)**(-1.0*gamma)*(1.0 + (r/r_break)**delta)**((gamma-alpha)/delta)*r**2

      #  elif (rhoModel=='dehnen'):
       #     try:
                #(3.0 - gamma)*r**2*scaling_a / (r**gamma*(r+scaling_a)**(4.0-gamma))
         #       r**(-1.0*gamma+2.0)*((r/scaling_a)+1.0)**(-4.0+gamma)
                #print scaling_a
         #   except:
          #      pdb.set_trace()
            #return (3.0 - gamma)*r**2*scaling_a / (r**gamma*(r+scaling_a)**(4.0-gamma))
          #  return r**(-1.0*gamma+2.0)*((r/scaling_a)+1.0)**(-4.0+gamma)

    if (solver==True):
        const1 = integrate.quad(density,min_r,max_r)
        const1 = const1[0]
    else:
    #if not using solver
        delta_r = (max_r - min_r)/numt
        test_r = np.array([delta_r * i for i in range(numt)])
        test_r += min_r
        if (test_r[0]==0.0):
            test_r[0] += 1.0e-100
        test_c = np.array([density(test_r[i]) for i in range(numt)])
        const1 = np.sum(test_c)
        inverse_c = np.zeros(numt)
        for i in range(numt):
            tdex = np.arange(i+1)
            inverse_c[i] = np.sum(test_c[tdex]) / const1

    #pdb.set_trace()

    def CDF(rprime,*needed):
        const1,yprime = needed
        value = np.array([])
        for val in rprime:
            if (val < min_r):
                value = np.append(value, -1.0)
            else:
                tmp = integrate.quad(density,min_r,val)
                value = np.append(value,abs(tmp[0]/const1 - yprime))
        return value

    #rand01 = np.random.rand(numStars)

    radius = np.zeros(numStars)
    R2d = np.zeros(numStars)
    randFinal = np.zeros(numStars)
    ii = 0
    while(ii < numStars):
        #pdb.set_trace()
        #yprime = rand01[i]

        rand01 = np.random.rand(1)
        rand_angle = math.acos(np.random.rand(1))

        if (solver==True):
            rand_r = -1.0
            while (rand_r < min_r):
                rand_r = fsolve(CDF, 0.1,args=(const1,rand01)) #in pc
                if (rand_r==0.1):
                    tmp = minimize(CDF,0.1,args=(const1,rand01))
                    rand_r = tmp.x
                    if (rand_r==0.1):
                        tmp = minimize(CDF,0.1,method='nelder-mead',args=(const1,rand01),options={'maxiter':1000})
                        rand_r = tmp.x
                        if (rand_r==0.1):
                            pdb.set_trace()
        else: 
        #if not running solver
            tmpdex = np.argmin(abs(rand01 - inverse_c))
            rand_r = test_r[tmpdex]

        rand_R2d = np.sin(rand_angle) * rand_r
        if (rand_R2d <= R2d_cut):
            if (nonRadial==0):
                radius[ii] = rand_r
                R2d[ii] = rand_R2d
                randFinal[ii] = rand01
                print ii
                ii += 1
            else:
                rand_angle2 = math.acos(np.random.rand(1)*2.0 - 1.0)
                rand_x_pixels = (rand_R2d * math.cos(rand_angle2))/plate_scale + gcows_zero[0]
                rand_y_pixels = (rand_R2d * math.sin(rand_angle2))/plate_scale + gcows_zero[1]
                if((rand_x_pixels >= 0.0) & (rand_y_pixels >= 0.0) & (gcows[rand_x_pixels,rand_y_pixels]>0.0)):
                    radius[ii] = rand_r
                    R2d[ii] = rand_R2d
                    randFinal[ii] = rand01
                    print ii
                    ii += 1
        #print ii

    #pdb.set_trace()
    radius *= cm_in_pc #in cm
    R2d *= cm_in_pc #in cm
    #pdb.set_trace()

    #WRONG!!!!R2d = np.cos(np.random.rand(numStars)*pi/2.0) * radius
    accel = -1.0 * GM * R2d / radius**3 #in cm/s^2

    randex = np.random.randint(len(are_act),size=numStars)
    are = are_act[randex]
    pOld = pOld_act[randex]

    if (perturb == True):
        pert = np.random.rand(numStars)*2.0 - 1.0
        accel = accel + pert*are #perturb true accelerations within one sigma
    
    np.savetxt('/u/schappell/code/c/stars_mn.dat',np.transpose([R2d,accel,are,pOld]),delimiter=' ')

    #cutUse = R2d_cut
    if (situation > 2):
        #cutUse = outerCut * dist / au_in_pc #if using maser stars, max R2d is that of outer maser cut
        #otherwise max R2d is the cut in the central pointing

      #  if (solver==True):
      #      const1 = integrate.quad(density,min_r,max_r)
      #      const1 = const1[0]
      #  else:
      #      #if not using solver
      #      delta_r = (max_r - min_r)/numt
      #      test_r = np.array([delta_r * i for i in range(numt)])
      #      test_r += min_r
      #      if (test_r[0]==0.0):
      #          test_r[0] += 1.0e-100
      #      test_c = np.array([density(test_r[i]) for i in range(numt)])
      #      const1 = np.sum(test_c)
      #      inverse_c = np.zeros(numt)
      #      for i in range(numt):
      #          tdex = np.arange(i+1)
      #          inverse_c[i] = np.sum(test_c[tdex]) / const1


        radiusM = np.zeros(numStarsM)
        R2dM = np.zeros(numStarsM)
        randFinalM = np.zeros(numStarsM)
        ii = 0
        while(ii < numStarsM):

            rand01 = np.random.rand(1)
            rand_angle = math.acos(np.random.rand(1))

            if (solver==True):
                rand_r = -1.0
                while (rand_r < min_r):
                    rand_r = fsolve(CDF, 0.1,args=(const1,rand01)) #in pc
                    if (rand_r==0.1):
                        tmp = minimize(CDF,0.1,args=(const1,rand01))
                        rand_r = tmp.x
                        if (rand_r==0.1):
                            tmp = minimize(CDF,0.1,method='nelder-mead',args=(const1,rand01),options={'maxiter':1000})
                            rand_r = tmp.x
                            if (rand_r==0.1):
                                pdb.set_trace()
            else: 
                #if not running solver
                tmpdex = np.argmin(abs(rand01 - inverse_c))
                rand_r = test_r[tmpdex]

            rand_R2d = np.sin(rand_angle) * rand_r
            if ((rand_R2d < outerCut) & (rand_R2d > innerCut)):
                radiusM[ii] = rand_r
                R2dM[ii] = rand_R2d
                randFinalM[ii] = rand01
                print ii
                ii += 1

        R2dM *= cm_in_pc #in cm
        #read in actual stars
        #actual = np.loadtxt('/u/schappell/code/c/maser_mn_actual.dat')
        actual = np.loadtxt('/u/schappell/code/c/stars_schodel_for_mock.dat')
        pOldM_act = actual[:,1]
        randex = np.random.randint(len(pOldM_act),size=numStarsM)
        pOldM = pOldM_act[randex]
        np.savetxt('/u/schappell/code/c/maser_mn.dat',np.transpose([R2dM,pOldM]),delimiter=' ')
        radius = np.append(radius,radiusM)
        np.savetxt('/u/schappell/code/c/stars_r.dat',np.transpose([radius]),delimiter=' ')

    #pdb.set_trace()

    mnRoot = '/u/schappell/pmnOld/'+label+'_'
    mnAlpha = alpha + offset_al
    mnDelta = delta + offset_de
    mnBreak = r_break + offset_br
    #mnScaling = scaling_a + offset_a

    outCommand = './sc_mn '+mnRoot+' '

   # if (rhoModel == 'broken'):
    os.system('g++ sc_mn_maser.cpp gauss_legendre.c -o sc_mn -I/u/schappell/code/c -I/u/schappell/code/c/boost/config -gge -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest')
    outCommand += str(mnAlpha)+' '+str(mnDelta)+' '+str(mnBreak)+' '+str(max_r)+' '+str(R2d_cut)+' '+str(situation)
    outCommand += ' '+str(innerCut)+' '+str(outerCut)+' '+str(nonRadial)
    os.system(outCommand)
  #  elif (rhoModel == 'dehnen'):
   #     os.system('g++ sc_mn_fix_dehnen.cpp gauss_legendre.c -o sc_mn -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest')
    #    outCommand += str(mnScaling)+' '+str(max_r)+' '+str(R2d_cut)
     #   os.system(outCommand)



def testCMN(numStars=3705,min_r=0.0,max_r=5.0,gamma=-1.0,alpha=4.0,delta=4.0,r_break=0.5,perturb=False,R2d_cut=5.0,
            label='test_r2donly',offset_al=0.0,offset_de=0.0,offset_br=0.0,solver=True,situation=4,innerCut=5.0,outerCut=15.0,
            nonRadial=1,resume=False,cfact=0.5,numTest=0):


    #max_r in pc
    #R2d_cut in arcsec
    #r_break in pc
    #scaling_a in pc
    #offset parameters ADD that offset to the true value (alpha, delta, etc) and gives that new value to MN
    #perturb adds random perturbation (+ or -) within one sigma to all accelerations
    #rhoModel sets which model for density used, broken=borken power law (gamma, alpha, delta, and r_break)
    #         dehnen=dehnen (gamma and scaling radius)

    #situation = 1, alpha, delta, and r_break free, stars in maser mosaic not included
               # 2, only gamma free, stars in maser mosaic not included
               # 3, alpha, delta, and r_break free, stars in maser mosaic included
               # 4, only gamma free, stars in maser mosaic included

    if ((resume == False) | (numTest > 1)):
        if ((situation != 1) & (situation != 2) & (situation != 3) & (situation != 4)):
            sys.exit("Did not pick available situation, try 1, 2, 3, or 4")
        

        if (solver==False):
    #if not using solver
            numt = 100000

        R2d_cut *= dist / au_in_pc #R2d_cut now in pc
        innerCut *= dist / au_in_pc
        outerCut *= dist / au_in_pc

    #read in actual stars
        actual = np.loadtxt('/u/schappell/code/c/stars_gcows_for_mock.dat')
        are_act = actual[:,2]
        pOld_act = actual[:,3]
        
        if (nonRadial>0):
            plate_scale=0.00995 #arcsec/pixel
            gcows_zero = np.array([1500.0,1500.0])
            gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')

            xPIX = np.array([])
            yPIX = np.array([])

    #density profile
        def density(r,gamma=gamma,alpha=alpha,delta=delta,r_break=r_break):
            try:
                (r/r_break)**(-1.0*gamma)*(1.0 + (r/r_break)**delta)**((gamma-alpha)/delta)*r**2
            except:
                pdb.set_trace()
            return (r/r_break)**(-1.0*gamma)*(1.0 + (r/r_break)**delta)**((gamma-alpha)/delta)*r**2


        if (solver==True):
            const1 = integrate.quad(density,min_r,max_r)
            const1 = const1[0]
        else:
    #if not using solver
            delta_r = (max_r - min_r)/numt
            test_r = np.array([delta_r * i for i in range(numt)])
            test_r += min_r
            if (test_r[0]==0.0):
                test_r[0] += 1.0e-100
            test_c = np.array([density(test_r[i]) for i in range(numt)])
            const1 = np.sum(test_c)
            inverse_c = np.zeros(numt)
            for i in range(numt):
                tdex = np.arange(i+1)
                inverse_c[i] = np.sum(test_c[tdex]) / const1

    #pdb.set_trace()

        def beta_CDF(rprime,*needed):
            const1,yprime = needed
            value = np.array([])
            for val in rprime:
                if (val < min_r):
                    value = np.append(value,-1.0)
                else:
                    atmp = (3.0 - gamma) / delta
                    btmp = -1.0*(3.0 - alpha) / delta
                    xtmp = val**delta / (r_break**delta + val**delta)
                    tmp = scipy.special.betainc(atmp,btmp,xtmp) / scipy.special.beta(atmp,btmp)
                    value = np.append(value,abs(tmp - yprime))
            return value

        def CDF(rprime,*needed):
            const1,yprime = needed
            value = np.array([])
            for val in rprime:
                if (val < min_r):
                    value = np.append(value, -1.0)
                else:
                    tmp = integrate.quad(density,min_r,val)
                    value = np.append(value,abs(tmp[0]/const1 - yprime))
            return value

    #rand01 = np.random.rand(numStars)

        radius = np.array([])
        R2d = np.array([])
        randFinal = np.array([])
        accel = np.array([])
        are = np.array([])
        
        if (situation > 2):
            radiusM = np.array([])
            R2dM = np.array([])
            randFinalM = np.array([])

        ii = 0
        if (numTest > 1):
            numStars *= numTest
        while(ii < numStars):
        #pdb.set_trace()
        #yprime = rand01[i]

            cval = np.random.rand(1)

            rand01 = np.random.rand(1)
            rand_angle = math.acos(np.random.rand(1))

            if (solver==True):
                rand_r = -1.0
                while (rand_r < min_r):
                    rand_r = fsolve(CDF, 0.1,args=(const1,rand01)) #in pc
                    if (rand_r==0.1):
                        tmp = minimize(CDF,0.1,args=(const1,rand01))
                        rand_r = tmp.x
                        if (rand_r==0.1):
                            tmp = minimize(CDF,0.1,method='nelder-mead',args=(const1,rand01),options={'maxiter':1000})
                            rand_r = tmp.x
                            if (rand_r==0.1):
                                pdb.set_trace()
            else: 
        #if not running solver
                tmpdex = np.argmin(abs(rand01 - inverse_c))
                rand_r = test_r[tmpdex]

            rand_R2d = np.sin(rand_angle) * rand_r #in pc
            #if (cval <= cfact):
            if ((rand_R2d <= R2d_cut) & (cval < cfact)):
                if (nonRadial==0):
                    rand_r *= cm_in_pc
                    rand_R2d *= cm_in_pc
                    radius = np.append(radius,rand_r)
                    R2d = np.append(R2d,rand_R2d)
                    randex = np.random.randint(len(are_act),size=1)
                    are = np.append(are,are_act[randex])
                    accel_abs = -1.0 * GM * rand_R2d / rand_r**3 #in cm/s^2
                    accel = np.append(accel,np.random.normal(loc=accel_abs,scale=are_act[randex],size=1))
                    randFinal = np.append(randFinal,rand01)
                    print ii
                    ii += 1
                else:
                    rand_angle2 = np.random.rand(1)*2.0*pi
                    rand_x_pixels = int(round((rand_R2d * au_in_pc * math.cos(rand_angle2))/(plate_scale * dist) + gcows_zero[0]))
                    rand_y_pixels = int(round((rand_R2d * au_in_pc * math.sin(rand_angle2))/(plate_scale * dist) + gcows_zero[1]))
                #pc to arcsec, then to pixels
                    if((rand_x_pixels >= 0.0) & (rand_y_pixels >= 0.0) & (gcows[rand_y_pixels,rand_x_pixels] > 0.0)):
                        rand_r *= cm_in_pc
                        rand_R2d *= cm_in_pc
                        radius = np.append(radius,rand_r)
                        R2d = np.append(R2d,rand_R2d)
                        randex = np.random.randint(len(are_act),size=1)
                        are = np.append(are,are_act[randex])
                        accel_abs = -1.0 * GM * rand_R2d / rand_r**3 #in cm/s^2
                        accel = np.append(accel,np.random.normal(loc=accel_abs,scale=are_act[randex],size=1))
                        randFinal = np.append(randFinal,rand01)
                        xPIX = np.append(xPIX,rand_x_pixels)
                        yPIX = np.append(yPIX,rand_y_pixels)
                        print ii
                        ii += 1

            #else:
            if ((rand_R2d < outerCut) & (rand_R2d > innerCut)):
                if ((situation > 2) & (cval > cfact)):
                    radiusM = np.append(radiusM,rand_r)
                    R2dM = np.append(R2dM,rand_R2d)
                    randFinalM = np.append(randFinalM,rand01)
                    print ii
                    ii += 1
                    
    #radius *= cm_in_pc #in cm
    #R2d *= cm_in_pc #in cm

    #WRONG!!!!R2d = np.cos(np.random.rand(numStars)*pi/2.0) * radius
    #accel_mu = -1.0 * GM * R2d / radius**3 #in cm/s^2

        randex = np.random.randint(len(are_act),size=len(radius))
    #are = are_act[randex]
        pOld = pOld_act[randex]

        if (perturb == True):
            pert = np.random.rand(numStars)*2.0 - 1.0
            accel = accel + pert*are #perturb true accelerations within one sigma
    
    #pdb.set_trace()
        np.savetxt('/u/schappell/code/c/stars_mn'+label+'.dat',np.transpose([R2d,accel,are,pOld]),delimiter=' ')

        if (situation > 2):
            radiusM *= cm_in_pc
            R2dM *= cm_in_pc #in cm
        #read in actual stars
        #actual = np.loadtxt('/u/schappell/code/c/maser_mn_actual.dat')
            actual = np.loadtxt('/u/schappell/code/c/stars_schodel_for_mock.dat')
            pOldM_act = actual[:,1]
            randex = np.random.randint(len(pOldM_act),size=len(radiusM))
            pOldM = pOldM_act[randex]
            np.savetxt('/u/schappell/code/c/maser_mn'+label+'.dat',np.transpose([R2dM,pOldM]),delimiter=' ')
            radius = np.append(radius,radiusM)

    #pdb.set_trace()
        np.savetxt('/u/schappell/code/c/stars_r'+label+'.dat',np.transpose([radius]),delimiter=' ')
        np.savetxt('/u/schappell/code/c/xy_pixels'+label+'.dat',np.transpose([xPIX,yPIX]),delimiter=' ')


    if (numTest <= 1):
        mnRoot = '/u/schappell/pmnOld/'+label+'_'
        mnAlpha = alpha + offset_al
        mnDelta = delta + offset_de
        mnBreak = r_break + offset_br

        outCommand = './sc_mn '+mnRoot+' '

        pdb.set_trace()
        os.system('g++ sc_mn_mthread.cpp gauss_legendre.c -o sc_mn -std=c++11 -lpthread -I/u/schappell/code/c -I/u/schappell/code/c/boost/config -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest')
    #os.system('g++ sc_mn_maser_accel.cpp gauss_legendre.c -o sc_mn -I/u/schappell/code/c -I/u/schappell/code/c/boost/config -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest')
    
        outCommand += str(mnAlpha)+' '+str(mnDelta)+' '+str(mnBreak)+' '+str(max_r)+' '+str(R2d_cut)+' '+str(situation)
        outCommand += ' '+str(innerCut)+' '+str(outerCut)+' '+str(nonRadial)+' '+label
        os.system(outCommand)




def testModel(numStars=27,min_r=0.0,max_r=5.0,gamma=-1.0,alpha=4.0,delta=4.0,r_break=0.5,perturb=False,
              offset_al=0.0,offset_de=0.0,offset_br=0.0,numTest=100,R2d_cut=1.7,solver=False,situation=1,
              innerCut=5.0,outerCut=10.0,numStarsM=1856,nonRadial=1,label='',resume=False,cfact=0.5):

    #variance for broken power law
    #numStars is number of stars in each run, numTest is the number of iterations
    #min and max_r are the 3D radius ranges considered for the synthetic data
    #perturb adds random error within one sigma to data
    #offset_* ADDS this offset to the real model parameter and gives new value to MultiNest

    label = 'tmp'+label

    yesGamma = np.zeros(numTest)
    tot_r3d = np.array([])
    tot_r2d = np.array([])
    tot_ar = np.array([])
    totpost = np.array([])
    meanGamma = np.zeros(numTest)
    maxLGamma = np.zeros(numTest)
    mnEvid = np.zeros(numTest)

    for i in range(numTest):
        #pdb.set_trace()
        #os.system('rm /u/schappell/pmnOld/tmp_*')
        testCMN(numStars=numStars,min_r=min_r,max_r=max_r,gamma=gamma,alpha=alpha,delta=delta,r_break=r_break,
                perturb=perturb,label=label,offset_al=offset_al,offset_de=offset_de,offset_br=offset_br,
                R2d_cut=R2d_cut,solver=solver,situation=situation,innerCut=innerCut,outerCut=outerCut,
                nonRadial=nonRadial,resume=resume,cfact=cfact)
                #numStarsM=numStarsM)

        pdb.set_trace()
        tmpPost = np.loadtxt('/u/schappell/pmnOld/'+label+'_post_equal_weights.dat')
        #pdb.set_trace()
        tmphist,tmpbins = np.histogram(tmpPost[:,0]*5-3,bins=100,weights=abs(tmpPost[:,1]))
        levels = getContourLevels(tmphist)
        tmpdex = np.where(tmphist >= levels[0])[0]
        totpost = np.append(totpost,tmpPost[:,0]*5-3)
        if ((gamma >= np.min(tmpbins[tmpdex])) & (gamma <= np.max(tmpbins[tmpdex+1]))):
            yesGamma[i] = 1.0
        tmpDat = np.loadtxt('/u/schappell/code/c/stars_mn'+label+'.dat')
        #tmpOld = tmpDat[:,3]
        #tmpR2d = tmpDat[:,0]
        #tmpdex = np.where(tmpR2d <= (0.07*dist*cm_in_au))[0]
        #if (np.sum(tmpOld[tmpdex]) >= 8.48):
        #    yesStars[i] = 1.0
        tot_r2d = np.append(tot_r2d,tmpDat[:,0])
        tot_ar = np.append(tot_ar,tmpDat[:,1])
        tmpDat = np.loadtxt('/u/schappell/code/c/stars_r'+label+'.dat')
        tot_r3d = np.append(tot_r3d,tmpDat[:,0])
        mnOut = np.loadtxt('/u/schappell/pmnOld/'+label+'_summary.txt')
        meanGamma[i] = mnOut[0]*5-3
        maxLGamma[i] = mnOut[2]*5-3
        mnEvid[i] = mnOut[4]


    outFile = '/u/schappell/pmnOld/SchCor_'
    if (situation==1):
        outFile += 'free_'
    elif (situation==2):
        outFile += 'fixed_'
    elif (situation==3):
        outFile += 'maser_free_'
    elif (situation==4):
        outFile += 'maser_fixed_'
    outFile += 'g'+str(gamma)+'_a'+str(alpha)+'_d'+str(delta)+'_br'
    outFile += str(r_break)+'_nStar'+str(numStars)+'_Rcut'+str(R2d_cut)
    if (situation > 2):
        outFile += '_'+str(innerCut)+'_'+str(outerCut)
    outFile += 'arcsec'+'_meanG_maxLG_evid.dat'
    np.savetxt(outFile,np.transpose([meanGamma,maxLGamma,mnEvid]),delimiter=' ')


   # if (whatReturn=='mean'):
   #     return np.mean(mnGamma)
   # elif (whatReturn=='median'):
   #     return np.median(mnGamma)
   # elif (whatReturn=='mode'):
   #     mnGamma = np.round(mnGamma,2)
   #     tmp = stats.mode(mnGamma)
   #     return tmp[0]
   # elif (whatReturn=='pass?'):
    ghist, gbins = np.histogram(totpost,bins=100)
    r2hist, r2bins = np.histogram(tot_r2d,bins=100)
    r3hist, r3bins = np.histogram(tot_r3d,bins=100)
    arhist, arbins = np.histogram(np.log10(-1.0*tot_ar),bins=100)
    tmpout = np.zeros([numTest,8])
    tmpout[:,1] = ghist
    tmpout[:,3] = r2hist
    tmpout[:,5] = r3hist
    tmpout[:,7] = arhist
    tmpout[:,0] = np.array([(gbins[tt+1]+gbins[tt])/2.0 for tt in range(numTest)])
    tmpout[:,2] = np.array([(r2bins[tt+1]+r2bins[tt])/2.0 for tt in range(numTest)])
    tmpout[:,4] = np.array([(r3bins[tt+1]+r3bins[tt])/2.0 for tt in range(numTest)])
    tmpout[:,6] = np.array([(arbins[tt+1]+arbins[tt])/2.0 for tt in range(numTest)])
    outFile = '/u/schappell/pmnOld/SchCor_'
    if (situation==1):
        outFile += 'free_'
    elif (situation==2):
        outFile += 'fixed_'
    elif (situation==3):
        outFile += 'maser_free_'
    elif (situation==4):
        outFile += 'maser_fixed_'
    outFile += 'g'+str(gamma)+'_a'+str(alpha)+'_d'+str(delta)+'_br'
    outFile += str(r_break)+'_nStar'+str(numStars)+'_Rcut'+str(R2d_cut)
    if (situation > 2):
        outFile += '_'+str(innerCut)+'_'+str(outerCut)
    outFile += 'arcsec'+'_hist_gamma_r2d_3d_logar_bins.dat'
    np.savetxt(outFile,tmpout,delimiter=' ')
    return np.sum(yesGamma)/numTest


#def testDehnenPL(numStars=27,min_r=0.0,max_r=1.0,perturb=False,offset_a=0.0,numTest=100,whatReturn='mode',R2d_cut=1.7,solver=False):

    #tests range of dehnen power law models, range in gamma and a (scaling radius)
    #gamma is always the free parameter in MultiNest

 #   gtest = np.array([-1.0,-0.5,0.5,1.0,1.5])
 #   atest = np.array([0.7,0.8,0.9,1.0,1.1])

#    found_gamma = np.zeros((len(gtest),len(atest)))
    
 #   for gg in range(len(gtest)):
  #      for aa in range(len(atest)):
   #         tmp=testModel(numStars=numStars,min_r=min_r,max_r=max_r,gamma=gtest[gg],scaling_a=atest[aa],R2d_cut=R2d_cut,
    #                      perturb=perturb,rhoModel='dehnen',offset_a=offset_a,numTest=numTest,whatReturn=whatReturn,solver=solver)
     #       found_gamma[gg,aa] = tmp


    #pdb.set_trace()



def testBrokenPL(numStars=3705,min_r=0.0,max_r=5.0,perturb=False,offset_al=0.0,offset_de=0.0,offset_br=0.0,numTest=20,
                 R2d_cut=5.0,solver=True,situation=4,innerCut=5.0,outerCut=15.0,numStarsM=1856,nonRadial=1,label='',
                 resume=False,cfact=0.5):

    #tests range of broken power law models, range in gamma, alpha, delta, and r_break
    #gamma is always the free parameter in MultiNest

    #situation = 1, alpha, delta, and r_break free, stars in maser mosaic not included
               # 2, only gamma free, stars in maser mosaic not included
               # 3, alpha, delta, and r_break free, stars in maser mosaic included
               # 4, only gamma free, stars in maser mosaic included

    gtest = np.array([0.01])
    atest = np.array([2.5])
    dtest = np.array([3.0])
    btest = np.array([0.5])

    #if (whatReturn=='pass?'):
    percent_gamma = np.zeros((len(gtest),len(atest),len(dtest),len(btest)))
    #percent_stars = np.zeros((len(gtest),len(atest),len(dtest),len(btest)))
    #else:
     #   found_return = np.zeros((len(gtest),len(atest),len(dtest),len(btest)))
    
    for gg in range(len(gtest)):
        for aa in range(len(atest)):
            for dd in range(len(dtest)):
                for bb in range(len(btest)):
                    tmp=testModel(numStars=numStars,min_r=min_r,max_r=max_r,gamma=gtest[gg],alpha=atest[aa],
                                  delta=dtest[dd],r_break=btest[bb],perturb=perturb, offset_al=offset_al,
                                  offset_de=offset_de,offset_br=offset_br,numTest=numTest,R2d_cut=R2d_cut,
                                  solver=solver,situation=situation,innerCut=innerCut,outerCut=outerCut,
                                  numStarsM=numStarsM,nonRadial=nonRadial,label=label,resume=resume,cfact=cfact)
                    percent_gamma[gg,aa,dd,bb] = tmp
                    #percent_stars[gg,aa,dd,bb] = tmp[1]
                    #else:
                     #   found_return[gg,aa,dd,bb] = tmp


    pdb.set_trace()



def getContourLevels(probDist):
#From Breann Sitarski, corrected version of Sylvana Yelda's code
    """
    If we want to overlay countours, we need to figure out the
    appropriate levels. The algorithim is:
        1. Sort all pixels in the 2D histogram (largest to smallest)
        2. Make a cumulative distribution function
        3. Find the level at which 68% of trials are enclosed.
    """
    # Get indices for sorted pixel values (smallest to largest)
    sid0 = probDist.flatten().argsort()
    # Reverse indices, now largest to smallest
    sid = sid0[::-1]
    # Sort the actual pixel values
    pixSort = probDist.flatten()[sid]
    
    # Make a cumulative distribution function starting from the
    # highest pixel value. This way we can find the level above
    # which 68% of the trials will fall.
    cdf = np.cumsum(pixSort)
    cdf = cdf/max(cdf)
    
    # Determine point at which we reach 68% level
    percents = np.array([0.6827, .95, .997])
    levels = np.zeros(len(percents), dtype=float)
    for ii in range(len(levels)):
        # Get the index of the pixel at which the CDF
        # reaches this percentage (the first one found)
        idx = (np.where(cdf < percents[ii]))[0]
        
        # Now get the level of that pixel
        levels[ii] = pixSort[idx[-1]]
    return levels



def maserStars(maserDate='13jullgs1',innerCut=5.0,outerCut=15.0,magCut=15.5,schodel=False,lmagCut=0.0,
               onlySchodel=False,mnAlpha=4.0,mnDelta=4.0,mnBreak=0.5,max_r=5.0,situation=1,label=''):


    if (onlySchodel==True):
        schodel=True

    if (schodel==False):
        scale = 0.00995
        maserFile = np.loadtxt('/u/ghezgroup/data/gc/'+maserDate+'/combo/starfinder/tables/mag'+
                               maserDate+'_msr_kp_rms_named.lis',usecols=(1,3,4))
        mag = maserFile[:,0]
        x_as = maserFile[:,1] * scale
        y_as = maserFile[:,2] * scale

        names = np.genfromtxt('/u/ghezgroup/data/gc/'+maserDate+'/combo/starfinder/tables/mag'+
                              maserDate+'_msr_kp_rms_named.lis',usecols=(0),dtype='str')

        R2d = np.sqrt(x_as**2 + y_as**2)

        oldProb = np.ones(len(mag)) #assume old unless otherwise stated

        dbfile = '/g/ghez/data/gc/database/stars.sqlite'
    # Create a connection to the database file
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
            if ((tmpName == 'S0-38') | (tmpName == 'S0-49') | (tmpName == 'S0-35') | (tmpName == 'S1-32')):
                oldProb[i] = 1.0
                if (tmpName == 'S0-61'):
                    oldProb[i] = 0.0

        for i in range(len(names)):
            tmpName = str(names[i])
            try:
                cur.execute('SELECT name,young,old FROM stars WHERE name=?', [tmpName])
                for row in cur:
                    if ((row[1] == 'T') | (row[2] == 'F')):
                        oldProb[i] = 0.0
                    
            except:
                continue

    else:
        R2d = np.array([])
        mag = np.array([])
        oldProb = np.array([])

        ext_scale = 0.02705 #arcsec/pixel
        ext_center = [808,611]
        extMap = pyfits.getdata('/u/schappell/Downloads/AKs_fg6.fits')

        dbfile = '/g/ghez/data/gc/database/stars.sqlite'
    # Create a connection to the database file
        connection = sqlite.connect(dbfile)
    # Create a cursor object
        cur = connection.cursor()

        #cur.execute('SELECT young,old,r2d,kext FROM schodel2009')
        cur.execute('SELECT r2d,k,x,y FROM schoedel2010')
        for row in cur:
            #R2d = np.append(R2d,row[0])
            #mag = np.append(mag,row[1])
            #oldProb = np.append(oldProb,1.0)
            #if ((row[0]=='T') & (row[1]=='F')):
             #   oldProb = np.append(oldProb,0.0)
            #else:
             #   oldProb = np.append(oldProb,1.0)

            Xext = int(round(row[2]/ext_scale))+ext_center[0]
            Yext = int(round(row[3]/ext_scale))+ext_center[0]
            if ((Xext < 0) | (Yext < 0)):
                print 'Something is wrong, star is calculated as being off extinction map'
                return
            else:
                if ((Xext < 1600) & (Yext < 1600)):
                    oldProb = np.append(oldProb,1.0)
                    R2d = np.append(R2d,row[0])
                    mag = np.append(mag,row[1] + 2.7 - extMap[Xext,Yext])

    mdex = np.where((R2d > innerCut) & (R2d < outerCut) & (mag < magCut) & (oldProb > 0.0) & (mag > lmagCut))[0]
    R2d *= dist * cm_in_au
    if (onlySchodel==False):
        np.savetxt('/u/schappell/code/c/maser_mn'+label+'.dat',np.transpose([R2d[mdex],oldProb[mdex]]),delimiter=' ')
    else:
        np.savetxt('/u/schappell/code/c/onlySchodel_mn.dat',np.transpose([R2d[mdex],oldProb[mdex]]),delimiter=' ')
        mnRoot = '/u/schappell/pmnOld/Schodelonly_'+label+'_'
        
        outCommand = './sc_mn '+mnRoot+' '

        os.system('g++ sc_mn_schodel.cpp gauss_legendre.c -o sc_mn -L/u/schappell/multinest/MultiNest_v2.18_CMake/lib -lmultinest')
        outCommand += str(mnAlpha)+' '+str(mnDelta)+' '+str(mnBreak)+' '+str(max_r)+' '+str(situation)+' '+str(innerCut)
        outCommand += ' '+str(outerCut)
        os.system(outCommand)




