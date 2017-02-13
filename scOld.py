import matplotlib
matplotlib.use('PDF')
from gcwork import objects
from gcwork import starset
import sc_accel_class as acc
from gcwork import util
from gcwork import orbits
from gcwork import young
from pysqlite2 import dbapi2 as sqlite
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
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter 
import corner as cor
# Several functions were taken from:
# /ghezgroup/code/python/gcwork/polyfit/accel.py,
# velocity.py and residuals.py, written by JLu.

def usetexTrue():
    rc('text', usetex=True)
    rc('font', **{'family':'sans-serif', 'size':16})
    rc('axes', titlesize=16, labelsize=16)
    rc('xtick', labelsize=12)
    rc('ytick', labelsize=12)

def usetexFalse():
    rc('text', usetex=False)
    rc('font', family='sans-serif', size=14)
    rc('axes', titlesize=16, labelsize=16)
    rc('xtick', labelsize=14)
    rc('ytick', labelsize=14)

#root = '/u/syelda/research/gc/aligndir/'
root = '/g/ghez/align/'
home='/u/schappell/'

pi = math.pi

# Mass and Ro from S0-2 - Ghez et al. 2008
mass = 4.07e6
masse = 0.6e6
dist = 7960.0
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


def accelInfo(alnDir='14_06_18/', root_tmp='/g/ghez/align/',
              align = 'align/align_d_rms_1000_abs_t', updateErr = True,Rcut = 1.7,
              poly='polyfit_nz/fit', points='points_nz/', polyj='polyfit_nzj/fit',
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',nonRadial=0,
              starlist='all', magCut=22, nEpochs=14,ak_correct=True):

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
        idx = np.where((cnt > nEpochs) & (mag < magCut) & (r < Rcut))[0]
    else:
        in_gcows = np.zeros(len(mag))
        plate_scale = 0.00995
        gcows_zero = np.array([1500.0,1500.0])
        gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')
        gcows_dex = np.where(gcows > 0)
        np.savetxt('/u/schappell/code/c/gcows_field.dat',np.transpose([gcows_dex[1],gcows_dex[0]]),delimiter=' ')
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
        idx = np.where((cnt > nEpochs) & (mag < magCut) & (in_gcows==1) & (r < Rcut))[0]
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

    return names, r_cm, ar_cmss, are_cmss, mag



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
    origin_val = asciidata.open('/g/ghez/align/' + alnDir  + chainsDir + 'efit_summary.txt')
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



def histAccel(alnDir='14_06_18/', root_tmp='/g/ghez/align/',
              align = 'align/align_d_rms_1000_abs_t', updateErr = True,radCut = 100.0,
              poly='polyfit_nz/fit', points='points_nz/', sigma=3.0, polyj='polyfit_nzj/fit',
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains_S0-2_newRV2/',outTolatex=False,outToAL=False,
              starlist='all', plotSigAcc=False, magCut=22, nEpochs=14, verbose = True, likelihood=False):
    """
    Make a histogram of the accelerations in the radial/tangential,
    X/Y, and inline/perp. direction of motion. Also, the large outliers
    (> sigma) are printed out.

    Inputs:
    alnDir   = The root directory of an astrometry analysis
                (e.g. '08_02_16/' or './' if you are in the directory).
    align     = The align root file name (including the directory relative
                to root). Make sure that polyfit was run on this align
		output.
    poly      = The polyfit root file name (including the directory relative
                to root). This should be run on the same align as above.
    points    = The points directory.
    starlist  = Only plot specific subset of stars. Must be 'oldstars', 
                'yngstars', 'all'.
    plotSigAcc = Set to True to plot the fits for stars with significant radial
                 accelerations.
    nEpochs    = Number of epochs a star must have been detected in for this analysis 

    Output:
    plots/polyfit_hist_accel.eps (and png)
    -- Contains the histograms of the accelerations.

    plots/polyfit_hist_accel_nepochs.eps (and png)
    -- Contains a plot of number of stars vs. acceleration significance.
    """
#    outdir = root + alnDir + 'plots/'
    outdir = home + 'plots/'

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
    xchi2r = s.getArray('xchi2r')
    ychi2r = s.getArray('ychi2r')
        
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
    jxchi2 = _fit[9].tonumpy()
#    jxq = _fit[10].tonumpy()

    jy0 = _fit[11].tonumpy()
    jvy = _fit[12].tonumpy()
    jay = _fit[13].tonumpy()
    jy = _fit[14].tonumpy()
    jy0e = _fit[15].tonumpy()
    jvye = _fit[16].tonumpy()
    jaye = _fit[17].tonumpy()
    jye = _fit[18].tonumpy()
    jychi2 = _fit[19].tonumpy()
#    jyq = _fit[20].tonumpy()

#    jmag = np.zeros(len(jnames))
#    for i in range(len(jnames)):
#        try:
#            jdex=names.index(jnames[i])
#            jmag[i]=mag[jdex]
#        except:
#            jmag[i] = 0.0

    #Some of these stars have jerk fits, if they do, use those fits
#    for i in range(len(names)):
#        try:
#            jndex = np.where(jnames == names[i])
#            x0[i] = jx0[jndex]
#            y0[i] = jy0[jndex]
#            x0e[i] = jx0e[jndex]
#            y0e[i] = jy0e[jndex]
#            vx[i] = jvx[jndex]
#            vy[i] = jvy[jndex]
#            vxe[i] = jvxe[jndex]
#            vye[i] = jvye[jndex]
#            ax[i] = jax[jndex]
#            ay[i] = jay[jndex]
#            axe[i] = jaxe[jndex]
#            aye[i] = jaye[jndex]
#        except:
#            continue
    #Update errors in position


    # Make an epochs cut
    r = np.sqrt(x0**2 + y0**2)
    idx = np.where((cnt > nEpochs) & (mag < magCut) & (r < radCut))[0]
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
    xchi2r = xchi2r[idx]
    ychi2r = ychi2r[idx]
    t0x = t0x[idx]
    t0y = t0y[idx]
    names = [names[nn] for nn in idx]

    #for AG
    py.clf()
    pdb.set_trace()
    r = r[idx]
    ag1 = np.where((r < 1.7) & (cnt > 14))[0]
    ag2 = np.where((r < 1.7) & (cnt > 14) & (mag < 15.5))[0]
    co1 = plt.scatter(cnt,mag,facecolors='none',label='Entire Catalogue')
    co2 = py.scatter(cnt[ag1],mag[ag1],label='Projected Radius and Epochs Cut')
    py.scatter(mag[ag2],cnt[ag2],'ok')
    py.xlabel('Number Epochs')
    py.ylabel('K Mag')
    py.axis([0,60,22,8])
    py.legend((co1,co2),['Entire Catalogue','Proj. R + Epochs Cut'])
    py.show()
    pdb.set_trace()

#    co1 = py.scatter(cnt[ag1],mag[ag1],facecolors='none',label='Projected Radius and Epcohs Cut')
#    co2 = py.scatter(mag[ag2],cnt[ag2],'ok',label='+K Magnitude Cut')
#    py.xlabel('Number Epochs')
#    py.ylabel('K Mag')
#    py.show()

    pdb.set_trace()


    chiBins=py.arange(0.0,50,0.1)
    py.clf()
    py.hist(xchi2r,chiBins)
    py.xlabel('Reduced X Chi^2')
    py.ylabel('Number')
    py.savefig('/u/schappell/distXchi2.png')

    py.clf()
    py.hist(ychi2r,chiBins)
    py.xlabel('Reduced Y Chi^2')
    py.ylabel('Number')
    py.savefig('/u/schappell/distYchi2.png')

    ###########
    #
    # F Test -- Accel or Velocity fit?
    #
    ###########
    if f_test == True:
        # Pass in the star names and run the F test
        pass_f,xprob,yprob = run_f_test(names,pvalue,root_tmp,alnDir,align,poly,points,verbose=verbose)
        pass_fj = run_fjerk_test(names,pvalue,root_tmp,alnDir,align,poly,points,verbose=verbose)
    else:
        pass_f = np.ones(len(names)) # check everything
        pass_fj = np.ones(len(names))

    idex = np.where(pass_fj ==1)[0]

    for i in idex:
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
        xchi2r[i] = jxchi2[jdex]/(cnt[i] - 4.0)
        ychi2r[i] = jychi2[jdex]/(cnt[i] - 4.0)
        t0x[i] = jt0x[jdex]
        t0y[i] = jt0y[jdex]
#        ar[i] = jar[jdex]
#        at[i] = jat[jdex]
#        are[i] = jare[jdex]
#        ate[i] = jate[jdex]
#    sigmaR = ar/are
#    sigmaT = at/ate


    # Get accelerations in the radial/tangential direction
    r = np.hypot(x0, y0)
    ar = ((ax*x0) + (ay*y0)) / r
    at = ((ax*y0) - (ay*x0)) / r
    are =  (axe*x0/r)**2 + (aye*y0/r)**2
    are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    are =  np.sqrt(are)
    ate =  (axe*y0/r)**2 + (aye*x0/r)**2
    ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
    ate =  np.sqrt(ate)


    writeJLatex = np.where((pass_fj ==1) & (pass_f == 1) & (r < 1.0))[0]

#    jr = np.hypot(jx0, jy0)
#    jar = ((jax*jx0) + (jay*jy0)) / jr
#    jat = ((jax*jy0) - (jay*jx0)) / jr
#    jare =  (jaxe*jx0/jr)**2 + (jaye*jy0/jr)**2
#    jare += (jy0*jx0e*jat/jr**2)**2 + (jx0*jy0e*jat/jr**2)**2
#    jare =  np.sqrt(jare)
#    jate =  (jaxe*jy0/jr)**2 + (jaye*jx0/jr)**2
#    jate += (jy0*jx0e*jar/jr**2)**2 + (jx0*jy0e*jar/jr**2)**2
#    jate =  np.sqrt(jate)

#    fsp = np.where((r < 1.7) & (cnt >= 15))[0]
#    print 'For star program:'
#    print len(fsp)
#    for i in fsp:
#        print names[i]

    py.clf()
    tmpBins=py.arange(-15,15,0.5)
    py.subplot(1,2,1)
    py.hist(at/ate,tmpBins)
    py.ylim([0,300])
    py.xlabel('$\sigma_{T}$')
    py.title('Original')

    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, alnDir=alnDir, chainsDir = chainsDir)
#        jx0e,jy0e,jvxe,jvye = nzErr(jx0e, jy0e, jvxe, jvye, jt0x, jt0y, jmag, alnDir=alnDir, chainsDir = chainsDir)

    if updateErr:
        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17.1])
        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5
       # atBins=np.array([9,13.39,14.32,15.19,15.44,15.64,15.86,16.265,18])
       # deltaAt=np.array([1.455,2.227,3.4624,4.036,3.8197,4.7166,5.7052,5.7413])*1e-5
#        are = np.sqrt(are**2 + (3.1585e-5)**2)
#        ate = np.sqrt(ate**2 + (3.6118e-5)**2)
#        are = np.sqrt(are**2 + (2.0e-4)**2)
#        ate = np.sqrt(ate**2 + (1.0e-4)**2)
#        atBins=np.array([9,11,12,13,14,15,16,17,18,19,20,23])
#        deltaArr=np.array([3.5,71.0,58.0,210.0,300.0,650.0,700.0,1100.0,1900.0,2200.0,3000.0])*1e-6
       # radBins=np.array([0.0,1.2267,1.9434,2.3978,2.9501,3.3867,3.8977,4.4208,6.5])
       # deltaRa=np.array([3.6889,3.4419,2.6797,2.7718,2.6686,2.4569,1.2631,0.9888])*1e-5
        

        #Finding the delta to add in quad to the error depending on the radius and mag, whichever one, radius or mag
        #error is larger
        delta = mag*0.0
#        jdelta = jmag*0.0
#        for i in range(len(mag)):
#            for j in range(len(deltaAt)):
#                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
#                    delta[i]=deltaAt[j]
#                    for k in range(len(deltaRa)):
#                        if ((r[i] > radBins[k]) & (r[i] <= radBins[k+1])):
#                            if deltaRa[k] > delta[i]:
#                                delta[i]=deltaRa[k]
        
#        magBins = np.array([9,14.32,15.44,15.86,18])
#        radBins = np.array([0.0,1.9434,2.9501,3.8977,6.5])
#        deltaArr = np.array([[14.601,36.187,10.737,19.217],[24.742,21.776,35.729,13.194],[22.587,41.706,44.864,12.631],
#                             [49.032,46.48,50.812,28.274]])*1e-6
#        deltaArr = np.array([[11.49,2.85,11.09,28.78],[22.64,15.93,37.53,5.35],[22.59,25.97,44.86,42.86],
#                             [49.43,51.14,6.54,28.28]])*1e-6
#        for i in range(len(mag)):
#            for j in range(len(magBins)-1):
#                if ((mag[i] > magBins[j]) & (mag[i] <= magBins[j+1])):
#                    idxm = j
#            for p in range(len(radBins)-1):
#                if ((r[i] > radBins[p]) & (r[i] <= radBins[p+1])):
#                    idxr = p
#            delta[i] = deltaArr[idxm,idxr]
#            print mag[i],r[i],delta[i]

        for i in range(len(mag)):
            for j in range(len(atBins)-1):
                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
                    delta[i] = deltaArr[j]
#        for i in range(len(jmag)):
#            for j in range(len(atBins)-1):
#                if ((jmag[i] > atBins[j]) & (jmag[i] <= atBins[j+1])):
#                    jdelta[i] = deltaArr[j]                

        ate = np.sqrt(ate**2 + delta**2)
        are = np.sqrt(are**2 + delta**2)
#        jate = np.sqrt(jate**2 + jdelta**2)
#        jare = np.sqrt(jare**2 + jdelta**2)

#    tmpBins=py.arange(-10,10,0.1)
    py.subplot(1,2,2)
    py.hist(at/ate,tmpBins)
    py.xlabel('$\sigma_{T}$')
    py.title('After')
    py.ylim([0,300])
    py.savefig('/u/schappell/plots/histatan.png')


    # Lets also do parallel/perpendicular to velocity
    v = np.sqrt(vx**2 + vy**2)
    am = ((ax*vx) + (ay*vy)) / v
    an = ((ax*vy) - (ay*vx)) / v
    ame = np.sqrt((axe*vx)**2 + (aye*vy)**2) / v
    ane = np.sqrt((axe*vy)**2 + (aye*vx)**2) / v

    # Total acceleration
    atot = py.hypot(ax, ay)
    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

    # Print out some info
    foo = np.where((mag < 13.5) & (cnt > 30))[0]
    print 'For stars brighter than K=13.5 (N=%s)' % str(len(foo))
    print 'Average radial acceleration error: %5.3f +- %5.3f mas/yr^2' % \
          (are[foo].mean()*1.e3, are[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average tangential acceleration error: %5.3f +- %5.3f mas/yr^2' % \
          (ate[foo].mean()*1.e3, ate[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average X velocity error: %5.3f +- %5.3f mas/yr' % \
          (vxe[foo].mean()*1.e3, vxe[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average Y velocity error: %5.3f +- %5.3f mas/yr' % \
          (vye[foo].mean()*1.e3, vye[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average X position error: %5.3f +- %5.3f mas' % \
          (x0e[foo].mean()*1.e3, x0e[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average Y position error: %5.3f +- %5.3f mas' % \
          (y0e[foo].mean()*1.e3, y0e[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print ''

    foo = np.where((mag < 15.5) & (cnt > 30))[0]
    print 'For stars brighter than K=15.5 (N=%s)' % str(len(foo))
    print 'Average radial acceleration error: %5.3f +- %5.3f mas/yr^2' % \
          (are[foo].mean()*1.e3, are[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average tangential acceleration error: %5.3f +- %5.3f mas/yr^2' % \
          (ate[foo].mean()*1.e3, ate[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average X velocity error: %5.3f +- %5.3f mas/yr' % \
          (vxe[foo].mean()*1.e3, vxe[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average Y velocity error: %5.3f +- %5.3f mas/yr' % \
          (vye[foo].mean()*1.e3, vye[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average X position error: %5.3f +- %5.3f mas' % \
          (x0e[foo].mean()*1.e3, x0e[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)
    print 'Average Y position error: %5.3f +- %5.3f mas' % \
          (y0e[foo].mean()*1.e3, y0e[foo].std()*np.sqrt(len(foo)/(len(foo)-1.0))*1.e3)

    ##########
    #
    # KS Test... Are the distributions Normal?
    #
    ##########
    sigmaR = ar/are
    sigmaT = at/ate
 #   jsigmaR = jar/jare
 #   jsigmaT = jat/jate

    sigmaAll = concatenate([sigmaR, sigmaT])
    (ksdR, kspR) = stats.stats.kstest(sigmaR, 'norm', N=len(sigmaR))
    (ksdT, kspT) = stats.stats.kstest(sigmaT, 'norm', N=len(sigmaT))
    (ksdA, kspA) = stats.stats.kstest(sigmaAll, 'norm', N=len(sigmaAll))
    #print 'KS Test for normality (prob that observed is gaussian):'
    #print '\tRadial Acc. KS Prob. = %5.3f' % (kspR)
    #print '\tTangen Acc. KS Prob. = %5.3f' % (kspT)
    #print '\tCombo Acc. KS Prob. = %5.3f' % (kspA)
    print ''


    def makePlots(val, acc, r, v, pass_f, pos, label):
        py.figure(1)
        py.subplot(3, 2, pos)
        py.hist(val, bins=range(-25, 25, 1), color='b')
        #py.axis([-8, 8, 0, (len(val)/3) + 2])
        py.title(label)

#        fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %2d epochs  %5s  %5s'
        fmt = '%15s  %5.2f   %5.3f"  %8.3f mag/yr^2  %8.2f sigma  %2d epochs  %5s  %5s  %8.3f  %8.3f'
        hdr = '%15s  %5s  %6s  %17s  %14s  %8s  %8s %8s %7s %7s  '
        #hdr = '%15s  %5s  %6s  %17s  %14s  %17s  %14s  %8s  %8s %8s  %7s  %7s'
        #fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %8.3f mas/yr^2  %8.2f sigma  %2d epochs  %5s  %5s  %8.3f  %8.3f'
        if (label == 'Radial'):
            # Get significant radial accelerations (physical)
            py.figure(2)
            py.clf()
            py.subplots_adjust(wspace=0.0, hspace=0.2, left=0.12, right=0.88, top=0.9, bottom=0.1)
            idx = (np.where((val < -sigma) & (cnt > nEpochs)))[0]
            print '%s Significant %s (physical) Accelerations' % (len(idx), label)
            print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)', 'Nepochs', 'A/V Fit', 'J/A Fit','X PVal','Y PVal')
            if len(idx) > 0:
                for i in idx:
                    if pass_f[i] == 1:
                        fit = 'Acc'
                    elif pass_f[i] == 0:
                        fit = 'Vel'
                    if pass_fj[i] == 1:
                        jfit = 'Jrk'
                    elif pass_fj[i] == 0:
                        jfit = 'Acc'
                    print fmt  % (names[i], mag[i], r[i], acc[i]*1e3, val[i], cnt[i], fit, jfit,xprob[i],yprob[i])
                    py.subplot(1,2,1)
                    py.plot(r[i],val[i],'k.')
                    py.subplot(1,2,2)
                    py.plot(mag[i],val[i],'k.')

                    if plotSigAcc == True:
                        py.figure(2)
                        py.clf()
                        plotStar(names[i].strip(),align=align,poly=poly,points=points)

                py.subplot(1,2,1)
                py.xlabel('R (arcsec)')
                py.ylabel('a_rad (sigma)')
#                py.savefig(outdir + 'physical_radial.png')
#                py.close(2)
                
#                for i in idx:
#                    py.figure(11)
#                    py.plot(mag[i],val[i],'k.')
                rght_plt=py.subplot(1,2,2)
                py.xlabel('K')
                rght_plt.yaxis.tick_right()
                rght_plt.yaxis.set_label_position("right")
                py.ylabel('a_rad (sigma)')
                py.savefig(outdir + 'physical_radial_mag.png')
                py.close(2)
#                py.close(11)

            # Get significant unphysical accelerations (positive radial)
            idx = (np.where((val > sigma) & (cnt > nEpochs)))[0]
            print
            print '%s Significant Positive %s (unphysical) Accelerations' % \
                  (len(idx), label)
            print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)',
                         'Nepochs', 'A/V Fit', 'J/A Fit','X PVal','Y PVal')
            if len(idx) > 0:
                for i in idx:
                    if pass_f[i] == 1:
                        fit = 'Acc'
                    elif pass_f[i] == 0:
                        fit = 'Vel'
                    if pass_fj[i] == 1:
                        jfit = 'Jrk'
                    elif pass_fj[i] == 0:
                        jfit = 'Acc'
                    print fmt % (names[i], mag[i], r[i], acc[i]*1e3, val[i], cnt[i], fit, jfit,xprob[i],yprob[i])
                    py.figure(3)
                    py.plot(r[i],acc[i]*1e3,'k.')

                    if plotSigAcc == True:
                        py.figure(2)
                        py.clf()
                        plotStar(names[i].strip(),align=align,poly=poly,
                                 points=points,radial=True,suffix='_unph')

                py.xlabel('R (arcsec)')
                py.ylabel('Positive (unphysical) Radial Acceleration (mas/yr/yr)')
                py.savefig(outdir + 'unphysical_radial.png')
                print ''
                print 'Median positive radial acceleration = %5.3f mas/yr/yr' % \
                      np.median(acc[idx]*1.e3)
                    
                py.close(3)    

        if (label == 'Tangential'):
            # Get significant unphysical accelerations (tangential)
            idx = (np.where((np.abs(val) > sigma) & (cnt > nEpochs)))[0]
            print
            print '%s Significant %s (unphysical) Accelerations' % (len(idx), label)
            print hdr % ('Name','K','r (")', 'a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs', 'A/V Fit', 'J/A Fit','X PVal','Y PVal')
            if len(idx) > 0:
                for i in idx:
                    if pass_f[i] == 1:
                        fit = 'Acc'
                    elif pass_f[i] == 0:
                        fit = 'Vel'
                    if pass_fj[i] == 1:
                        jfit = 'Jrk'
                    elif pass_fj[i] == 0:
                        jfit = 'Acc'
                    print fmt % (names[i], mag[i], r[i], acc[i]*1e3, val[i], cnt[i], fit, jfit,xprob[i],yprob[i])
                    py.figure(4)
                    py.plot(r[i],acc[i]*1e3,'k.')

                    if plotSigAcc == True:
                        py.figure(2)
                        py.clf()
                        plotStar(names[i].strip(),align=align,poly=poly,points=points,
                                 radial=True,suffix='_unph')

                py.axis([0,4,-0.3,0.3])
                py.xlabel('R (arcsec)')
                py.ylabel('Tangential (unphysical) Acceleration (mas/yr/yr)')
                py.savefig(outdir + 'unphysical_tangential.png')
                print ''
                print 'Median tangential acceleration = %5.3f mas/yr/yr' % \
                      np.median(acc[idx]*1.e3)
                py.close(4)

    tag = ''
    if starlist=='all':
        tag = '_all'
    if starlist=='oldstars':
        tag = '_old'
    if starlist=='yngstars':
        tag = '_yng'

   # py.figure(3, figsize=(7.5,10))
   # py.clf()
   # makePlots(atot / atoterr, 1, 'Total')

    py.clf()

    for ii in range(4):
        py.figure(ii)
        py.figure(ii,figsize=(10,10))
        py.subplots_adjust(wspace=0.3, hspace=0.3, right=0.95, top=0.95)
        py.clf()

    # XY (in units of sigma)
    makePlots(ax / axe, ax, r, v, pass_f, 1, 'X')
    makePlots(ay / aye, ay, r, v, pass_f, 2, 'Y')
    
    # Radial/Tangential (in units of sigma)
    makePlots(ar / are, ar, r, v, pass_f, 3, 'Radial')
    makePlots(at / ate, at, r, v, pass_f, 4, 'Tangential')

    # In line of Motion and out (in units of sigma)
    makePlots(an / ane, an, r, v, pass_f, 5, 'Perpendicular to v')
    makePlots(am / ame, am, r, v, pass_f, 6, 'Parallel to v')

    #py.savefig(outdir + 'polyfit_hist_accel%s.eps' % tag)
    py.savefig(outdir + 'polyfit_hist_accel%s.png' % tag)

    # Analyze the non-central arcsec sources for Nepochs threshold
    idx = (np.where(r > 0.8))[0]
    #idx = (np.where(r > 0.0))[0]

    py.close(1)

    py.figure(5)
    py.clf()
    py.subplots_adjust(hspace=0.2, left=0.15, right=0.85,
                       top=0.9, bottom=0.1)
    py.subplot(2, 1, 1)
    py.plot(ar[idx] / are[idx], cnt[idx], 'k.')
    py.axis([-20, 20, 0, 45])
    py.xlabel('Radial Acc. Sig. (sigma)')
    py.ylabel('Number of Epochs Detected')

    py.subplot(2, 1, 2)
    py.plot(at[idx] / ate[idx], cnt[idx], 'k.')
    py.axis([-20, 20, 0, 45])
    py.xlabel('Tangent. Acc. Sig. (sigma)')
    py.ylabel('Number of Epochs Detected')

    py.savefig(outdir + 'polyfit_hist_accel_nepochs%s.png' % tag)
    py.close(5)


    ##############################################################
    #SC added
    ##############################################################
    hdr = '%15s  %5s  %6s  %17s  %14s  %17s  %14s  %8s  %8s %8s  %7s  %7s'
    fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %8.3f mas/yr^2  %8.2f sigma  %2d epochs  %5s  %5s  %8.3f  %8.3f'
#    idsig = (np.where((sigmaR < -sigma) & (np.abs(sigmaT) < sigma)))[0] #sig accel
    idnp = (np.where(((sigmaR > sigma) | (np.abs(sigmaT) > sigma)) & (cnt > 23)))[0] #non physical accel
#    writeLatex = (np.where((sigmaR < -sigma) & (np.abs(sigmaT) < sigma) & (pass_f == 1)))[0]
    idsig = (np.where(((sigmaR < -22.0) | ((cnt > 24) & (sigmaR < -sigma))) & (np.abs(sigmaT) < sigma)))[0] #sig accel
    writeLatex = (np.where((((sigmaR < -22.0) | ((cnt > 24) & (sigmaR < -sigma))) & (np.abs(sigmaT) < sigma) & (pass_f==1)) | ((pass_fj == 1) & (pass_f == 1) & (r < 1.0))))[0]
#    for i in writeLatex:
#        print names[i],sigmaR[i],sigmaT[i]


    magB = np.zeros(len(deltaArr))
    medB = np.zeros(len(deltaArr))
    medA = np.zeros(len(deltaArr))
    for i in range(len(deltaArr)):
        magB[i] = (atBins[i] + atBins[i+1])/2.0
        inhere = np.where((mag > atBins[i]) & (mag <= atBins[i+1]))[0]
        medB[i] = median(ate[inhere])
        medA[i] = median(ate[inhere])

    
    sig3 = np.delete(medA,4)
    deleteMag = np.delete(magB,4)
    deleteMag[0] = 9
    deleteMag[6] = 18
    py.clf()
#    py.plot(magB,medB*1e6,'.',label='Before')
#    py.plot(magB,deltaArr*1e6,'o',label='Delta')
#    py.plot(magB,medA*1e6,'.',label='After')
#    py.legend(loc=2)

#    errpatch = np.array([[densityBin[0]-1.0,error_density[1],error_density[2]],error_density])
    lowerErr = are*1e6
    for i in range(len(lowerErr)):
        if ((-1.0*ar[i]*1e6 - lowerErr[i]) < 1.0):
            lowerErr[i] = -1.0*ar[i]*1e6 - 1.0


    py.plot(deleteMag,sig3*1e6,label='Error')
    py.errorbar(mag[writeLatex],ar[writeLatex]*1e6*(-1),yerr=np.array([lowerErr[writeLatex],are[writeLatex]*1e6]),fmt='.',label='Radial Accel')
    py.legend(loc=2,numpoints=1)
    py.yscale('log')
    py.axis([9.0,18.0,8.,10000.])
    py.xlabel('K Mag')
    py.ylabel('Radial Acceleration ($\mu$as/yr$^2$)')
    py.savefig('/u/schappell/sig3_sample.png')

    pdb.set_trace()

#    sigmaR = ar/are
#    py.clf()
#    py.plot(cnt[writeLatex],sigmaR[writeLatex],'o')
#    py.xlabel('Epochs (#)')
#    py.ylabel('Radial Acceleration (sigma)')
#    py.savefig('/u/schappell/sigrCnt.png')

    py.clf()
    toplot = (np.where(cnt > 24))[0]
    py.errorbar(1e3*at[toplot],-1e3*ar[toplot],xerr=1e3*ate[toplot],yerr=1e3*are[toplot],fmt='.')
    py.xlabel('Tan Accel (mas/yr^2)')
    py.ylabel('Rad Accel (mas/yr^2)')
    py.savefig('/u/schappell/plots/epochCut_ar_at.png')

    py.clf()
    toplot = (np.where(((sigmaR < -22.0) | ((cnt > 24) & (sigmaR < -sigma))) & (np.abs(sigmaT) < sigma)))[0]
    py.errorbar(1e3*at[toplot],-1e3*ar[toplot],xerr=1e3*ate[toplot],yerr=1e3*are[toplot],fmt='.')
    py.xlabel('Tan Accel (mas/yr^2)')
    py.ylabel('Rad Accel (mas/yr^2)')
    py.savefig('/u/schappell/plots/epochSigCut_ar_at.png')

    py.figure(6)
    py.clf()
    py.plot(ar*1e3,at*1e3,'k.')
    leg1=py.errorbar(ar[idnp]*1e3,at[idnp]*1e3,xerr=are[idnp]*3e3,yerr=ate[idnp]*3e3,fmt='.',label='Non Phys. Accel.')
    leg2=py.errorbar(ar[idsig]*1e3,at[idsig]*1e3,xerr=are[idsig]*3e3,yerr=ate[idsig]*3e3,fmt='.',label='Sig. Accel.')
    py.plot([-10,10],[0,0],'k')
    py.plot([0,0],[-10,10],'k')
    py.axis([-1.2,4,-1.2,1.2])
    py.xlabel('Radial Acceleration (mas/yr/yr)')
    py.ylabel('Tangent Acceleration (mas/yr/yr)')
    py.legend()
    py.savefig(outdir + 'tan_rad.png')



#    py.figure(7)
#    py.clf()
#    are_bins=[np.min(are),np.max(are),(np.max(are)-np.min(are))/18.]
#    ate_bins=[np.min(ate),np.max(ate),(np.max(ate)-np.min(ate))/18.]

    print
    print '%s SC Significant physical Accelerations' % (len(idsig))
    print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)','a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs', 'A/V Fit', 'J/A Fit', 'X PVal', 'Y PVal')
    if len(idsig) > 0:
        for i in idsig:
            if pass_f[i] == 1:
                fit = 'Acc'
            elif pass_f[i] == 0:
                fit = 'Vel'
            if pass_fj[i] == 1:
                jfit = 'Jrk'
            elif pass_fj[i] == 0:
                jfit = 'Acc'
            print fmt % (names[i], mag[i], r[i], ar[i]*1e3, sigmaR[i], at[i]*1e3, sigmaT[i], cnt[i], fit, jfit, xprob[i], yprob[i])

    idx = np.where((pass_f == 1) & (sigmaR <= 0.0))[0]
    print
    print '%s Pass F Test' % (len(idx))
    print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)','a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs', 'A/V Fit', 'J/A Fit', 'X PVal', 'Y PVal')
    if len(idx) > 0:
        for i in idx:
            if pass_f[i] == 1:
                fit = 'Acc'
            elif pass_f[i] == 0:
                fit = 'Vel'
            if pass_fj[i] == 1:
                jfit = 'Jrk'
            elif pass_fj[i] == 0:
                jfit = 'Acc'
            print fmt % (names[i], mag[i], r[i], ar[i]*1e3, sigmaR[i], at[i]*1e3, sigmaT[i], cnt[i], fit, jfit, xprob[i], yprob[i])
    idx = np.where(pass_fj ==1)[0]
    hdr = '%15s  %5s  %6s  %14s  %14s  %8s  %14s  %14s '
    fmt = '%15s  %5.2f   %5.3f  %8.2f sigma  %8.2f sigma  %2d epochs  %8.2f sigma  %8.2f sigma '
    print ''
    print '%s Pass F Test for Jerk' % (len(idx))
    print hdr % ('Name','K','r (")', 'a_rad (sigma)','a_tan (sigma)','Nepochs','j_x (sigma)','j_y (sigma)')
    sig_jx = jx/jxe
    sig_jy = jy/jye

    out_jx = []
    out_jy = []
    out_jxe = []
    out_jye = []
    out_jx_sig = []
    out_jy_sig = []

    for i in writeJLatex:
        jdex = np.where(jnames == names[i])
        out_jx.append(jx[jdex])
        out_jy.append(jy[jdex])
        out_jxe.append(jxe[jdex])
        out_jye.append(jye[jdex])
        out_jx_sig.append(sig_jx[jdex])
        out_jy_sig.append(sig_jy[jdex])
    
#        x0[i] = jx0[jdex]
#        y0[i] = jx0[jdex]
#        x0e[i] = jx0e[jdex]
#        y0e[i] = jy0e[jdex]
#        vx[i] = jvx[jdex]
#        vy[i] = jvy[jdex]
#        vxe[i] = jvxe[jdex]
#        vye[i] = jvye[jdex]
#        ax[i] = jax[jdex]
#        ay[i] = jay[jdex]
#        axe[i] = jaxe[jdex]
#        aye[i] = jaye[jdex]
#        ar[i] = jar[jdex]
#        at[i] = jat[jdex]
#        are[i] = jare[jdex]
#        ate[i] = jate[jdex]
        
#        test = (GM_as_yr*(ax[i]*vy[i]/ay[i] - vx[i])/(jx[jdex] - ax[i]*jy[jdex]/ay[i]))**(1.0/3.0)
#        xbh = -test**3*ax[i]/GM_as_yr
#        ybh = -test**3*ay[i]/GM_as_yr

#        jerkt=(jx[jdex] - ax[i]*jy[jdex]/ay[i])
#        p_vxe = -GM_as_yr/(3.0*test**2*jerkt)*vxe[i]
#        p_vye = GM_as_yr/(3.0*test**2*jerkt)*(ax[i]/ay[i])*vxe[i]
#        p_axe = GM_as_yr/(3.0*test**2*ay[i]*jerkt)*(vy[i]+jy[jdex]*(ax[i]*vy[i]/ay[i] - vx[i])/(jx[jdex] - ax[i]*jy[jdex]/ay[i]))*axe[i]
#        p_aye = -GM_as_yr*ax[i]/(3.0*test**2*ay[i]**2*jerkt)*(vy[i]+jy[jdex]*(ax[i]*vy[i]/ay[i] - vx[i])/(jx[jdex] - ax[i]*jy[jdex]/ay[i]))*aye[i]
#        p_jxe = -test*jxe[jdex]/(3.0*jerkt)
#        p_jye = test*ax[i]*jye[jdex]/(3.0*ay[i]*jerkt)
        if (pass_f[i] == 1):
            print fmt % (names[i],mag[i],r[i],sigmaR[i],sigmaT[i],cnt[i],sig_jx[jdex],sig_jy[jdex])
#        print test,x0[i],y0[i],p_vxe,p_vye,p_axe,p_aye,p_jxe,p_jye
    totalSample=np.where((r < 1.7) & (cnt > 14))[0]
    print 'Total in sample (<1.7" and >14 epochs): %s' % (len(totalSample))


#    pdb.set_trace()
    #Calculating z0 aand a_z and their errors
    x0_pc = x0 * dist / au_in_pc
    y0_pc = y0 * dist / au_in_pc
    x0_cm = x0_pc * cm_in_pc
    y0_cm = y0_pc * cm_in_pc
    x0e_pc = x0e * dist / au_in_pc
    y0e_pc = y0e * dist / au_in_pc
    x0e_cm = x0e_pc * cm_in_pc
    y0e_cm = y0e_pc * cm_in_pc
#    r2d = np.sqrt(x0**2 + y0**2) # arcsec
#    r2d_pc = r2d * dist / au_in_pc
    r_cm = r * dist * cm_in_pc / au_in_pc
    re_cm = np.sqrt((x0_cm*x0e_cm/r_cm)**2 + (y0_cm*y0e_cm/r_cm)**2)
    ar_cmss = ar * asy_to_kms * 1e5 / sec_in_yr
    are_cmss = are * asy_to_kms * 1e5 / sec_in_yr
    vx_cms = vx * asy_to_kms * 1e5
    vy_cms = vy * asy_to_kms * 1e5
    vxe_cms = vxe * asy_to_kms * 1e5
    vye_cms = vye * asy_to_kms * 1e5
    vproj_cms = np.sqrt(vx_cms**2 + vy_cms**2)
    vproje_cms = np.sqrt((vx_cms*vxe_cms)**2 + (vy_cms*vye_cms)**2)/vproj_cms
    z0_cm = np.sqrt(abs(((GM * r_cm / ar_cmss)**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)

    az_cmss = GM * z0_cm / ((x0_cm**2 + y0_cm**2 + z0_cm**2)**(1.5))

    z0e_cm = (abs(G*r_cm/ar_cmss))**(4.0/3.0)*(mass_g)**(-2.0/3.0)*masse_g**2 + ((abs(GM/ar_cmss))**(2.0/3.0)*(abs(r_cm))**(-1.0/3.0)-r_cm)**2*re_cm**2
    z0e_cm += (abs(GM*r_cm))**(4.0/3.0)*(abs(ar_cmss))**(-10.0/3.0)*are_cmss**2
    z0e_cm = np.sqrt(z0e_cm)/(3.0 * z0_cm)

    r3d_cm = np.sqrt(x0_cm**2 + y0_cm**2 + z0_cm**2)
    r3de_cm = np.sqrt(((x0_cm*x0e_cm)**2 + (y0_cm*y0e_cm)**2 + (z0_cm*z0e_cm)**2)/r3d_cm**2)
    r3d_pc = r3d_cm / cm_in_pc
    r3de_pc = r3de_cm / cm_in_pc
    mext = 10.0*pi*density0*(r3d_pc**2)
    mext_smbh = mext / mass

    aze_cmss = (z0_cm*masse_g)**2 + ((3.0/2.0)*mass_g*z0_cm*r_cm*re_cm/(r3d_cm**2))**2 + (mass_g*z0e_cm*(1.0-(3.0*z0_cm**2/r3d_cm**2)))**2
    aze_cmss = np.sqrt(aze_cmss) * G / r3d_cm**3

    az_kmsyr = az_cmss * sec_in_yr / 1e5
    aze_kmsyr = aze_cmss * sec_in_yr / 1e5

    z0_pc = z0_cm / cm_in_pc
    z0e_pc = z0e_cm / cm_in_pc

    vz_cms = (GM/r3d_cm) - vproj_cms**2


    if (outTolatex == True):
        out = open(home + 'tables/totsample.tex','w')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{singlespace} \n')
        out.write('\\begin{deluxetable}{lccccccc} \n')
        out.write('\\rotate \n')
        out.write('\\tabletypesize{\\tiny} \n')
        out.write('\\setlength{\\tabcolsep}{1.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        out.write('\\begin{document} \n')
        out.write('\\tablecaption{Total Sample}\n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write('  \\colhead{$Kp$} & \n')
        out.write('  \\colhead{Epochs} & \n')
        out.write('  \\colhead{2D r} & \n')
        out.write('  \\colhead{$X_0$} & \n')
        out.write('  \\colhead{$Y_0$} & \n')
        out.write('  \\colhead{$V_x$} & \n')
        out.write('  \\colhead{$V_y$} & \n')
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(mag)} & \n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{(mas/yr)} & \n')
        out.write('  \\colhead{(mas/yr)} & \n')
        out.write('} \n')
        out.write('\\startdata \n')

        fmt = '%15s  %1s  %5.2f  %1s  %2d  %1s  %5.3f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %4s\n'
        for ii in range(len(totalSample)):

            out.write(fmt % (names[ii], '&', mag[ii], '&', cnt[ii], '&', r[ii], '&', x0[ii], '$\pm$', x0e[ii], '&', y0[ii], '$\pm$', y0e[ii], '&', vx[ii]*1e3, '$\pm$', vxe[ii]*1e3, '&', vy[ii]*1e3, '$\pm$', vye[ii]*1e3, '\\\\'))
        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\end{singlespace} \n')
        out.write('\\end{document} \n')
        out.close()


    #For BS
#    for i in range(len(names)):
#        if (names[i] == '50starg1_13'):
#            arho_zero_cmss = GM/r_cm[i]**2
#            arhoe_zero_cmss = ((G*masse_g/r_cm[i]**2)**2 + (2*GM*re_cm[i]/r_cm[i]**3)**2)**(0.5)
#            print 'For BS'
#            print arho_zero_cmss
#            print arhoe_zero_cmss
    print "Stars that are unbound:"
    for jj in writeLatex:
        if (vz_cms[jj] < 0.0):
            print names[jj]
    vz_cms = np.sqrt(abs(vz_cms))
    vze_cms = np.sqrt((G*masse_g/(2*r_cm))**2 + (GM*re_cm/(2*r_cm**2))**2 + (vproj_cms*vproje_cms)**2)/vz_cms
    if (outTolatex == True):
        jnames = [names[mm] for mm in writeJLatex]
        names = [names[nn] for nn in writeLatex]

        return names,x0[writeLatex],y0[writeLatex],mag[writeLatex],r[writeLatex],ar[writeLatex]*1e3,are[writeLatex]*1e3,sigmaR[writeLatex],z0_pc[writeLatex],z0e_pc[writeLatex],az_kmsyr[writeLatex],aze_kmsyr[writeLatex],r3d_pc[writeLatex],r3de_pc[writeLatex],mext_smbh[writeLatex], cnt[writeLatex],vz_cms[writeLatex],vze_cms[writeLatex],pass_fj[writeLatex],jnames,out_jx,out_jy,out_jxe,out_jye,out_jx_sig,out_jy_sig,xchi2r[writeJLatex],ychi2r[writeJLatex]

    elif (outToAL == True):
        return names,x0,y0,x0e,y0e,vx,vy,vxe,vye,ax,ay,axe,aye,ar,are,at,ate,t0x
    elif (likelihood == True):
        return names,r_cm,ar_cmss,are_cmss
    else:
        print "Not outputting anything"

#    return ate, are, at, ar
#    print ''
#    print np.min(x0e)
#    print np.max(x0e)
#    print ''
#    print np.min(y0e)
#    print np.max(y0e)

    #####################################################################

def run_f_test(stars, pvalue, root_tmp, alnDir, align, poly, points, returnAcc=False,
               verbose=True):
    """
    Send in list of names and check whether a velocity or accel fit is best based
    on the F test. Calls Tuan's accelclass code to do the F test.

    Input:
    stars -- List of names to run F test on
    pvalue -- Significance threshold for F test

    """
    outdir = home + 'plots/'
    
    signif = scipy.special.erfc(pvalue/np.sqrt(2.0))
    print 'Significance value: %8.3e' % signif

    import sc_accel_class as acc ##########CHANGED! root
    data = acc.accelClass(rootDir=root_tmp+alnDir,align=align,poly=poly,points=points,verbose=verbose)
#    data.computeAccel() ####CHANGED!
    data.computeFTest()
    data.computeFTestJerk()

    xFprob = data.xFProb # F test p-value for accel vs. vel
    yFprob = data.yFProb

    # Round to nearest decimal
    xFprob = np.around(xFprob,decimals=4)
    yFprob = np.around(yFprob,decimals=4)

    # Other info from this align
    names = data.names
    nEpochs = data.nEpochs # original number of epochs -- but not necessarily correct
    mag = data.mag
    x = data.x
    y = data.y
    xe = data.xerr
    ye = data.yerr
    vx = data.vx
    vy = data.vy
    vxe = data.vxe
    vye = data.vye
    ax = data.ax
    ay = data.ay
    axe = data.axe
    aye = data.aye
    r2d = data.r2d
    cnt = data.nEpochs

    acc = np.zeros(len(stars))
    accnames = []
    accXFp = []
    accYFp = []
    xprob = np.zeros(len(stars))
    yprob = np.zeros(len(stars))

    for ss in range(len(stars)):
        idx = np.where(names == stars[ss])[0]
        xFp = xFprob[idx]
        yFp = yFprob[idx]
        xprob[ss] = xFprob[idx]
        yprob[ss] = yFprob[idx]

#        print "Star, prob in X, prob in Y"
#        print stars[ss], xFp, yFp

        if ((xFp < signif) | (yFp < signif)):
            acc[ss] = 1 # Accel fit
            accnames = np.concatenate([accnames, [stars[ss]]])
            accXFp = np.concatenate([accXFp, xFp])
            accYFp = np.concatenate([accYFp, yFp])
        else:
            acc[ss] = 0 # Velocity fit
    
    # Pass back arrays specifying accelerating names and F probs
    if returnAcc == True:
        return accnames, accXFp, accYFp
    # Pass back an array specifying if accel or velocity fit is best
    else:
        return acc, xprob, yprob

#####SC acced
def run_fjerk_test(stars, pvalue, root_tmp, alnDir, align, poly, points, verbose=True, returnJrk = False):
    """
    Send in list of names and check whether a velocity or accel fit is best based
    on the F test. Calls Tuan's accelclass code to do the F test.

    Input:
    stars -- List of names to run F test on
    pvalue -- Significance threshold for F test

    """
    outdir = home + 'plots/'
    
    signif = scipy.special.erfc(pvalue/np.sqrt(2.0))
    print 'Significance value: %8.3e' % signif

#    pdb.set_trace()

    import sc_accel_class as acc ##########CHANGED!
    data = acc.accelClass(rootDir=root_tmp+alnDir,align=align,poly=poly,points=points,verbose=verbose)
#    data.computeAccel() ####CHANGED!
#    data.computeFTest()
    data.computeFTestJerk()

    # Other info from this align
    names = data.names
    nEpochs = data.nEpochs # original number of epochs -- but not necessarily correct
    mag = data.mag
    x = data.x
    y = data.y
    xe = data.xerr
    ye = data.yerr
    vx = data.vx
    vy = data.vy
    vxe = data.vxe
    vye = data.vye
    ax = data.ax
    ay = data.ay
    axe = data.axe
    aye = data.aye
    r2d = data.r2d
    cnt = data.nEpochs

    xFprobj = data.xFProbJrk
    yFprobj = data.yFProbJrk

    # Round to nearest decimal
    xFprobj = np.around(xFprobj,decimals=4)
    yFprobj = np.around(yFprobj,decimals=4)

    jerk = np.zeros(len(stars))
    jrknames = []
    jrkXFp = []
    jrkYFp = []

    for ss in range(len(stars)):
        idx = np.where(names ==stars[ss])[0]
        xFpj = xFprobj[idx]
        yFpj = yFprobj[idx]

        if ((xFpj < signif) | (yFpj < signif)):
            jerk[ss] = 1 #Jerk Fit
            jrknames = np.concatenate([jrknames, [stars[ss]]])
            jrkXFp = np.concatenate([jrkXFp, xFpj])
            jrkYFp = np.concatenate([jrkYFp, yFpj])
        else:
            jerk[ss] = 0 #Accel Fit


    # Pass back arrays specifying accelerating names and F probs
    if returnJrk == True:
        return jrknames, jrkXFp, jrkYFp
    # Pass back an array specifying if accel or velocity fit is best
    else:
        return jerk


def accelLimit(alnDir='14_06_18/', root_tmp ='/g/ghez/align/',
              align = 'align/align_d_rms_1000_abs_t', updateErr = True,
              poly='polyfit_nz/fit',polyj='polyfit_nzj/fit', points='points_nz/',
              starlist='all',epochCut=23):
    """
    Plot plane-of-the-sky acceleration limits as a function of projected distance.

    Inputs:
    alnDir   = The root directory of an astrometry analysis
                (e.g. '08_02_16/' or './' if you are in the directory).
    align     = The align root file name (including the directory relative
                to root). Make sure that polyfit was run on this align
		output.
    poly      = The polyfit root file name (including the directory relative
                to root). This should be run on the same align as above.
    points    = The points directory.
    starlist  = Only plot specific subset of stars. Must be 'oldstars', 
                'yngstars', 'all'. 
    nEpochs   = Minimum number of epochs a star must be detected in (default=20).

    Output:

    """
    outdir = home + 'plots/'

    if starlist=='all':
        tag = '_all'
        lbl = 'All'
    if starlist=='oldstars':
        tag = '_old'
        lbl = 'Old'
    if starlist=='yngstars':
        tag = '_yng'
        lbl = 'Young'

    # Load data
#    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,align=align,poly=poly,points=points,starlist=starlist)
    names,x0,y0,x0e,y0e,vx,vy,vxe,vye,ax,ay,axe,aye,ar,are,at,ate,t0=histAccel(alnDir=alnDir,root_tmp=root_tmp,align=align,
                                                                               updateErr=updateErr,poly=poly,points=points,
                                                                               polyj=polyj,starlist=starlist,outToAL=True)
#    names = s.getArray('name')
#    mag = s.getArray('mag')

#    x0 = s.getArray('x0')
#    y0 = s.getArray('y0')
#    x0e = s.getArray('x0e')
#    y0e = s.getArray('y0e')
#    vx = s.getArray('vx')
#    vy = s.getArray('vy')
#    vxe = s.getArray('vxe')
#    vye = s.getArray('vye')
#    ax = s.getArray('ax')
#    ay = s.getArray('ay')
#    axe = s.getArray('axe')
#    aye = s.getArray('aye')
#    t0 = s.getArray('t0x')

    x0_pc = x0 * dist / au_in_pc
    y0_pc = y0 * dist / au_in_pc
    x0_km = x0_pc * km_in_pc
    y0_km = y0_pc * km_in_pc
    r2d = np.sqrt(x0**2 + y0**2) # arcsec
    r2d_pc = r2d * dist / au_in_pc
    rcgs = r2d * dist * cm_in_pc / au_in_pc

    # Put velocities in km/s and accelerations in km/s/yr
    vx_kms = vx * asy_to_kms
    vy_kms = vy * asy_to_kms
    vxe_kms = vxe * asy_to_kms
    vye_kms = vye * asy_to_kms
    ax *= asy_to_kms
    ay *= asy_to_kms
    axe *= asy_to_kms
    aye *= asy_to_kms
    ar *= asy_to_kms
    at *= asy_to_kms
    are *= asy_to_kms
    ate *= asy_to_kms


    # Determine # of data points per star
    pntcnt = np.zeros(len(names))
    for ii in range(len(names)):
        pntFileName = '%s%s%s%s.points' % (root_tmp, alnDir, points, names[ii])
        pntFile = open(pntFileName)
        data = pntFile.readlines()
        pntcnt[ii] = len(data)

    # Only plot acceleration upper limits for stars detected in Nepochs:
    idx = np.where(pntcnt >= epochCut)[0]

    names = [names[ii] for ii in idx]
    x0 = x0[idx]
    y0 = y0[idx]
    x0_pc = x0_pc[idx]
    y0_pc = y0_pc[idx]
    x0e = x0e[idx]
    y0e = y0e[idx]
    vx = vx[idx]
    vy = vy[idx]
    vxe = vxe[idx]
    vye = vye[idx]
    vx_kms = vx_kms[idx]
    vy_kms = vy_kms[idx]
    vxe_kms = vxe_kms[idx]
    vye_kms = vye_kms[idx]
    ax = ax[idx]
    ay = ay[idx]
    axe = axe[idx]
    aye = aye[idx]
    ar = ar[idx]
    at = at[idx]
    are = are[idx]
    ate = ate[idx]
    r2d = r2d[idx]
    r2d_pc = r2d_pc[idx]
    rcgs = rcgs[idx]
    t0 = t0[idx]

    # Total acceleration
#    atot = py.hypot(ax, ay)
#    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot
    atot = py.hypot(ar, at)
    atoterr = np.sqrt((ar*are)**2 + (at*ate)**2) / atot

    # Construct an array of radii out to 4 arcsec in steps of 0.1''
    r = (arange(10 * 5) * 0.1) + 0.1
    r_au = r * dist
    r_pc = r_au / au_in_pc
    r_cm = r_au * cm_in_au

    # Determine the theoretical curve for a vs. r
    a_cm_s2 = -GM/ r_cm**2
    a_km_s_yr = a_cm_s2 * sec_in_yr / 1.0e5


    ########################################
    #SC added
#    r_au_s134 = 0.051 * au_in_pc
#    r_cm_s134 = r_au_s134 * cm_in_au
#    a_cm_s2_s134 = -GM / r_cm_s134**2
#    a_km_s_yr_s134 = a_cm_s2_s134 * sec_in_yr / 1.0e5
#    print 'Acceleration with z=0 at projected radius of S1-34 (0.051 pc): %5.3f' % (a_km_s_yr_s134)
 ####################################################   

    # What is the polyfit acceleration upper limit
    # along the radial component? Convert 2D cartesian to circular:
#    (ar, at, are, ate) = util.xy2circErr(x0, y0, ax, ay,
#                                         x0e, y0e, axe, aye)
    
    # What is the radial component limit (3 sigma)
    accLimRadial = ar - (3.0 * are) # km/s/yr
    alimcgs = accLimRadial * 1.0e5 / sec_in_yr
    arcgs = ar *1.0e5/sec_in_yr

    # also look at the 1 sigma value
    accLimRadial1 = ar - (1.0 * are) # km/s/yr
    alimcgs1 = accLimRadial1 * 1.0e5 / sec_in_yr

    # What is the z-lower limit as set by the accel. limits
    # for just those stars below the curve
    sig_id = []
    py.figure(1, figsize=(6,6))
    py.clf()
    _out = open(home + 'plots/'+tag+'.txt','w')
    hdr = '#%10s  %7s  %12s  %12s  %12s  %7s  %7s  %15s   %9s\n'
    fmt = '%10s  %7.2f  %12.4e  %12.4e  %12.4e  %7.2f  %7.2f  %15.2f  %11.3f\n'
    _out.write(hdr % ('Name','t_ref','X (km)', 'Y (km)', 'Zmin (km)','Vx (km/s)',\
                      'Vy (km/s)','a_rad (km/s/yr','Zmin (pc)'))

    print
    print 'Stars with (radial acc + 1 sigma) below curve:'
    print '%8s  %7s  %7s  %7s  %7s   %16s' % \
          ('Star','X (pc)', 'Y (pc)','R2D (pc)', 'Zmin (pc)', 'a_rad (km/s/yr)')
    print '-----------------------------------------------------------------'

    # Plot zmin vs. 2d radius
    # If stars fall below the curve in the above plot, then their 3D distance must
    # be greater than their projected distance.  Can therefore get z_min.
    min_cnt = 0
    for ii in range(len(rcgs)):
        # Upper limits to accelerations give lower limits to 3D radii
        # Figure out which stars are below the curve, use ar + 1 sigma to get sources (alimcgs)
        foo = (-GM * rcgs[ii] / alimcgs1[ii])**(2.0/3.0) # Gives (minimum r)-squared
        #foo = (-GM * rcgs[ii] / arcgs[ii])**(2.0/3.0) # Gives (minimum r)-squared
        if ((foo - rcgs[ii]**2) > 0): # Below the curve
            # Get the indices of these stars for later
            sig_id = concatenate([sig_id,[ii]])

            zmincgs = np.sqrt(foo - rcgs[ii]**2)
            r_p = rcgs / cm_in_pc
            zmin = zmincgs / (cm_in_pc) # pc
            zmin_km = zmincgs / 1.e5

            py.plot([r_p[ii]], [zmin], 'k.')
            py.xlabel('Projected Radius (pc)')
            py.ylabel('Minimum Z (pc)')

            print '%8s  %7.3f  %7.3f  %7.3f  %7.3f     %6.3f' % \
                  (names[ii], x0_pc[ii], y0_pc[ii], r_p[ii], zmin, ar[ii])

            # Write out the kinematic variables to a file
            _out.write(fmt % (names[ii], t0[ii],  x0_km[ii], y0_km[ii], zmin_km, \
                              vx_kms[ii], vy_kms[ii], ar[ii], zmin))

            min_cnt += 1
            
    _out.close()

    print
    print 'Number of stars below curve: %s' % min_cnt

    py.figure(1)
    py.savefig(outdir+'zmin_projR%s.png' % tag)
    py.close(1)

    # Also plot the significant radial acceleration detections,
    # with their 3 sigma error bars, for just those stars below the curve
    py.figure(2)
    py.figure(figsize=(6,6))
    py.clf()
    usetexTrue()
    dtct = []
    for ii in sig_id: # they're below the curve
        if (-ar[ii] - (3.*are[ii]) > 0.0): # and they're signif detections
            det = py.errorbar(r2d_pc[ii], -ar[ii], yerr=(3.*are[ii]), fmt='b.',ms=8)
            xt = r2d_pc[ii]
            yt = -ar[ii]
            nm = names[int(ii)]
            #text(xt,yt,nm,fontsize=9)
            dtct = np.concatenate([dtct,[ii]])
    #plot([0,0.1],[0,0],'k--')
    # Plot upper limits (from top of this program)
    for ii in range(len(r2d_pc)):
        if ii in dtct:
            continue
        else:
            uppr = py.plot([r2d_pc[ii]], [-accLimRadial[ii]], 'k_')
            py.arrow(r2d_pc[ii], -accLimRadial[ii], 0, -5, hold=True, color='black',\
                     width=0.0003, head_length=1, linewidth=0.001) #, \
            #          units='x',width=0.0005, color='black')
    py.plot(r_pc, -a_km_s_yr, '--k')
    lgd1 = 'Upper Limit'
    lgd2 = 'Detection'
    #lgd = py.legend((uppr[0], det[0]), (lgd1,lgd2), fancybox=True, numpoints=1)
    #lgdLines = lgd.get_lines()
    #lgdLines[0].set_marker('_')
    #lgdLines[1].set_marker('.')
    #lgdLines[0].set_mfc('k')
    #lgdLines[1].set_mfc('b')
    #lgdLines[0].set_ms(8)
    #lgdLines[1].set_ms(8)
    py.text(0.02, 80, r'{\bf $|$a$|_{max} = \frac{GM}{R^2}$}')
    py.xlabel(r'{\bf Projected Radius (pc)}')
    py.ylabel(r'{\bf $|$a$_R|$ (km/s/yr)}')
    py.title('Acceleration Detections \& Upper Limits')
    py.axis([0,0.1,-10,100])
    py.savefig(outdir+'r2d_accel_detectUpLim%s.png' % tag)
    py.savefig(outdir+'r2d_accel_detectUpLim%s.eps' % tag)
    py.close(2)

##############################################################################################################################
#######SC Part, plot just points I define to be sig with line of max accel
    sigmaR = ar/are
    sigmaT = at/ate
    sigma = 6.0
    idsig = (np.where((sigmaR < -sigma) & (np.abs(sigmaT) < sigma)))[0] #sig accel
    #Python having issues with this, some of the stars that are passing have names python can't deal with
#    py.figure(3)
#    py.clf()
#    sc1=py.errorbar(r2d_pc[idsig],ar[idsig],yerr=(3.*are[idsig]),fmt='b.')
#    py.plot(r_pc,a_km_s_yr,'k')
#    py.xlabel('Projected Radius (pc)')
#    py.ylabel('a$_R$ (km/s/yr)')
#    py.axis([0,0.08,-60,10])
#    for ii in idsig:
#        py.annotate(names[ii],(r2d_pc[ii],ar[ii]),fontsize='x-small')
#    for i, txt in enumerate(idsig):
#        ax.annotate(txt,(r2d_pc[idsig[i]],ar[idsig[i]]))
#    py.savefig(outdir+'sig_acc.png')
#    py.close(3)
    print '%8s  %7s  %16s' % \
          ('Star','R2D (pc)', 'a_rad (km/s/yr)')
    for ii in idsig:
        print '%8s  %7.3f    %6.3f' % \
            (names[ii], r2d_pc[ii], ar[ii])
###############################################################################################################################    

    # Find the significant accelerations above zero
    #print ''
    #for ii in sig_id:
    #    ii = int(ii)
    #    if ((ar[ii] + (3.*are[ii])) < 0.):
    #        print names[ii]

    #fig=figure(3)
    #clf()
    #ax = p3.Axes3D(fig)
    py.figure(4)
    py.clf()
    for ii in range(len(rcgs)):
        foo = (-GM * rcgs[ii] / alimcgs1[ii])**(2.0/3.0) # Gives (minimum r)-squared
        if ((foo - rcgs[ii]**2) > 0):
            if (-ar[ii] - (3.*are[ii]) > 0.0): # significant acceleration (above zero)
                zmincgs = np.sqrt(foo - rcgs[ii]**2)
                zmin = zmincgs / (cm_in_pc)

                #figure(3)
                #ax.scatter3D([x0_pc[ii]],[zmin],[y0_pc[ii]],facecolor='blue')
                #if tag == '_old':
                #    ax.plot3D([x0_pc[ii],x0_pc[ii]],[zmin,zmin],[y0_pc[ii],-0.06],\
                #              color='g',ls='--',lw='1')
                #else:
                #    ax.plot3D([x0_pc[ii],x0_pc[ii]],[zmin,zmin],[y0_pc[ii],-0.08],\
                #              color='g',ls='--',lw='1')
                rmin = np.sqrt(x0_pc[ii]**2+y0_pc[ii]**2+zmin**2)

                figure(4)
                plotVvsR(names[ii],rmin,vx_kms[ii],vy_kms[ii],vxe_kms[ii],vye_kms[ii])
                

    #figure(3)
    #ax.scatter3D([0.0],[0.0],[0.0],facecolor='orange')
    #if tag == '_old':
    #    ax.plot3D([0.0,0.0],[0.0,0.0],[0.0,-0.06],color='r',ls='--',lw='1')
    #else:
    #    ax.plot3D([0.0,0.0],[0.0,0.0],[0.0,-0.08],color='r',ls='--',lw='1')
    #ax.set_xlabel('X')
    #ax.set_ylabel('Z min')
    #ax.set_zlabel('Y')
    #savefig(outdir+'3d_xyzmin%s.png' % tag)

    py.figure(4)
    py.xlabel('3D Radius (pc)')
    py.ylabel('Proper Motion (km/s)')
    py.title('Proper Motion vs. Radius')
    py.savefig(outdir+'VvsR%s.png' % tag)

    usetexFalse()


def plotStar(starName, alnDir='13_08_21/',
             align='align/align_d_rms_1000_abs_t',
             poly='polyfit_c/fit', points='points_c/', suffix='', radial=False):

    outdir = home + 'plots/'

    s = starset.StarSet(root + alnDir + align)
    s.loadPolyfit(root + alnDir + poly, accel=0, arcsec=0)
    s.loadPolyfit(root + alnDir + poly, accel=1, arcsec=0)

    names = s.getArray('name')
    
    ii = names.index(starName)
    star = s.stars[ii]

    pointsTab = asciidata.open(root + alnDir + points + starName + '.points')

    time = pointsTab[0].tonumpy()
    x = pointsTab[1].tonumpy()
    y = pointsTab[2].tonumpy()
    xerr = pointsTab[3].tonumpy()
    yerr = pointsTab[4].tonumpy()

    fitx = star.fitXv
    fity = star.fitYv
    dt = time - fitx.t0
    fitLineX = fitx.p + (fitx.v * dt)
    fitSigX = np.sqrt( fitx.perr**2 + (dt * fitx.verr)**2 )

    fitLineY = fity.p + (fity.v * dt)
    fitSigY = np.sqrt( fity.perr**2 + (dt * fity.verr)**2 )

    if (radial == True):
        # Lets also do radial/tangential
        x0 = fitx.p
        y0 = fity.p
        vx = fitx.v
        vy = fity.v
        x0e = fitx.perr
        y0e = fity.perr
        vxe = fitx.verr
        vye = fity.verr
        
        r0 = np.sqrt(x0**2 + y0**2)

        vr = ((vx*x0) + (vy*y0)) / r0
        vt = ((vx*y0) - (vy*x0)) / r0
        vre =  (vxe*x0/r0)**2 + (vye*y0/r0)**2
        vre += (y0*x0e*vt/r0**2)**2 + (x0*y0e*vt/r0**2)**2
        vre =  np.sqrt(vre)
        vte =  (vxe*y0/r0)**2 + (vye*x0/r0)**2
        vte += (y0*x0e*vr/r0**2)**2 + (x0*y0e*vr/r0**2)**2
        vte =  np.sqrt(vte)

        r = ((x*x0) + (y*y0)) / r0
        t = ((x*y0) - (y*x0)) / r0
        rerr = (xerr*x0/r0)**2 + (yerr*y0/r0)**2
        rerr += (y0*x0e*t/r0**2)**2 + (x0*y0e*t/r0**2)**2
        rerr =  np.sqrt(rerr)
        terr =  (xerr*y0/r0)**2 + (yerr*x0/r0)**2
        terr += (y0*x0e*r/r0**2)**2 + (x0*y0e*r/r0**2)**2
        terr =  np.sqrt(terr)

        fitLineR = ((fitLineX*x0) + (fitLineY*y0)) / r0
        fitLineT = ((fitLineX*y0) - (fitLineY*x0)) / r0
        fitSigR = ((fitSigX*x0) + (fitSigY*y0)) / r0
        fitSigT = ((fitSigX*y0) - (fitSigY*x0)) / r0

        diffR = r - fitLineR
        diffT = t - fitLineT
        sigR = diffR / rerr
        sigT = diffT / terr

        idxR = np.where(abs(sigR) > 4)
        idxT = np.where(abs(sigT) > 4)
        

    diffX = x - fitLineX
    diffY = y - fitLineY
    #diff = np.hypot(diffX, diffY)
    diff = (diffX + diffY) / 2.0
    rerr = np.sqrt((diffX*xerr)**2 + (diffY*yerr)**2) / diff
    sigX = diffX / xerr
    sigY = diffY / yerr
    sig = diff / rerr

    
    # Determine if there are points that are more than 5 sigma off
    idxX = np.where(abs(sigX) > 4)
    idxY = np.where(abs(sigY) > 4)
    idx = np.where(abs(sig) > 4)

    print 'Star:        ', starName
    print '\tX Chi^2 = %5.2f (%6.2f for %2d dof)' % \
          (fitx.chi2red, fitx.chi2, fitx.chi2/fitx.chi2red)
    print '\tY Chi^2 = %5.2f (%6.2f for %2d dof)' % \
          (fity.chi2red, fity.chi2, fity.chi2/fity.chi2red)
    print 'X  Outliers: ', time[idxX]
    print 'Y  Outliers: ', time[idxY]
    if (radial):
        print 'R  Outliers: ', time[idxX]
        print 'T  Outliers: ', time[idxY]
    print 'XY Outliers: ', time[idx]
    print '*****************'
    print


    close(2)
    figure(2, figsize=(7, 8))
    clf()

    dateTicLoc = MultipleLocator(3)
    #dateTicRng = [2006, 2010]
    dateTicRng = [1995, 2010]

    maxErr = np.array([xerr, yerr]).max()
    resTicRng = [-3*maxErr, 3*maxErr]

    from matplotlib.ticker import FormatStrFormatter
    fmtX = FormatStrFormatter('%5i')
    fmtY = FormatStrFormatter('%6.2f')

    paxes = subplot(3, 2, 1)
    plot(time, fitLineX, 'b-')
    plot(time, fitLineX + fitSigX, 'b--')
    plot(time, fitLineX - fitSigX, 'b--')
    errorbar(time, x, yerr=xerr, fmt='k.')
    rng = axis()
    axis(dateTicRng + [rng[2], rng[3]])
    xlabel('Date (yrs)')
    ylabel('X (pix)')
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)
    
    paxes = subplot(3, 2, 2)
    plot(time, fitLineY, 'b-')
    plot(time, fitLineY + fitSigY, 'b--')
    plot(time, fitLineY - fitSigY, 'b--')
    errorbar(time, y, yerr=yerr, fmt='k.')
    rng = axis()
    axis(dateTicRng + [rng[2], rng[3]])
    xlabel('Date (yrs)')
    ylabel('Y (pix)')
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)
    
    paxes = subplot(3, 2, 3)
    plot(time, np.zeros(len(time)), 'b-')
    plot(time, fitSigX, 'b--')
    plot(time, -fitSigX, 'b--')
    errorbar(time, x - fitLineX, yerr=xerr, fmt='k.')
    axis(dateTicRng + resTicRng)
    xlabel('Date (yrs)')
    ylabel('X Residuals (pix)')
    paxes.get_xaxis().set_major_locator(dateTicLoc)

    paxes = subplot(3, 2, 4)
    plot(time, np.zeros(len(time)), 'b-')
    plot(time, fitSigY, 'b--')
    plot(time, -fitSigY, 'b--')
    errorbar(time, y - fitLineY, yerr=yerr, fmt='k.')
    axis(dateTicRng + resTicRng)
    xlabel('Date (yrs)')
    ylabel('Y Residuals (pix)')
    paxes.get_xaxis().set_major_locator(dateTicLoc)

    bins = np.arange(-7, 7, 1)
    subplot(3, 2, 5)
    (n, b, p) = hist(sigX, bins)
    setp(p, 'facecolor', 'k')
    axis([-5, 5, 0, 20])
    xlabel('X Residuals (sigma)')
    ylabel('Number of Epochs')

    subplot(3, 2, 6)
    (n, b, p) = hist(sigY, bins)
    axis([-5, 5, 0, 20])
    setp(p, 'facecolor', 'k')
    xlabel('Y Residuals (sigma)')
    ylabel('Number of Epochs')

    subplots_adjust(wspace=0.4, hspace=0.4, right=0.95, top=0.95)
    #savefig(outdir+'plotStar_' + starName + '.eps')
    savefig(outdir+'plotStar_' + starName + suffix + '.png')
    #show()

    ##########
    #
    # Also plot radial/tangential
    #
    ##########
    if (radial == True):
        clf()

        dateTicLoc = MultipleLocator(3)
        dateTicRng = [1995, 2008]
        
        maxErr = np.array([rerr, terr]).max()
        resTicRng = [-3*maxErr, 3*maxErr]
        
        from matplotlib.ticker import FormatStrFormatter
        fmtX = FormatStrFormatter('%5i')
        fmtY = FormatStrFormatter('%6.2f')
        
        paxes = subplot(3, 2, 1)
        plot(time, fitLineR, 'b-')
        plot(time, fitLineR + fitSigR, 'b--')
        plot(time, fitLineR - fitSigR, 'b--')
        errorbar(time, r, yerr=rerr, fmt='k.')
        rng = axis()
        axis(dateTicRng + [rng[2], rng[3]])
        xlabel('Date (yrs)')
        ylabel('R (pix)')
        paxes.xaxis.set_major_formatter(fmtX)
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        paxes.yaxis.set_major_formatter(fmtY)
        
        paxes = subplot(3, 2, 2)
        plot(time, fitLineT, 'b-')
        plot(time, fitLineT + fitSigT, 'b--')
        plot(time, fitLineT - fitSigT, 'b--')
        errorbar(time, t, yerr=terr, fmt='k.')
        rng = axis()
        axis(dateTicRng + [rng[2], rng[3]])
        xlabel('Date (yrs)')
        ylabel('T (pix)')
        paxes.xaxis.set_major_formatter(fmtX)
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        paxes.yaxis.set_major_formatter(fmtY)
        
        paxes = subplot(3, 2, 3)
        plot(time, np.zeros(len(time)), 'b-')
        plot(time, fitSigR, 'b--')
        plot(time, -fitSigR, 'b--')
        errorbar(time, r - fitLineR, yerr=rerr, fmt='k.')
        axis(dateTicRng + resTicRng)
        xlabel('Date (yrs)')
        ylabel('R Residuals (pix)')
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        
        paxes = subplot(3, 2, 4)
        plot(time, np.zeros(len(time)), 'b-')
        plot(time, fitSigT, 'b--')
        plot(time, -fitSigT, 'b--')
        errorbar(time, t - fitLineT, yerr=terr, fmt='k.')
        axis(dateTicRng + resTicRng)
        xlabel('Date (yrs)')
        ylabel('T Residuals (pix)')
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        
        bins = np.arange(-7, 7, 1)
        subplot(3, 2, 5)
        (n, b, p) = hist(sigR, bins)
        setp(p, 'facecolor', 'k')
        axis([-5, 5, 0, 20])
        xlabel('T Residuals (sigma)')
        ylabel('Number of Epochs')
        
        subplot(3, 2, 6)
        (n, b, p) = hist(sigT, bins)
        axis([-5, 5, 0, 20])
        setp(p, 'facecolor', 'k')
        xlabel('Y Residuals (sigma)')
        ylabel('Number of Epochs')
        
        subplots_adjust(wspace=0.4, hspace=0.4, right=0.95, top=0.95)
        #savefig(outdir+'plotStarRadial_' + starName + '.eps')
        savefig(outdir+'plotStarRadial_' + starName + suffix + '.png')


def fitfuncPL(p, fjac=None, e=None, cdf=None, err=None):
    """Find residuals of fit.

    For power-law, p should be list of form [A, alpha],
    while data should be list of form [e, cdf, err].
    """
    num = len(e)
    model = zeros(num, dtype=float)
    devs = zeros(num, dtype=float)

    # Set parameters
    A = p[0]
    alpha = p[1]

    model = A*(e**alpha)
    residuals = (cdf - model)/err
    status = 0

    return [status, residuals]

def fitPowerLawMP(p0=None,data=None,quiet=0):
    """Fits power law using mpfit.

    Inputs should be given as p0=[A, alpha] and
    data=[e,cdf,err]. Returns object of class
    mpfit.
    """

    print 'Initial Guess:'
    print '   A     = %6.2f' % p0[0]
    print '   alpha = %5.3f' % p0[1]

    # Set data to be passed to fit
    functargs = {'e':data[0],'cdf':data[1],'err':data[2]}

    # Set initial values and limits (no limits on parameters)
    pinfo = [{'value':0,'fixed':0,'limited':[0,0],
	      'limits':[0,0]}]*len(p0)
    for ii in range(len(p0)):
        pinfo[ii]['value'] = p0[ii]

    m = nmpfit_sy.mpfit(fitfuncPL, p0, functkw=functargs, parinfo=pinfo,
		    quiet=quiet)
    if (m.status <= 0):
        print 'Error message = ', m.errmsg

    p = m.params                # Best-fit parameters
    perr = m.perror             # Error in parameter fits
                                # from covariance matrix

    m.dof = len(data[0])-len(p) # Number of degrees of freedom
    Rchi2 = m.fnorm/m.dof       # Reduced Chi^2 statistic

    p = m.params          # Best-fit parameters
    perr = m.perror       # Error in parameter fits from covariance matrix
    #x   = (data[0])[nonzero(data[2])]  
    #m.dof = len(x)-len(p) # Number of degrees of freedom

    print 'Final Solution:'
    print '   A      = %6.2f +/- %5.2f' % (p[0],perr[0])
    print '   alpha  = %5.3f +/- %5.3f' % (p[1],perr[1])
    print '   chi^2  = %7.4f' % m.fnorm
    print '   Rchi^2 = %7.4f' % Rchi2

    eline = arange(.01,1.0,.1)
    cdfline = p[0]*eline**p[1]

    return m


################################################################################################
########SC added, taken from syYoung
def plotXandY(starName,alnDir='13_08_21/',align='align/align_d_rms_1000_abs_t',
             poly='polyfit_c/fit', points='points_c/'):

    outdir = home + 'plots/'
    s = starset.StarSet(root + alnDir + align)
    names = s.getArray('name')
    ii = names.index(starName)
    star = s.stars[ii]
    pointsTab = asciidata.open(root + alnDir + points + starName + '.points')

    time = pointsTab[0].tonumpy()
    x = pointsTab[1].tonumpy()
    y = pointsTab[2].tonumpy()
    xerr = pointsTab[3].tonumpy()
    yerr = pointsTab[4].tonumpy()

#Velocity residuals
    fitx = star.fitXv
    fity = star.fitYv
    dt = time - fitx.t0
    fitLineX = fitx.p + (fitx.v * dt)
    fitSigX = np.sqrt( fitx.perr**2 + (dt * fitx.verr)**2 )
    fitLineY = fity.p + (fity.v * dt)
    fitSigY = np.sqrt( fity.perr**2 + (dt * fity.verr)**2 )

    diffX = x - fitLineX
    diffY = y - fitLineY
    #diff = np.hypot(diffX, diffY)
    diff = (diffX + diffY) / 2.0
    rerr = np.sqrt((diffX*xerr)**2 + (diffY*yerr)**2) / diff
    sigX = diffX / xerr
    sigY = diffY / yerr
    sig = diff / rerr

    
    # Determine if there are points that are more than 5 sigma off
    idxX = np.where(abs(sigX) > 4)
    idxY = np.where(abs(sigY) > 4)
    idx = np.where(abs(sig) > 4)

    print 'Star:        ', starName
    print '\tX Chi^2 = %5.2f (%6.2f for %2d dof)' % \
          (fitx.chi2red, fitx.chi2, fitx.chi2/fitx.chi2red)
    print '\tY Chi^2 = %5.2f (%6.2f for %2d dof)' % \
          (fity.chi2red, fity.chi2, fity.chi2/fity.chi2red)
    print 'X  Outliers: ', time[idxX]
    print 'Y  Outliers: ', time[idxY]
    if (radial):
        print 'R  Outliers: ', time[idxX]
        print 'T  Outliers: ', time[idxY]
    print 'XY Outliers: ', time[idx]
    print '*****************'
    print

    close(2)
    figure(2, figsize=(7, 8))
    clf()

    dateTicLoc = MultipleLocator(3)
    dateTicRng = [2006, 2011]
    #dateTicRng = [1995, 2011]

    maxErr = np.array([xerr, yerr]).max()
    resTicRng = [-3*maxErr, 3*maxErr]

    from matplotlib.ticker import FormatStrFormatter
    fmtX = FormatStrFormatter('%5i')
    fmtY = FormatStrFormatter('%6.2f')

    paxes = subplot(2, 1, 1)
    plot(time, fitLineX, 'b-')
    plot(time, fitLineX + fitSigX, 'b--')
    plot(time, fitLineX - fitSigX, 'b--')
    errorbar(time, x, yerr=xerr, fmt='k.')
    rng = axis()
    axis(dateTicRng + [rng[2], rng[3]])
    xlabel('Date (yrs)')
    ylabel('X (arcsec)')
#    title('X chi2_red = %4.2f' % fitx.chi2red)
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)
    
    paxes = subplot(2, 1, 2)
    plot(time, fitLineY, 'b-')
    plot(time, fitLineY + fitSigY, 'b--')
    plot(time, fitLineY - fitSigY, 'b--')
    errorbar(time, y, yerr=yerr, fmt='k.')
    rng = axis()
    axis(dateTicRng + [rng[2], rng[3]])
    xlabel('Date (yrs)')
    ylabel('Y (arcsec)')
#    title('Y chi2_red = %4.2f' % fity.chi2red)
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)

    subplots_adjust(wspace=0.4, hspace=0.4, right=0.95, top=0.95)
    savefig(outdir+'plotXandY_' + starName + '.png')


#Acceleration residuals
    fitxa = star.fitXa
    fitya = star.fitYa
    dta = time - fitxa.t0
    fitLineXa = fitxa.p + (fitxa.v * dta)
    fitSigXa = np.sqrt( fitxa.perr**2 + (dta * fitxa.verr)**2 )
    fitLineYa = fitya.p + (fitya.v * dta)
    fitSigYa = np.sqrt( fitya.perr**2 + (dta * fitya.verr)**2 )

    diffXa = x - fitLineXa
    diffYa = y - fitLineYa
    #diff = np.hypot(diffX, diffY)
    diffa = (diffXa + diffYa) / 2.0
    rerra = np.sqrt((diffXa*xerr)**2 + (diffYa*yerr)**2) / diff
    sigXa = diffXa / xerr
    sigYa = diffYa / yerr
    siga = diffa / rerr




########################################################################
####SC added######################
def starInfo(starName, alnDir='14_06_18/',
              align = 'align/align_d_rms_1000_abs_t',
              poly='polyfit_nz/fit', points='points_nz/',
              starlist='all', nEpochs=10):
    """
    Finds info about given star

    Inputs:
    starName  = name of the star
    alnDir   = The root directory of an astrometry analysis
                (e.g. '08_02_16/' or './' if you are in the directory).
    align     = The align root file name (including the directory relative
                to root). Make sure that polyfit was run on this align
		output.
    poly      = The polyfit root file name (including the directory relative
                to root). This should be run on the same align as above.
    points    = The points directory.
    starlist  = Only plot specific subset of stars. Must be 'oldstars', 
                'yngstars', 'all'. 
    nEpochs   = Minimum number of epochs a star must be detected in (default=20).

    Output:

    """
    outdir = home + 'plots/'

    if starlist=='all':
        tag = '_all'
        lbl = 'All'
    if starlist=='oldstars':
        tag = '_old'
        lbl = 'Old'
    if starlist=='yngstars':
        tag = '_yng'
        lbl = 'Young'

    # Load data
    s = loadPop(alnDir=alnDir,align=align,poly=poly,points=points,starlist=starlist)

    names = s.getArray('name')
    mag = s.getArray('mag')
    ii = names.index(starName)
    star = s.stars[ii]

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
    t0 = s.getArray('t0x')
    cnt = s.getArray('cnt')

    x0_pc = x0 * dist / au_in_pc
    y0_pc = y0 * dist / au_in_pc
    x0_km = x0_pc * km_in_pc
    y0_km = y0_pc * km_in_pc
    r2d = np.sqrt(x0**2 + y0**2) # arcsec
    r2d_pc = r2d * dist / au_in_pc
    rcgs = r2d * dist * cm_in_pc / au_in_pc

    # Put velocities in km/s and accelerations in km/s/yr
    vx_kms = vx * asy_to_kms
    vy_kms = vy * asy_to_kms
    vxe_kms = vxe * asy_to_kms
    vye_kms = vye * asy_to_kms
    ax *= asy_to_kms
    ay *= asy_to_kms
    axe *= asy_to_kms
    aye *= asy_to_kms

    (ar, at, are, ate) = util.xy2circErr(x0, y0, ax, ay, x0e, y0e, axe, aye)

    print 'Projected radius of %15s:  %5.3f pc' % (starName, r2d_pc[ii])
    print 'Radial acceleration of %15s:  %5.3f +- %5.3f km/s/yr' % (starName, ar[ii], 3.*are[ii])
    print 'Tangent acceleration of %15s:  %5.3f +- %5.3f km/s/yr' % (starName, at[ii], 3.*ate[ii])

    return x0[ii],y0[ii],mag[ii],cnt[ii]


def scPopCut(alnDir='11_10_26/',
              align = 'align/align_d_rms_1000_abs_t',
              poly='polyfit_c/fit', points='points_c/',
              nEpochs=30, magCut=16):

    outdir = home + 'plots/'
    alllist=loadPop(starlist='all')
    oldlist=loadPop(starlist='oldstars')

    amag=alllist.getArray('mag')
    omag=oldlist.getArray('mag')
    acnt=alllist.getArray('cnt')
    ocnt=oldlist.getArray('cnt')

    aid_e=np.where(acnt > nEpochs)[0]
    oid_e=np.where(ocnt > nEpochs)[0]
    aid_m=np.where(amag < magCut)[0]
    oid_m=np.where(omag < magCut)[0]
    aid=np.where((acnt > nEpochs) & (amag < magCut))[0]
    oid=np.where((ocnt > nEpochs) & (omag < magCut))[0]

    print '%d of all stars above epoch cut' % (len(aid_e))
    print ' '
    print '%d of old stars above epoch cut' % (len(oid_e))
    print ' '
    print '%d of all stars below mag cut' % (len(aid_m))
    print ' '
    print '%d of old stars below mag cut' % (len(oid_m))
    print ' '
    print '%d of all stars with both cuts' % (len(aid))
    print ' '
    print '%d of old stars with both cuts' % (len(oid))
    print ' '

    apl=py.plot(acnt,amag,'g.')
    opl=py.plot(ocnt,omag,'b.')
    py.xlabel('# Epochs')
    py.ylabel('K mag')
    py.legend(('Young','Not Young'))
#    py.show()
    py.savefig(outdir + 'all_magvepoch.png')


#    paxes=subplot(2,2,1)
#    errorbar(time, diffXv,yerr=sigXv,fmt='b.')

def plotpos_mag(alnDir='13_08_21/', sigma=5,
                align='align/align_d_rms_1000_abs_t',
                poly='polyfit_c/fit', points='points_c/',nEpochs=30):

    #Bring in image of central 10", taken from sythesis
    image_file='/u/ghezgroup/data/gc/14maylgs2/combo/mag14maylgs2_kp.fits'
    img_scale=0.00993
    sgra=[624.5,726.3]
    img=pyfits.getdata(image_file)
    imgsize=(img.shape)[0]
    pixL = np.arange(0,imgsize)
    xL = [-1*(xpos - sgra[0])*img_scale for xpos in pixL]
    yL = [(ypos - sgra[1])*img_scale for ypos in pixL]

    outdir = home +'plots/'
    f=open('output.csv')
    sq_known=f.read().splitlines() #reading in stars that are known to be old from sqlite table
    old_names=[]
    for i in range(len(sq_known)):
        old_names.append(sq_known[i][1:-1])

    not_young=loadPop(alnDir=alnDir,starlist='oldstars',align=align,poly=poly,points=points)
    ny_names=not_young.getArray('name')
    id_old=[]
    for nm in old_names:
        try:
            idx = ny_names.index(nm)
            id_old.append(idx)
        except ValueError:
            continue

    id_unk=np.array([i for i in range(len(ny_names)) if i not in id_old])

    ny_cnt=not_young.getArray('cnt')
    ny_mag=not_young.getArray('mag')
    ny_x=not_young.getArray('x0')
    ny_y=not_young.getArray('y0')
    ny_xe = not_young.getArray('x0e')
    ny_ye = not_young.getArray('y0e')
    ny_ax = not_young.getArray('ax')
    ny_ay = not_young.getArray('ay')
    ny_axe = not_young.getArray('axe')
    ny_aye = not_young.getArray('aye')

    ny_r = np.hypot(ny_x, ny_y)
    ny_ar = ((ny_ax*ny_x) + (ny_ay*ny_y)) / ny_r
    ny_at = ((ny_ax*ny_y) - (ny_ay*ny_x)) / ny_r
    ny_are =  (ny_axe*ny_x/ny_r)**2 + (ny_aye*ny_y/ny_r)**2
    ny_are += (ny_y*ny_xe*ny_at/ny_r**2)**2 + (ny_x*ny_ye*ny_at/ny_r**2)**2
    ny_are =  np.sqrt(ny_are)
    ny_ate =  (ny_axe*ny_y/ny_r)**2 + (ny_aye*ny_x/ny_r)**2
    ny_ate += (ny_y*ny_xe*ny_ar/ny_r**2)**2 + (ny_x*ny_ye*ny_ar/ny_r**2)**2
    ny_ate =  np.sqrt(ny_ate)

    ny_sigmaR = ny_ar/ny_are
    ny_sigmaT = ny_at/ny_ate
    ny_scale=6*(17-ny_mag)**2
#mag, x, y, sigmaT, sigmaR, cnt, scale

    omag=ny_mag[id_old]
    ox=ny_x[id_old]
    oy=ny_y[id_old]
    osr=ny_sigmaR[id_old]
    ost=ny_sigmaT[id_old]
    ocn=ny_cnt[id_old]
    osc=ny_scale[id_old]
    onm = [ny_names[nn] for nn in id_old]
#    print onm

    umag=ny_mag[id_unk]
    ux=ny_x[id_unk]
    uy=ny_y[id_unk]
    usr=ny_sigmaR[id_unk]
    ust=ny_sigmaT[id_unk]
    ucn=ny_cnt[id_unk]
    usc=ny_scale[id_unk]
    unm=[ny_names[nn] for nn in id_unk]

    oa_id = (np.where((osr < -sigma) & (np.abs(ost) < sigma) & (ocn > nEpochs)))[0]
    on_id = (np.where(((osr >= -sigma) | (np.abs(ost) >= sigma)) & (ocn > nEpochs)))[0]
    ua_id = (np.where((usr < -sigma) & (np.abs(ust) < sigma) & (ucn > nEpochs)))[0]
    un_id = (np.where(((usr >= -sigma) | (np.abs(ust) >= sigma)) & (ucn > nEpochs)))[0]

    py.clf()
    py.figure()
    py.imshow(np.log10(img+1), aspect='equal', interpolation='bicubic',
              extent=[max(xL), min(xL), min(yL), max(yL)],vmin=2.2,vmax=5,
              origin='lowerleft', cmap=py.cm.gray_r)
    cmap = py.cm.spectral
    norm = py.normalize(0.0, 1.0)
    py.plot([0],[0],'k+',ms=7,mew=2)
    na_u_pl=py.scatter(ux[un_id],uy[un_id],s=usc[un_id],facecolors='none',edgecolors='g')
    a_u_pl=py.scatter(ux[ua_id],uy[ua_id],s=usc[ua_id],facecolors='g',edgecolors='g')
    na_o_pl=py.scatter(ox[on_id],oy[on_id],s=osc[on_id],facecolors='none',edgecolors='r')
    a_o_pl=py.scatter(ox[oa_id],oy[oa_id],s=osc[oa_id],facecolors='r',edgecolors='r')
    py.xlabel('X (arcsec)')
    py.ylabel('Y (arcsec)')
    py.axis([5,-4,-4,4])
#    py.legend((na_o_pl,a_o_pl,na_u_pl,a_u_pl),('Old','Old Accelerating','Unknown','Unknown Accelerating'),scatterpoints=1)
#    py.savefig(outdir + 'test.png')
    py.savefig(outdir + 'unk_old_a_na.png')



#    Star   X (pc)   Y (pc)  Zmin (pc)    a_rad (km/s/yr)
#--------------------------------------------------------
#    S0-3    0.013    0.005    0.013     -34.670
#    S0-5    0.007   -0.014    0.007     -53.029
#  * S0-26    0.013    0.008    0.026     -8.939
#  *  S0-8   -0.009    0.006    0.007     -88.541
#    S0-4    0.016   -0.012    0.012     -27.715
#    S0-7    0.019    0.004    0.029     -7.195
#    S0-9    0.008   -0.023    0.093     -0.111
#   S0-30   -0.019   -0.015    0.024     -9.248
#  * S0-31    0.021    0.017    0.070     -0.167
#  * S0-14   -0.030   -0.011    0.051     -2.275
#  *  S1-3    0.014    0.034    0.037     -4.243
#  *  S1-1    0.040    0.001    0.056     -1.862
#  *  S1-2    0.002   -0.039    0.056     -1.823
#    S1-4    0.033   -0.026    0.023     -6.586
#  *  S1-8   -0.023   -0.035    0.069     -1.062
#  * irs16C    0.042    0.021    0.045     -2.726
#  * S1-33   -0.048   -0.001    0.041     -2.623
#  * S1-12   -0.030   -0.040    0.085     -0.570
#  * S1-14   -0.052   -0.014    0.025     -4.304
# irs16SW    0.042   -0.037    0.034     -3.211
#  * S1-21   -0.064    0.004    0.015     -3.668
#  * S1-22   -0.062   -0.020    0.019     -3.416
#    S2-4    0.058   -0.057    0.049     -1.338
#    S2-6    0.064   -0.052    0.009     -2.243
#  * irs29N   -0.061    0.054    0.022     -1.840
#   S2-17    0.051   -0.073    0.098     -0.310




def maserRes(alnDir='13_08_21/',
              align = 'align/align_d_rms_1000_abs_t', nEpochs=22,
              poly='polyfit_c/fit', points='points_c/',starlist='all'):
    """
    Calculates the RMS between given directory and that where one maser is
    taken out of setting zero point.

    Inputs:
    alnDir = root directory of the "absolute" analysis where all 7 masers are
             considered
    alig = align root file name, make sure polyfit was run on this align
    poly = the polyfit root file name for the "absolute" rest frame, or when
           all 7 masers are considered
    points = correspondign points directory
    starlist = only calculates the RMS between different reference frames for
               this set of star, must be 'oldstars,' 'yngstars,' or 'all'
    """

    outdir = home + 'plots/'


    name_tag=np.array(['all7','no10ee','no12n','no15ne','no17','no28','no7','no9'])
    root_tmp=root
    to_print=[]
    for i in range(8):
        if i > 0:
            root_tmp='/u/syelda/research/gc/aligndir/drop_maser_tests/'
            alnDir='13_10_16_'+name_tag[i]+'/'
    #NOTE: changed the loadPop such that I could change the root directory
        s=loadPop(alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points,root_tmp=root_tmp)
        names = s.getArray('name')

    # In arcsec, mas/yr, mas/yr^2
        x0 = s.getArray('x0')
        y0 = s.getArray('y0')
        x0e = s.getArray('x0e')
        y0e = s.getArray('y0e')
#        vx = s.getArray('vx')
#        vy = s.getArray('vy')
#        vxe = s.getArray('vxe')
#        vye = s.getArray('vye')
        ax = s.getArray('ax')
        ay = s.getArray('ay')
        axe = s.getArray('axe')
        aye = s.getArray('aye')
        cnt = s.getArray('cnt')
        mag = s.getArray('mag')

    # Get accelerations in the radial/tangential direction
        r = np.hypot(x0, y0)
        ar = ((ax*x0) + (ay*y0)) / r
        at = ((ax*y0) - (ay*x0)) / r
        are =  (axe*x0/r)**2 + (aye*y0/r)**2
        are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
        are =  np.sqrt(are)
        ate =  (axe*y0/r)**2 + (aye*x0/r)**2
        ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
        ate =  np.sqrt(ate)
        sigmaT=at/ate
        sigmaR=ar/are
        
#        to_print.append(names[2240])

        #Make an array for given maser directory
        if i == 0:
            all7=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_all7=names
        elif i==1:
            no10ee=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no10ee=names
        elif i==2:
            no12n=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no12n=names
        elif i==3:
            no15ne=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no15ne=names
        elif i==4:
            no17=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no17=names
        elif i==5:
            no28=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no28=names
        elif i==6:
            no7=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no7=names
        elif i==7:
            no9=np.array([x0,x0e,y0,y0e,r,cnt,ar,are,at,ate])
            n_no9=names

#    print to_print
    index=[]
    names=[]
    rms_at=[]
    rms_ar=[]
    radius=[]
    at_err=[]
    ar_err=[]
    
    #Stars aren't matched between directories, find the stars and gather info for stars that are in all directories
    for obj in n_all7:
        iabs=n_all7.index(obj)
        print iabs
        print all7[5,iabs]
        if int(all7[5,iabs]) > nEpochs:
            try:
                ione=n_no10ee.index(obj)
                itwo=n_no12n.index(obj)
                ithree=n_no15ne.index(obj)
                ifour=n_no17.index(obj)
                ifive=n_no28.index(obj)
                isix=n_no7.index(obj)
                iseven=n_no9.index(obj)
                names.append(obj)
                index.append(iabs)
                radius.append(all7[4,iabs])
                tarray1=np.array([no10ee[6,ione],no12n[6,itwo],no15ne[6,ithree],no17[6,ifour],no28[6,ifive],
                                  no7[6,isix],no9[6,iseven]])-all7[6,iabs]
                tarray2=np.array([no10ee[8,ione],no12n[8,itwo],no15ne[8,ithree],no17[8,ifour],no28[8,ifive],
                                  no7[8,isix],no9[8,iseven]])-all7[8,iabs]
                tmp1=np.sqrt(sum(tarray1**2)/7.0)*1e3
                tmp2=np.sqrt(sum(tarray2**2)/7.0)*1e3
#                tarray1=np.array([all7[6,iabs],no10ee[6,ione],no12n[6,itwo],no15ne[6,ithree],no17[6,ifour],no28[6,ifive],
#                                  no7[6,isix],no9[6,iseven]])
#                tarray2=np.array([all7[8,iabs],no10ee[8,ione],no12n[8,itwo],no15ne[8,ithree],no17[8,ifour],no28[8,ifive],
#                                  no7[8,isix],no9[8,iseven]])
#                tmp1=np.std(tarray1)
#                tmp2=np.std(tarray2)
                rms_ar.append(tmp1)
                rms_at.append(tmp2)
                at_err.append(tmp2/(all7[9,iabs]*1e6))
                ar_err.append(tmp1/(all7[7,iabs]*1e6))
            except ValueError:
                continue
            

    print names
    py.figure(1)
    py.clf()
    pl_rad=py.plot(radius,rms_ar,'.',label='RMS radial accel')
    pl_tan=py.plot(radius,rms_at,'.',label='RMS tangen accel')
    py.xlabel('Radius (")')
    py.ylabel('RMS (micro-as/yr)')
    py.legend()
    py.savefig(outdir+'rms_rad_tan.png')

    py.figure(2)
    py.clf()
    npl_rad=py.plot(radius,ar_err,'.',label='Radial Accel')
    npl_tan=py.plot(radius,at_err,'.',label='Tangen Accel')
    py.xlabel('Radius (")')
    py.ylabel('RMS (accel error)')
    py.legend()
    py.savefig(outdir+'rms_radtan_err.png')




def calcAccel(s,names):

    """
    Computes accelerations using the origin from S0-2 orbit fit. Takes in x and y
    values for all stars given in s, corrects the x and y for new zero, reruns
    acceleration fits and returns values

    Inputs:
    s      = name of loadPop
    names  = array of names in s

    """

    #Values for zero from S0-2 orbit, taken from efit, efit in camera coordinates, x0 and vx already corrected
    ori_x0= 0.398237949233159853E-02
    ori_y0= -0.852560358066134505E-02
    ori_vx= -0.184690930460452768E-03
    ori_vy= 0.548734449435500692E-03
    ori_t0= 0.200235378468537306E+04
    ori_xerr= 0.632109315557945173E-03
    ori_yerr= 0.136635578849187800E-02
    ori_vxerr= 0.378284674963935368E-04
    ori_vyerr= 0.378284674963935368E-04
    ori_terr= 0.771053652701554799E-02

    for i in range(len(names)):
        ind=names.index(names[i])

        x=s.stars[ind].getArrayAllEpochs('x')
        y=s.stars[ind].getArrayAllEpochs('y')

        xerr_a = s.stars[ind].getArrayAllEpochs('xerr_a')
        yerr_a = s.stars[ind].getArrayAllEpochs('yerr_a')
        xerr_p = s.stars[ind].getArrayAllEpochs('xerr_p')
        yerr_p = s.stars[ind].getArrayAllEpochs('yerr_p')

#        xerr = np.sqrt(xerr_p**2 + xerr_a**2)
#        yerr = np.sqrt(yerr_p**2 + yerr_a**2)

        t0 = s.stars[ind].fitXa.t0    # use the same T0 as polyfit
        years = np.array(s.stars[ind].years)

        #Update X and Y positions with new origin
        x -= ((years-ori_t0)*ori_vx + ori_x0)
        y -= ((years-ori_t0)*ori_vy - ori_y0)
        
        #Update X and Y errors
        xerr = np.sqrt(xerr_p**2 + xerr_a**2 + ori_xerr**2 + ((years-ori_t0)*ori_vxerr)**2 + (ori_vx*ori_terr)**2)
        yerr = np.sqrt(yerr_p**2 + yerr_a**2 + ori_yerr**2 + ((years-ori_t0)*ori_vyerr)**2 + (ori_vy*ori_terr)**2)

        #Run accel fit
        xfit=fitAccel(years - t0, x, xerr)
        yfit=fitAccel(years - t0, y, yerr)

        #Update values
        s.x[i] = xfit.params[0]
        s.xe[i] = xfit.perror[0]
        s.ax[i] = xfit.params[2]
        s.axe[i] = xfit.perror[2]
        s.vx[i] = xfit.params[1]
        s.vxe[i] = xfit.perror[1]

        s.y[i] = yfit.params[0]
        s.ye[i] = yfit.perror[0]
        s.ay[i] = yfit.params[2]
        s.aye[i] = yfit.perror[2]
        s.vy[i] = yfit.params[1]
        s.vye[i] = yfit.perror[1]


def fitAccel(x, y,  yerr):
     """ takes in an array of x and y position as well as time and
     fit for the acceleration. Will return an mpfit object
     """
     guessY =  [y[0],(np.amax(y) - np.amin(y))/(np.amax(x)-np.amin(x)), 0]
     functargsY = {'x':x, 'y':y, 'err': yerr}
     yfitLin = nmpfit_sy.mpfit(fitfunPoly2, guessY, functkw=functargsY,quiet=1)

     return yfitLin



def fitfunPoly2(p, fjac=None, x=None, y=None, err=None):
    # second order polynomial
    fun = p[0] + p[1]*x + 0.5*p[2]*x**2
        
    # deviations from the model
    deviates = (y - fun)/err
    return [0, deviates]



def getSigmas(alnDir='14_06_18/', align = 'align/align_d_rms_1000_abs_t', 
              updateErr = True, poly='polyfit_nz/fit', points='points_nz/',
              chainsDir = 'efit/chains_S0-2_newRV2/', starlist='all', 
              plotSigAcc=False, magCut=16, nEpochs=30, pvalue=3.0,root_tmp='/g/ghez/align/'):
    """
    Just returns ar, are, at, and ate

    Inputs:
    alnDir   = The root directory of an astrometry analysis
                (e.g. '08_02_16/' or './' if you are in the directory).
    align     = The align root file name (including the directory relative
                to root). Make sure that polyfit was run on this align
		output.
    poly      = The polyfit root file name (including the directory relative
                to root). This should be run on the same align as above.
    points    = The points directory.
    starlist  = Only plot specific subset of stars. Must be 'oldstars', 
                'yngstars', 'all'.
    plotSigAcc = Set to True to plot the fits for stars with significant radial
                 accelerations.
    nEpochs    = Number of epochs a star must have been detected in for this analysis 

    """
#    outdir = root + alnDir + 'plots/'
    outdir = home + 'plots/'

    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)

    names = s.getArray('name')

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
        
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

    #Update errors in position
    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, alnDir=alnDir, chainsDir = chainsDir)

    # Make an epochs cut
    idx = np.where((cnt > nEpochs) & (mag < magCut))[0]
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
    names = [names[nn] for nn in idx]

    projV = np.sqrt(vx**2 + vy**2)

    # Get accelerations in the radial/tangential direction
    r = np.hypot(x0, y0)
    ar = ((ax*x0) + (ay*y0)) / r
    at = ((ax*y0) - (ay*x0)) / r
    are =  (axe*x0/r)**2 + (aye*y0/r)**2
    are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    are =  np.sqrt(are)
    ate =  (axe*y0/r)**2 + (aye*x0/r)**2
    ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
    ate =  np.sqrt(ate)

    pass_f = run_f_test(names,pvalue,alnDir,align,poly,points,verbose=True)
    pass_fj = run_fjerk_test(names,pvalue,alnDir,align,poly,points,verbose=True)

    return ar, are, at, ate, mag, r, pass_f, pass_fj, projV
#    return x0, x0e, y0, y0e, ax, axe, ay, aye, mag, r, pass_f, pass_fj

def quick_delta(accel, accel_err, guess):
    for ess in guess:
        tmp_sigma=accel*0.0
        for j in range(len(accel)):
            tmp_sigma[j]=accel[j]/np.sqrt(accel_err[j]**2+ess**2)
        print np.std(tmp_sigma)


def hazGauss(xdat, ydat, guess=[1,0,1]):
    #Because python doesn't have a standard Gauss fit function
    #pdb.set_trace()
    def gauss(params, x, y):
        return y-(params[0]*np.exp((x-params[1])**2/(-2*params[2]**2)))

    return opter.leastsq(gauss, guess, args=(xdat, ydat))

def histGauss(data, binsize=0.1, guess=[1,0,1]):
    #Will output best fit values for gauss of hist
    #binsize is size of the bins of histogram, can't make smaller than 1
    #guess is best guess for constant, x_0, and sigma in that order

    bottom = int(min(data))-binsize
    top = int(max(data))+binsize
    bins = np.arange(bottom, top, binsize)
    (n, b, p) = py.hist(data,bins)
    yout = n
    xout = n*0.0
    for i in range(len(xout)):
        xout[i] = (b[i]+b[i+1])/2.0
    
    output = hazGauss(xout, yout, guess)
    return output

def fitWeird(accel, error, mag, minD = 0.0, maxD = 5e-5, size = 5e-8, magRange=[9,17]):
#    num = int((maxD - minD)/size)
#    delta = np.array([minD+(size * i) for i in range(num)])
    delta = np.arange(minD,maxD,size)
    bf_out = delta*0.0
    sigmas = bf_out*0.0
    idx = np.where((mag > magRange[0]) & (mag <= magRange[1]) & (error < 6e-5))[0]
    accel = accel[idx]
    error = error[idx]
    for j in range(len(delta)):
        tmp_data = accel/np.sqrt(error**2+delta[j]**2)
        results = histGauss(tmp_data)
        sigmas[j]=results[0][2]
        bf_out[j]=abs(1-results[0][2])
    
    want = np.argmin(bf_out)
    print delta[want]
    print len(idx)
    return delta, sigmas

def fitDelta(x0,y0,x0e,y0e,ax,axe,ay,aye,mag,radius,minD=0.0,maxD=1e-4,size=1e-6,magRange=[9,17],radRange=[0.0,7.0]):
    
#    num = int((maxD - minD)/size)
#    delta = np.array([minD+(size * i) for i in range(num)])
    delta = np.arange(minD,maxD,size)
    bf_out = delta*0.0
    sigmas = bf_out*0.0
    idx = np.where((mag>magRange[0]) & (mag<=magRange[1]) & ((radius>radRange[0]) & (radius<=radRange[1])) & (x0e<8.3e-4) & (y0e<1.75e-3))[0]
    x0 = x0[idx]
    y0 = y0[idx]
    x0e = x0e[idx]
    y0e = y0e[idx]
    ax = ax[idx]
    ay = ay[idx]
    axe = axe[idx]
    aye = aye[idx]
    for j in range(len(delta)):
        x0e_t = np.sqrt(x0e**2 + delta[j]**2)
        y0e_t = np.sqrt(y0e**2 + (10.0*delta[j])**2)
        r = np.hypot(x0, y0)
        ar = ((ax*x0) + (ay*y0)) / r
        at = ((ax*y0) - (ay*x0)) / r
    #are =  (axe*x0/r)**2 + (aye*y0/r)**2
    #are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    #are =  np.sqrt(are)
        ate =  (axe*y0/r)**2 + (aye*x0/r)**2
        ate += (y0*x0e_t*ar/r**2)**2 + (x0*y0e_t*ar/r**2)**2
        ate =  np.sqrt(ate)
        tmp_data = at/ate
        results = histGauss(tmp_data)
        sigmas[j] = results[0][2]
        bf_out[j] = abs(1-results[0][2])

    want = np.argmin(bf_out)
    print delta[want]
    print len(idx)
    return delta, sigmas    


def fitWeirdXY(accel,error,mag,radius,minD = 0.0,maxD = 5e-5,size = 5e-8,magRange=[9,17],radRange=[0.0,100.]):
#    num = int((maxD - minD)/size)
#    delta = np.array([minD+(size * i) for i in range(num)])
    delta = np.arange(minD,maxD,size)
    bf_out = delta*0.0
    sigmas = bf_out*0.0
    idx = np.where((mag>magRange[0]) & (mag<=magRange[1]) & (error<6e-5) & ((radius>radRange[0]) & (radius<=radRange[1])))[0]
    accel = accel[idx]
    error = error[idx]
    for j in range(len(delta)):
        tmp_data = accel/np.sqrt(error**2+delta[j]**2)
        results = histGauss(tmp_data)
        sigmas[j]=results[0][2]
        bf_out[j]=abs(1-results[0][2])
    
    want = np.argmin(bf_out)
    print delta[want]
    print len(idx)
    return delta, sigmas



def scLeast2(x, y, error, slope=[-10.0,10.0,0.1], inter=[-10.0,10.0,0.1],guess=[1,0,1]):
    #takes the x and y you give it, takes the log of both, and fits a log/log y=mx+b
    #using least residuals

    error = np.log10((x+error)/x)
    x = np.log10(x)
    y = np.log10(y)

#    if (error == 0):
#        error = np.zeros(len(x))+1.0

    sRange=np.arange(slope[0],slope[1],slope[2])
    iRange=np.arange(inter[0],inter[1],inter[2])
    residuals=np.zeros((len(sRange),len(iRange)))
    sVal=np.zeros((len(sRange),len(iRange)))
    iVal=np.zeros((len(sRange),len(iRange)))
    for i in range(len(sRange)):
        for j in range(len(iRange)):
            sVal[i,j]=sRange[i]
            iVal[i,j]=iRange[j]
            for k in range(len(x)):
                residuals[i,j] += (y[k]-(x[k]*sRange[i]+iRange[j]))**2
    want = np.where(residuals == np.min(residuals))
    #pdb.set_trace()
    print "Best Fit Slope: %5.5f" %sVal[want]
    print "Best Fit Intercept: %5.5f" %iVal[want]

    probsi = np.exp(-residuals/2.0)
    Sprob = probsi.sum(axis=1)

#    py.clf()
#    py.plot(sRange,Sprob)
#    py.show()
#    pdb.set_trace()
    gaussInfo = hazGauss(sRange,Sprob,guess=guess)
    print 'Gauss Fit'
    print gaussInfo

    return sVal[want],iVal[want]

def hazGraphs():
    ar,are,at,ate,mag,r=getSigmas(alnDir='14_06_18/',poly='polyfit_nz/fit',points='points_nz/',magCut=22,nEpochs=22,updateErr=False)
    are_new = np.sqrt(are**2 + (3.1585e-5)**2)
    ate_new = np.sqrt(ate**2 + (3.6118e-5)**2)

    py.clf()
    bottom = int(min(at/ate))-1
    top = int(max(at/ate))+2
    bins = range(bottom, top, 1)
    xplot = np.array([bottom + 0.25*i for i in range(len(bins)*4)])
    bestp = histGauss(at/ate)
    yplot = np.array([bestp[0][0]*math.exp((xplot[i]-bestp[0][1])**2/(-2*bestp[0][2]**2)) for i in range(len(bins)*4)])
    py.hist(at/ate,bins)
    print bestp
    py.plot(xplot,yplot)
    py.xlabel('Sigma')
    py.title('Original Error')
    py.xlim([-30,30])
    py.savefig('/u/schappell/plots/sigma_tan_orig.png')

    py.clf()
    bottom_new = int(min(at/ate_new))-1
    top_new = int(max(at/ate_new))+2
    bins_new = range(bottom_new, top_new, 1)
    xplot_new = np.array([bottom_new + 0.25*i for i in range(len(bins_new)*4)])
    bestp_new = histGauss(at/ate_new)
    yplot_new = np.array([bestp_new[0][0]*math.exp((xplot_new[i]-bestp_new[0][1])**2/(-2*bestp_new[0][2]**2)) for i in range(len(bins_new)*4)])
    py.hist(at/ate_new,bins_new)
    print bestp_new
    py.plot(xplot_new,yplot_new)
    py.xlabel('Sigma')
    py.title('New Error')
    py.xlim([-30,30])
    py.savefig('/u/schappell/plots/sigma_tan_new.png')

    py.clf()
    py.plot(mag,ate,'.')
    py.plot([9,17],[3.6118e-5,3.6118e-5])
    py.xlabel('K mag')
    py.ylabel('Error in a_t')
    py.savefig('/u/schappell/plots/Tan_error_kmag.png')

#    py.clf()
#    delta_pl, sigma_pl = fitWeird(at,ate)
#    py.clf()
#    py.plot(delta_pl,sigma_pl)
#    py.xlabel('Delta Added')
#    py.ylabel('Sigma of Distribution')
#    py.savefig('/u/schappell/plots/sigma_v_delta.png')
    
def into8(mag):
    middle=np.median(mag)
    iL=np.where(mag < middle)[0]
    iR=np.where(mag >= middle)[0]

    left=mag[iL]
    right=mag[iR]
    Lmiddle=np.median(left)
    Rmiddle=np.median(right)
    iLL=np.where(left < Lmiddle)[0]
    iLR=np.where(left >= Lmiddle)[0]
    iRL=np.where(right < Rmiddle)[0]
    iRR=np.where(right >= Rmiddle)[0]

    leftLeft=left[iLL]
    leftRight=left[iLR]
    rightLeft=right[iRL]
    rightRight=right[iRR]
    LLmiddle=np.median(leftLeft)
    LRmiddle=np.median(leftRight)
    RLmiddle=np.median(rightLeft)
    RRmiddle=np.median(rightRight)

    print LLmiddle, Lmiddle, LRmiddle, middle, RLmiddle, Rmiddle, RRmiddle


def makeNSFtable():
    
    names,x0,y0,mag,r,ar,are,sigmaR,z0,z0e,az,aze,r3d,r3de,mext,cnt,vz_cms,vze_cms,pass_fj,jnames,jx,jy,jxe,jye,jxsig,jysig,xchi2r,ychi2r=histAccel(verbose=False,outTolatex=True)
#    inG = ['S0-1','S0-3','S0-19','S0-20','S0-16','S0-8','S0-26','S0-7','S0-4','S0-5','S0-40',
#           'S0-103','S0-49','S0-28','S0-27','S0-52','S0-54','S0-61','S0-30','S0-38','S0-17',
#           'S0-68','S0-70','S0-45','S1-2','S1-3','S0-36','S0-15','S1-12','S1-13']
    inG = ['S0-26','S0-7','S0-4','S0-40','S0-1','S0-3','S0-19','S0-2','S0-16','S0-8','S0-20',
           'S0-103','S0-28','S0-27','S0-52','S0-61','S0-30','S0-38','S0-17',
           'S0-68','S0-70','S1-2','S1-3','S0-36','S1-12','S1-13']
    inY = ['S1-3','S0-15','irs16C','irs16SW','S1-12','S1-14']
#    inB = ['S1-3','S0-15','S1-12']
    inB = ['S1-3','S1-12']
#    py.clf()
#    py.plot(cnt,ar/are,'.')
#    py.show()

#    new = np.empty(len(names),dtype='string')
    new = np.array(['New' for i in range(len(names))])
    for i in range(len(names)):
#        new[i] = 'New'
        for j in range(len(inG)):
            if names[i] == inG[j]:
                new[i] = 'G'
        for j in range(len(inY)):
            if names[i] == inY[j]:
                new[i] = 'Y'
        for j in range(len(inB)):
            if names[i] == inB[j]:
                new[i] = 'G+Y'
#        if (names[i] == ''):
#            names[i] = 'New'

    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
    cur.execute('SELECT * FROM spectra')
    vz = np.zeros(len(names))
    vze = np.zeros(len(names))
    vzt = np.zeros(len(names))
#    field = np.empty(len(names),dtype='string')
#    type = np.empty(len(names),dtype='string')
    field = np.array(['NA' for i in range(len(names))])
    type = np.array(['U' for i in range(len(names))])
    starsInC = ['S0-40','S1-97','S0-49','S1-92','S0-68','S0-61','S0-33','S1-26','S1-32','S0-35','S1-51','S1-4','S1-40','S0-54','S0-38',
                'star_770']
    TDold = ['S1-32','S0-35','S0-38']
    knownO = ['S0-49']
    knownY = ['S0-61']

#    pdb.set_trace()
#    for i in range(len(names)):
#        accNames = str(names[i])
#        try:
#            cur.execute('SELECT name,ddate,vz,vz_err FROM spectra WHERE name=?', [accNames])
#            for row in cur:
#                vz[i] = row[2]
#                vze[i] = row[3]
#                vzt[i] = row[1]
#                if (vze[i] != vze[i]):
#                    vze[i] = vz[i]
#        except:
#            print "Not in SQLite:"
#            print accNames
#        if (vzt[i] == 0.0):
#            vze[i] = 0.0
#        field[i] = 'NA'
    for i in range(len(names)):
        accNames = str(names[i])
#        try:
#            cur.execute('SELECT name,field FROM spectra WHERE name=?', [accNames])
#            for row in cur:
#                field[i] = row[1]
#        except:
#            continue
        try:
            cur.execute('SELECT name,field FROM spectra WHERE name=?', [accNames])
            for row in cur:
                field[i] = row[1]
        except:
            print "Not in SQLite spectr:"
            print accNames
#        for j in range(len(starsInNE)):
#            if (names[i] == starsInNE[j]):
#                field[i] = 'NE'
#        for j in range(len(starsInE)):
#            if (names[i] == starsInE[j]):
#                field[i] = 'E'
        for j in range(len(starsInC)):
            if (names[i] == starsInC[j]):
                field[i] = 'C'      
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
    cur.execute('SELECT * FROM stars')
    for i in range(len(names)):
        accNames = str(names[i])
        try:
            cur.execute('SELECT name,young,old,vz,vz_err,vz_ddate FROM stars WHERE name=?', [accNames])
            for row in cur:
#                print row[0],row[1],row[2]
                if (row[1] == 'T') | (row[2] == 'F'):
                    type[i] = 'Y'
                if (row[1] == 'F') | (row[2] == 'T'):
                    type[i] = 'O'
#                else:
#                    type[i] = 'U'
                vz[i] = row[3]
                vze[i] = row[4]
                vzt[i] = row[5]
                if (vze[i] != vze[i]):
                    vze[i] = vz[i]
        except:
#            type[i] = 'U'
            print "Not in SQLite stars:"
            print accNames
#        pdb.set_trace()
        for j in range(len(TDold)):
            if (accNames == TDold[j]):
                type[i] = 'O'
        if (accNames == 'S0-38'):
            vzt[i] = 2013.36
            vz[i] = -111.00
            vze[i] = 25.00
        if (accNames == 'S0-61'):
            type[i] ='Y'
        if (accNames == 'S0-49'):
            type[i] = 'O'
#        for j in range(len(knownO)):
#            if (accNames == knownO[j]):
#                type[i] = 'O'
#        for j in range(len(knownY)):
#            if (accNames == knownY[j]):
#                type[i] ='Y'

#    for i in range(len(names)):
#        print 'Name: %15s' % (names[i])
#        print 'Mag: %5.2f' % (mag[i])
#        print 'Radius: %5.3f' % (r[i])
#        print 'Accel rad: %8.3f and error %8.3f and sigma %8.2f' % (ar[i], are[i],sigmaR[i])
#        print 'Z velocity: %8.3f and error %8.3f and time %8.2f' % (vz[i], vze[i], vzt[i])
#        print 'Type: %3s' % (type[i])
#        print 'Field: %5s' % (field[i])
#        print 'New?: %5s' % (new[i])

    # Create the latex file
    out = open(home + 'tables/ksm14.tex','w')
    out.write('\\documentclass{aastex} \n')
    out.write('\\begin{singlespace} \n')
    out.write('\\begin{deluxetable}{lcccccccccccccc} \n')
#    out.write('\\rotate')
    out.write('\\rotate \n')
#    out.write('\\tabletypesize{\\scriptsize} \n')
    out.write('\\tabletypesize{\\tiny} \n')
    out.write('\\setlength{\\tabcolsep}{1.0mm} \n')
    out.write('\\tablewidth{0pt} \n')
    out.write('\\begin{document} \n')
    out.write('\\tablecaption{Significant Accelerating Sources}\n')
    out.write('\\tablehead{ \n')
    out.write('  \\colhead{Star} & \n')
    out.write('  \\colhead{$Kp$} & \n')
    out.write('  \\colhead{Epochs} & \n')
    out.write('  \\colhead{2D r} & \n')
#    out.write('  \\colhead{$T_o$} & \n')
#    out.write('  \\colhead{$X_o$} & \n')
#    out.write('  \\colhead{$Y_o$} & \n')
    out.write('  \\colhead{$a_r$} & \n')
#    out.write('  \\colhead{Error} & \n')
    out.write('  \\colhead{$a_r$} & \n')
#    out.write('  \\colhead{$v_z$} & \n')
#    out.write('  \\colhead{Error} & \n')
#    out.write('  \\colhead{Time} & \n')
    out.write('  \\colhead{z} & \n')
    out.write('  \\colhead{3D r} & \n')
    out.write('  \\colhead{$a_z$} & \n')
    out.write('  \\colhead{$v_z$} & \n')
#    out.write('  \\colhead{$M_{ext}$} & \n')
    out.write('  \\colhead{Type} & \n')
#    out.write('  \\colhead{Field} & \n')
    out.write('  \\colhead{Notes} & \n')
    out.write('  \\colhead{Spec} & \n')
    out.write('  \\colhead{Jerk} & \n')
#    out.write('  \\colhead{$a_t$} & \n')
#    out.write('  \\colhead{$a_t$} & \n')
#    out.write('  \\colhead{Signif-} \\\\ \n')
    out.write('%\n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(mag)} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
#    out.write('  \\colhead{(yr)} & \n')
#    out.write('  \\colhead{(arcsec)} & \n')
#    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{(micro-as/yr$^2$)} & \n')
#    out.write('  \\colhead{(mas/yr$^2$)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
#    out.write('  \\colhead{(km/s)} & \n')
#    out.write('  \\colhead{(km/s)} & \n')
#    out.write('  \\colhead{(yr)} & \n')
#    out.write('  \\colhead{(mas/yr$^2$)} & \n')
#    out.write('  \\colhead{(sigma)} & \n')
#    out.write('  \\colhead{icant?} \n')
    out.write('  \\colhead{(mpc)} & \n')
    out.write('  \\colhead{(mpc)} & \n')
    out.write('  \\colhead{(km/s/yr)} & \n')
    out.write('  \\colhead{(km/s)} & \n')
#    out.write('  \\colhead{($mM_{SMBH}$)} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{} & \n')
    out.write('} \n')
    out.write('\\startdata \n')

#    fmt = '%15s  %1s  %5.2f  %1s  %5.3f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %3s  %1s  %4s  %1s  %5s  %4s\n'


#    fmt = '%15s  %1s  %5.2f  %1s  %7.2f  %1s  %6.2f  %1s  %6.2f  %1s  %7.3f  %5s  '
#    fmt += '%6.3f  %1s  %5.2f  %1s  %7.3f  %5s  %6.3f  %1s  %5.2f  %1s  %1s  %4s\n'

#    pdb.set_trace()

    for ii in range(len(names)):
#        if (vzt[ii] == 0.0) | (vzt[ii] != vzt[ii]):
#             fmt = '%15s  %1s  %5.2f  %1s  %5.3f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %3s  %1s  %3s  %1s  %6.2f  %5s  %6.2f  '
#             fmt += '%1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %3s  %1s  %4s  %1s  %5s  %4s\n'
#             out.write(fmt % (names[ii], '&', mag[ii], '&', r[ii], '&', ar[ii]*1e3, '$\pm$', are[ii]*1e3, '&', sigmaR[ii], '&', '-','&', '-', '&', z0[ii]*1e3, '$\pm$', 
#                              z0e[ii]*1e3, '&', r3d[ii]*1e3, '$\pm$', r3de[ii]*1e3, '&', az[ii], '$\pm$', aze[ii], '&', mext[ii]*1e3, '&', type[ii], '&', field[ii], '&', new[ii], '\\\\'))
#        else:
        fmt = '%15s  %1s  %5.2f  %1s  %2d  %1s  %5.3f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  '
        fmt += '%5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %3s  %1s  %5s  %1s  %3s  %1s  %3s  %4s\n'
        if (vzt[ii] != 0.0) & (vzt[ii] == vzt[ii]):
            haveSpec = 'X'
        else:
            haveSpec = ' '
        if (pass_fj[ii] == 1):
            haveJerk = 'X'
        else:
            haveJerk = ' '
        out.write(fmt % (names[ii], '&', mag[ii], '&', cnt[ii], '&', r[ii], '&', ar[ii]*1e3, '$\pm$', are[ii]*1e3, '&', sigmaR[ii], '&', z0[ii]*1e3, '$\pm$', z0e[ii]*1e3, '&', r3d[ii]*1e3, '$\pm$', r3de[ii]*1e3, '&', az[ii], '$\pm$', aze[ii], '&', vz_cms[ii]*1e-5, '$\pm$', vze_cms[ii]*1e-5, '&', type[ii], '&', new[ii], '&', haveSpec, '&', haveJerk, '\\\\'))

    out.write('\\\\\n')
    out.write('\\enddata \n')
    out.write('\\end{deluxetable} \n')
    out.write('\\end{singlespace} \n')
    out.write('\\end{document} \n')
    out.close()


    out = open(home + 'tables/jerkTab.tex','w')
    out.write('\\documentclass{aastex} \n')
    out.write('\\begin{singlespace} \n')
    out.write('\\begin{deluxetable}{lccccccccccc} \n')
    out.write('\\rotate \n')
    out.write('\\tabletypesize{\\tiny} \n')
    out.write('\\setlength{\\tabcolsep}{1.0mm} \n')
    out.write('\\tablewidth{0pt} \n')
    out.write('\\begin{document} \n')
    out.write('\\tablecaption{Sources With Significant Jerk Terms}\n')
    out.write('\\tablehead{ \n')
    out.write('  \\colhead{Star} & \n')
    out.write('  \\colhead{$Kp$} & \n')
    out.write('  \\colhead{Epochs} & \n')
    out.write('  \\colhead{2D r} & \n')
    out.write('  \\colhead{$a_r$} & \n')
    out.write('  \\colhead{$a_r$} & \n')
    out.write('  \\colhead{$J_x$} & \n')
    out.write('  \\colhead{$J_x$} & \n')
    out.write('  \\colhead{$J_y$} & \n')
    out.write('  \\colhead{$J_y$} & \n')
    out.write('  \\colhead{Red. $\chi_x^2$} & \n')
    out.write('  \\colhead{Red. $\chi_y^2$} & \n')
    out.write('%\n')
    out.write('  \\colhead{} & \n')
#    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(mag)} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{(micro-as/yr$^2$)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
    out.write('  \\colhead{(micro-as/yr$^3$)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
    out.write('  \\colhead{(micro-as/yr$^3$)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{} & \n')
    out.write('} \n')
    out.write('\\startdata \n')

    jdex = 0

    printJ = np.where(pass_fj == 1)[0]
    for ii in printJ:
        fmt = '%15s  %1s  %5.2f  %1s  %2d  %1s  %5.3f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %6.2f  %1s  %6.2f  %4s\n'

        out.write(fmt % (names[ii], '&', mag[ii], '&', cnt[ii], '&', r[ii], '&', ar[ii]*1e3, '$\pm$', are[ii]*1e3, '&', sigmaR[ii], '&', jx[jdex]*1e6, '$\pm$', jxe[jdex]*1e6, '&', jxsig[jdex], '&', jy[jdex]*1e6, '$\pm$', jye[jdex]*1e6, '&', jysig[jdex], '&', xchi2r[jdex], '&', ychi2r[jdex], '\\\\'))

        jdex += 1

    out.write('\\\\\n')
    out.write('\\enddata \n')
    out.write('\\end{deluxetable} \n')
    out.write('\\end{singlespace} \n')
    out.write('\\end{document} \n')
    out.close()

    py.clf()
    knownY = np.where((new!='New') & (type=='Y'))[0]
    knownO = np.where((new!='New') & (type=='O'))[0]
    knownU = np.where((new!='New') & (type=='U'))[0]
    newY = np.where((new=='New') & (type=='Y'))[0]
    newO = np.where((new=='New') & (type=='O'))[0]
    newU = np.where((new=='New') & (type=='U'))[0]
#    pdb.set_trace()
    py.plot(r[knownY],-1*ar[knownY],'co',label='Known Young')
    py.plot(r[knownO],-1*ar[knownO],'ro',label='Kn Old')
    py.plot(r[knownU],-1*ar[knownU],'go',label='Kn Unknown')
    py.plot(r[newY],-1*ar[newY],'cd',label='New Young')
    py.plot(r[newO],-1*ar[newO],'rd',label='New Old')
    py.plot(r[newU],-1*ar[newU],'gd',label='New Unknown')
    pl1=py.plot([],[],'ko')
    pl2=py.plot([],[],'kd')
    py.yscale('log')
#    pl3=py.plot([],[],'r*',ms=7)
#    pl4=py.plot([],[],'ro',ms=7)
    py.legend((pl2,pl1),['Known','New'],numpoints=1)
    py.text(1.5,3.0,'Young',color='c',fontsize=12)
    py.text(1.5,2.5,'Old',color='r',fontsize=12)
    py.text(1.5,2.1,'Unknown',color='g',fontsize=12)
#    py.legend(numpoints=1)
    py.axis([0.2,1.8,0.04,10])
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('Radial Acceleration (mas/yr$^2$)')
    py.savefig('/u/schappell/allAccel_raccel_r.png')

    py.clf()
    knownY = np.where((new!='New') & (type=='Y'))[0]
    knownO = np.where((new!='New') & (type=='O'))[0]
    knownU = np.where((new!='New') & (type=='U'))[0]
    newY = np.where((new=='New') & (type=='Y'))[0]
    newO = np.where((new=='New') & (type=='O'))[0]
    newU = np.where((new=='New') & (type=='U'))[0]
    py.plot(r[knownY],mag[knownY],'co',label='Known Young')
    py.plot(r[knownO],mag[knownO],'ro',label='Old')
    py.plot(r[knownU],mag[knownU],'go',label='Unknown')
    py.plot(r[newY],mag[newY],'cd',label='New Young')
    py.plot(r[newO],mag[newO],'rd',label='New Old')
    py.plot(r[newU],mag[newU],'gd',label='New Unknown')
    pl1=py.plot([],[],'ko')
    pl2=py.plot([],[],'kd')
#    pl3=py.plot([],[],'r*',ms=7)
#    pl4=py.plot([],[],'ro',ms=7)
    py.legend((pl2,pl1),['New','Known'],loc=3,numpoints=1)
    py.text(0.05,10.4,'Young',color='c',fontsize=12)
    py.text(0.05,10.8,'Old',color='r',fontsize=12)
    py.text(0.05,11.2,'Unknown',color='g',fontsize=12)
#    py.legend(numpoints=1)
#    py.axis([0,2,0,2500])
    py.xlabel('Projected Radius (arcsec)')
#    py.ylabel('Radial Acceleration (micro-as/yr^2)')
    py.ylabel('K (mag)')
    py.savefig('/u/schappell/allAccel_mag_r.png')

    imgFile='/u/ghezgroup/data/gc/09maylgs/combo/mag09maylgs_sgra_dim_kp.fits'
    sgra=[624.5,726.3]
    scale = 0.00995
    img = pyfits.getdata(imgFile)
    imgsize = (img.shape)[0]

    # Make axes for images in arcsec
    pixL = np.arange(0,imgsize)
    xL = [-1*(xpos - sgra[0])*scale for xpos in pixL]
    yL = [(ypos - sgra[1])*scale for ypos in pixL]

    # Plots
    fig = py.figure(figsize=(8,8))
    fig.subplots_adjust(left=0.1,right=0.95,top=0.95)
    fig.clf()
    ax = fig.add_subplot(111)

    ax.imshow(np.log10(img+1), aspect='equal', interpolation='bicubic',
              extent=[max(xL), min(xL), min(yL), max(yL)],vmin=2.2,vmax=5,
              origin='lowerleft', cmap=py.cm.gray_r)

    for ii in range(len(x0)):
        specT = str(type[ii])
#        p1 = ax.plot(x0[ii],y0[ii],color='g',marker='o',mfc='g',ms=5)
        if (specT == 'O'):
            p1 = ax.plot(x0[ii],y0[ii],color='r',marker='o',mfc='r',ms=7)
        if (specT == 'U'):
            p2 = ax.plot(x0[ii],y0[ii],color='g',marker='D',mfc='g',ms=7)
#        if ((vzt[ii] == 0.0) | (vzt[ii] != vzt[ii])):
#            p3 = ax.plot(x0[ii],y0[ii],color='r',marker='o',mfc='None',mec='r',mew=1,ms=11)
        if (specT == 'Y'):
            p3 = ax.plot(x0[ii],y0[ii],color='c',marker='o',mfc='c',ms=7)
    # Set up legend
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=10)

#    legend_items = ['O/Y','Unknown','No V_z']
    legend_items = ['Old','Unknown SpT','Young']
    ax.legend((p1,p2,p3),legend_items, numpoints=1, loc=3, prop=prop)
    py.plot([0],[0],'k+',ms=8,mew=2)
    lgdLines = ax.get_legend().get_lines()
    py.setp(lgdLines, visible=False)
    ax.set_xlabel('RA (arcsec)')
    ax.set_ylabel('Dec (arcsec)')
    ax.axis([2,-2,-2,2])
    fig.savefig('/u/schappell/unknown_sample.png')

#    outz = np.where(new == 'New')[0]
#    names = [names[nn] for nn in outz]

    bandaid = np.zeros(len(names))
    for i in range(len(names)):
        nstar=str(names[i])
        if ((nstar=='irs16C') | (nstar=='irs16SW') | (nstar=='S0-49') | (nstar=='S0-40') | (nstar=='S0-36') | (nstar=='S0-61')):
            bandaid[i] = 1
        if ((nstar=='S0-61')):
            type[i]='Y'
        if (nstar=='S0-49'):
            type[i]='O'

#    pdb.set_trace()
    py.clf()
    have_o = np.where((((vzt != 0.0) & (vzt == vzt)) | (bandaid == 1)) & (type=='O') )[0]
#    have_oBB = np.where((((vzt != 0.0) & (vzt == vzt)) | (bandaid == 1)) & (type=='O') & (mag < 15.7))[0]
    have_u = np.where((((vzt != 0.0) & (vzt == vzt)) | (bandaid == 1)) & (type=='U'))[0]
    have_y = np.where((((vzt != 0.0) & (vzt == vzt)) | (bandaid == 1)) & (type=='Y'))[0]
    no_o = np.where(((vzt == 0.0) | (vzt != vzt)) & (type=='O') & (bandaid != 1))[0]
    no_u = np.where(((vzt == 0.0) | (vzt != vzt)) & (type=='U') & (bandaid != 1))[0]
    no_y = np.where(((vzt == 0.0) | (vzt != vzt)) & (type=='Y') & (bandaid != 1))[0]
    gr=py.figure()
    ax=gr.add_subplot(1,1,1)
    py.plot(r[no_o],-1*ar[no_o],'go',mfc='None',mec='g',ms=7)
    py.plot(r[no_u],-1*ar[no_u],'go',mfc='None',mec='g',ms=7)
    py.plot(r[no_y],-1*ar[no_y],'go',mfc='None',mec='g',ms=7)
    py.plot(r[have_o],-1*ar[have_o],'ro',ms=7)
    py.plot(r[have_u],-1*ar[have_u],'go',ms=7)
#    py.plot(r[have_y],-1*ar[have_y],'bo',ms=5)
#    py.plot(r[have_oBB],-1*ar[have_oBB],'r*',ms=7)
    pl1=py.plot([],[],'ro',ms=7)
    pl2=py.plot([],[],'go',mfc='None',ms=7)
#    pl3=py.plot([],[],'r*',ms=7)
#    pl4=py.plot([],[],'ro',ms=7)
    py.legend((pl2,pl1),['Unknown SpT/$V_z$ (Proposed sample)','Late-type (old)'],numpoints=1)
#    py.text(1.4,2.3,'--Young',color='b',fontsize=14)
#    py.text(1.4,3.4,'--Old',color='r',fontsize=14)
#    py.text(1.4,2.8,'--Unknown',color='g',fontsize=14)
#    py.axis([0,1.8,0,10])
    py.yscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    py.xlabel('Projected Radius (")')
    py.ylabel('Acceleration in-the-plane-of-the-sky (mas/$yr^2$)')
    py.savefig('/u/schappell/vzAccel_raccel_r.ps')

    out = open(home + 'tables/rv_2015.5.tex','w')
    out.write('\\documentclass{aastex} \n')
    out.write('\\begin{singlespace} \n')
    out.write('\\begin{deluxetable}{lccccccr} \n')
#    out.write('\\rotate')
#    out.write('\\rotate \n')
    out.write('\\tabletypesize{\\scriptsize} \n')
#    out.write('\\tabletypesize{\\tiny} \n')
    out.write('\\setlength{\\tabcolsep}{1.0mm} \n')
    out.write('\\tablewidth{0pt} \n')
    out.write('\\begin{document} \n')
    out.write('\\tablecaption{Significant Accelerating Sources}\n')
    out.write('\\tablehead{ \n')
    out.write('  \\colhead{Star} & \n')
    out.write('  \\colhead{Astrometric} & \n')
    out.write('  \\colhead{Spectroscopy} & \n')
#    out.write('  \\colhead{$Kp$} & \n')
#    out.write('  \\colhead{2D r} & \n')
#    out.write('  \\colhead{$a_r$} & \n')
#    out.write('  \\colhead{$a_r$} & \n')
    out.write('  \\colhead{Time} & \n')
    out.write('  \\colhead{$v_z$} & \n')
#    out.write('  \\colhead{z} & \n')
#    out.write('  \\colhead{3D r} & \n')
    out.write('  \\colhead{$a_z$} & \n')
    out.write('  \\colhead{$\Delta$$v_z$ (2015.5)} & \n')
    out.write('  \\colhead{$\Delta$$v_z$} & \n')
#    out.write('  \\colhead{$M_{ext}$} & \n')
#    out.write('  \\colhead{Type} & \n')
#    out.write('  \\colhead{Field} & \n')
#    out.write('  \\colhead{Notes} & \n')
    out.write('%\n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{Epochs} & \n')
    out.write('  \\colhead{Epochs} & \n')
#    out.write('  \\colhead{(mag)} & \n')
#    out.write('  \\colhead{(arcsec)} & \n')
#    out.write('  \\colhead{(yr)} & \n')
#    out.write('  \\colhead{(arcsec)} & \n')
#    out.write('  \\colhead{(arcsec)} & \n')
#    out.write('  \\colhead{(micro-as/yr$^2$)} & \n')
#    out.write('  \\colhead{(mas/yr$^2$)} & \n')
#    out.write('  \\colhead{(sigma)} & \n')
    out.write('  \\colhead{(yr)} & \n')
    out.write('  \\colhead{(km/s)} & \n')
#    out.write('  \\colhead{(km/s)} & \n')
#    out.write('  \\colhead{(mas/yr$^2$)} & \n')
#    out.write('  \\colhead{(sigma)} & \n')
#    out.write('  \\colhead{icant?} \n')
#    out.write('  \\colhead{(mpc)} & \n')
#    out.write('  \\colhead{(mpc)} & \n')
    out.write('  \\colhead{(km/s/yr)} & \n')
    out.write('  \\colhead{(km/s)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
#    out.write('  \\colhead{($mM_{SMBH}$)} & \n')
#    out.write('  \\colhead{} & \n')
#    out.write('  \\colhead{} & \n')
#    out.write('  \\colhead{} & \n')
    out.write('} \n')
    out.write('\\startdata \n')

    cur = connection.cursor()
    cur.execute('SELECT * FROM spectra')
    for ii in range(len(names)):
        accNames = str(names[ii])
        if (vzt[ii] != 0.0) & (vzt[ii] == vzt[ii]):
            deltaV = az[ii]*(2015.5-vzt[ii])
            snr = deltaV / vze[ii]
            cur.execute('SELECT name,ddate,vz,vz_err FROM spectra WHERE name=?', [accNames])
            nRV = 0
            for row in cur:
                try:
                    row[3]/row[2]
                    nRV += 1
                except (ValueError,TypeError):
                    nRV += 0
            if (accNames == 'S1-13'):
                nRV = 1
            if (accNames == 'S0-38'):
                nRV = 2
            fmt = '%15s  %1s  %2d  %1s  %2d  %1s  %5.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %1s  %6.2f  %1s  %6.2f  %4s\n'
            out.write(fmt % (names[ii], '&', cnt[ii], '&', nRV, '&', vzt[ii], '&', vz[ii], '$\pm$', vze[ii], '&', az[ii], '&', deltaV, '&', snr, '\\\\'))

    out.write('\\\\\n')
    out.write('\\enddata \n')
    out.write('\\end{deluxetable} \n')
    out.write('\\end{singlespace} \n')
    out.write('\\end{document} \n')

    out.close()



    return names, az, aze



def zAccel():
    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
    cur.execute('SELECT * FROM spectra')
    namez,az,aze = makeNSFtable()
    for i in range(len(namez)):
        accNames = str(namez[i])
        time=[]
        vel=[]
        vele=[]
        cur.execute('SELECT name,ddate,vz,vz_err FROM spectra WHERE name=?', [accNames])
        for row in cur:
            print row[2], row[3]
#            if (row[3] != 'None') & (row[3] != '') & (row[2] != 'Nonw'):
            try:
                vele.append(float(row[3]))
                vel.append(float(row[2]))
                time.append(float(row[1]))
#            if (row[3] != row[3]):
#                vele.append(float(row[2]))
#            else:
#                vele.append(float(row[3]))
            except (ValueError,TypeError):
                print " "
        time = np.array(time)
        vel = np.array(vel)
        vele = np.array(vele)
#        print "Time"
#        print time
#        print "Velocity"
#        print vel
#        print "Error"
#        print vele
        deltaV = az[i] * (2015.5 - time)
#        finalV = deltaV * (2015.5 - time) + vel
#        snr = deltaV / vele
        print '%15s in 2015.5:' % (namez[i])
        print '%5s  %18s  %19s  %5s' % ('Year','Vel Error (km/s)', 'Vel 2015.5 (km/s)', 'SNR')
        for j in range(len(time)):
            finalvel = vel[j] + deltaV[j]
            sigma = deltaV[j] / vele[j]
            print '%6.2f  %8.2f  %20.2f  %10.2f' % (time[j], vele[j], finalvel, sigma)



def compGill():
    gAfile = '/u/schappell/gillessen.txt'
    aInfo = asciidata.open(gAfile)

    #pdb.set_trace()
    aname = aInfo[0].tonumpy()
    amag = aInfo[1].tonumpy()
    ax = aInfo[2].tonumpy()   #In milli-arcsec
    axe = aInfo[3].tonumpy()
    avx = aInfo[4].tonumpy()  #In milli-arcsec/yr
    avxe = aInfo[5].tonumpy()
    aax = aInfo[6].tonumpy()
    aax = aax * 2.0    #The 1/2 is in Gillessen reported values for accelerations
    aaxe = aInfo[7].tonumpy()
    aaxe = aaxe * 2.0
#    jx = aInfo[8].tonumpy()
#    jxe = aInfo[9].tonumpy()
    at = aInfo[10].tonumpy()   #In years
    ay = aInfo[11].tonumpy()
    aye = aInfo[12].tonumpy()
    avy = aInfo[13].tonumpy()
    avye = aInfo[14].tonumpy()
    aay = aInfo[15].tonumpy()
    aay = aay * 2.0
    aaye = aInfo[16].tonumpy()
    aaye = aaye * 2.0
#    jy = aInfo[17].tonumpy()
#    jye = aInfo[18].tonumpy()
    atv = aInfo[19].tonumpy()  #In years
    avz = aInfo[20].tonumpy()   #In km/s
    avze = aInfo[21].tonumpy()
    aaz = aInfo[22].tonumpy()
#    aaz = aaz * 2.0     #Not all stars have accel assosicated with their velocity in z, will be listed as 'n' as result, deal with this before use
#    aaze = aInfo[23].tonumpy()
#    aaze = aaze * 2.0
#    ajz = aInfo[24].tonumpy()  #STILL NEED to account for the 1/6 factor, do this before using
#    ajze = aInfo[25].tonumpy()
    ourAnames = aInfo[26].tonumpy()   #the names of these stars in our system

    r = np.hypot(ax,ay)
    aar = ((aax*ax) + (aay*ay)) / r
    aat = ((aax*ay) - (aay*ax)) / r
    aare =  (aaxe*ax/r)**2 + (aaye*ay/r)**2
    aare += (ay*axe*aat/r**2)**2 + (ax*aye*aat/r**2)**2
    aare =  np.sqrt(aare)
    aate =  (aaxe*ay/r)**2 + (aaye*ax/r)**2
    aate += (ay*axe*aar/r**2)**2 + (ax*aye*aar/r**2)**2
    aate =  np.sqrt(aate)

    hdr = '%15s  %7s  %5s  %6s  %17s  %14s  %17s  %14s'
    fmt = '%15s  %7s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %8.3f mas/yr^2  %8.2f sigma '
    print hdr % ('Name','G Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)','a_tan (mas/yr^2)','a_tan (sigma)')
    for i in range(len(aname)):
        print fmt % (ourAnames[i],aname[i],amag[i],r[i]*1e-3,aar[i],aar[i]/aare[i],aat[i],aat[i]/aate[i])

#    return aar,aare,aat,aate
    return ourAnames


def plotDelta(alnDir='14_06_18/',
              align = 'align/align_d_rms_1000_abs_t', updateErr = True,
              poly='polyfit_nz/fit', points='points_nz/', sigma=3.0,
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains_S0-2_newRV2/',
              starlist='all', plotSigAcc=False, magCut=22, nEpochs=24, verbose = True):
    outdir = home + 'plots/'

    s = loadPop(alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)

    names = s.getArray('name')
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
        
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

    #Update errors in position
    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, alnDir=alnDir, chainsDir = chainsDir)

    # Make an epochs cut
    idx = np.where((cnt > nEpochs) & (mag < magCut))[0]
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
    names = [names[nn] for nn in idx]

    # Get accelerations in the radial/tangential direction
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
        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17.1])
        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5
        delta = mag*0.0
        for i in range(len(mag)):
            for j in range(len(atBins)-1):
                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
                    delta[i] = deltaArr[j]

        ateN = np.sqrt(ate**2 + delta**2)
        areN = np.sqrt(are**2 + delta**2)
    magB = np.zeros(len(deltaArr))
    medB = np.zeros(len(deltaArr))
    medA = np.zeros(len(deltaArr))
    for i in range(len(deltaArr)):
        magB[i] = (atBins[i] + atBins[i+1])/2.0
        inhere = np.where((mag > atBins[i]) & (mag <= atBins[i+1]))[0]
        medB[i] = median(ate[inhere])
        medA[i] = median(ateN[inhere])

    py.clf()
#    py.plot(magB,medB*1e6,'.',label='Before')
#    py.plot(magB,deltaArr*1e6,'o',label='Delta')
    py.plot(magB,medA*1e6,'.',label='After')
#    py.legend(loc=2)
    py.axis([10,17,15,60])
    py.xlabel('K Mag')
    py.ylabel('Error (micro-as/yr^2)')
    py.savefig('/u/schappell/delta_before_after.png')


def calcInfo(starNames=['S0-1','S0-3','S0-19','S0-20','S0-16','S0-8','S0-26','S0-7','S0-4','S0-5','S0-40','S0-103','S0-49','S0-52','S0-28','S0-27','S0-54','S0-61','S0-30','S0-38','S0-45','S0-17','S0-68','S0-70','S1-2','S1-3','S0-36','S0-15','S1-12','S1-13'],
             alnDir='14_06_18/',root_tmp='/g/ghez/align/',
             align = 'align/align_d_rms_1000_abs_t', updateErr = True,
             poly='polyfit_nz/fit', points='points_nz/', sigma=3.0,
             f_test=True, pvalue=4.0, chainsDir = 'efit/chains_S0-2_newRV2/',
             starlist='all', plotSigAcc=False, verbose = True):
    """
    Gives info on the star(s) that input
    """
#    outdir = root + alnDir + 'plots/'
    outdir = home + 'plots/'

    s = loadPop(alnDir=alnDir,root_tmp=root_tmp,starlist=starlist,align=align,poly=poly,points=points)
    names = s.getArray('name')

#    pdb.set_trace()
    ii = np.zeros(len(starNames),dtype=int)
#    pdb.set_trace()
    for i in range(len(starNames)):
#    names = s.getArray('name'
        ii[i] = names.index(str(starNames[i]))
#    star = s.stars[ii]

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
        
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

    #Update errors in position
    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, alnDir=alnDir, chainsDir = chainsDir)

    # Get accelerations in the radial/tangential direction
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
        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17.1])
        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5
        delta = mag*0.0
        for i in range(len(mag)):
            for j in range(len(atBins)-1):
                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
                    delta[i] = deltaArr[j]

        ate = np.sqrt(ate**2 + delta**2)
        are = np.sqrt(are**2 + delta**2)


    # Lets also do parallel/perpendicular to velocity
#    v = np.sqrt(vx**2 + vy**2)
#    am = ((ax*vx) + (ay*vy)) / v
#    an = ((ax*vy) - (ay*vx)) / v
#    ame = np.sqrt((axe*vx)**2 + (aye*vy)**2) / v
#    ane = np.sqrt((axe*vy)**2 + (aye*vx)**2) / v

    # Total acceleration
#    atot = py.hypot(ax, ay)
#    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

    ##########
    #
    # KS Test... Are the distributions Normal?
    #
    ##########
    sigmaR = ar/are
    sigmaT = at/ate
#    sigmaAll = concatenate([sigmaR, sigmaT])
#    (ksdR, kspR) = stats.stats.kstest(sigmaR, 'norm', N=len(sigmaR))
#    (ksdT, kspT) = stats.stats.kstest(sigmaT, 'norm', N=len(sigmaT))
#    (ksdA, kspA) = stats.stats.kstest(sigmaAll, 'norm', N=len(sigmaAll))
    #print 'KS Test for normality (prob that observed is gaussian):'
    #print '\tRadial Acc. KS Prob. = %5.3f' % (kspR)
    #print '\tTangen Acc. KS Prob. = %5.3f' % (kspT)
    #print '\tCombo Acc. KS Prob. = %5.3f' % (kspA)
    print ''

    ###########
    #
    # F Test -- Accel or Velocity fit?
    #
    ###########
    if f_test == True:
        # Pass in the star names and run the F test
        pass_f,xprob,yprob = run_f_test(names,pvalue,root_tmp,alnDir,align,poly,points,verbose=verbose)
        pass_fj = run_fjerk_test(names,pvalue,root_tmp,alnDir,align,poly,points,verbose=verbose)
    else:
        pass_f = np.ones(len(names)) # check everything
        pass_fj = np.ones(len(names))


    #Calculating z0 aand a_z and their errors
    x0_pc = x0 * dist / au_in_pc
    y0_pc = y0 * dist / au_in_pc
    x0_cm = x0_pc * cm_in_pc
    y0_cm = y0_pc * cm_in_pc
    x0e_pc = x0e * dist / au_in_pc
    y0e_pc = y0e * dist / au_in_pc
    x0e_cm = x0e_pc * cm_in_pc
    y0e_cm = y0e_pc * cm_in_pc
#    r2d = np.sqrt(x0**2 + y0**2) # arcsec
#    r2d_pc = r2d * dist / au_in_pc
    r_cm = r * dist * cm_in_pc / au_in_pc
    re_cm = np.sqrt((x0_cm*x0e_cm/r_cm)**2 + (y0_cm*y0e_cm/r_cm)**2)
    ar_cmss = ar * asy_to_kms * 1e5 / sec_in_yr
    are_cmss = are * asy_to_kms * 1e5 / sec_in_yr
    vx_cms = vx * asy_to_kms * 1e5
    vy_cms = vy * asy_to_kms * 1e5
    vxe_cms = vxe * asy_to_kms * 1e5
    vye_cms = vye * asy_to_kms * 1e5
    vproj_cms = np.sqrt(vx_cms**2 + vy_cms**2)
    vproje_cms = np.sqrt((vx_cms*vxe_cms)**2 + (vy_cms*vye_cms)**2)/vproj_cms
    z0_cm = np.sqrt(abs(((GM * r_cm / ar_cmss)**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)

    az_cmss = GM * z0_cm / ((x0_cm**2 + y0_cm**2 + z0_cm**2)**(1.5))

    z0e_cm = (abs(G*r_cm/ar_cmss))**(4.0/3.0)*(mass_g)**(-2.0/3.0)*masse_g**2 + ((abs(GM/ar_cmss))**(2.0/3.0)*(abs(r_cm))**(-1.0/3.0)-r_cm)**2*re_cm**2
    z0e_cm += (abs(GM*r_cm))**(4.0/3.0)*(abs(ar_cmss))**(-10.0/3.0)*are_cmss**2
    z0e_cm = np.sqrt(z0e_cm)/(3.0 * z0_cm)

    r3d_cm = np.sqrt(x0_cm**2 + y0_cm**2 + z0_cm**2)
    r3de_cm = np.sqrt(((x0_cm*x0e_cm)**2 + (y0_cm*y0e_cm)**2 + (z0_cm*z0e_cm)**2)/r3d_cm**2)
    r3d_pc = r3d_cm / cm_in_pc
    r3de_pc = r3de_cm / cm_in_pc
    mext = 10.0*pi*density0*(r3d_pc**2)
    mext_smbh = mext / mass

    aze_cmss = (z0_cm*masse_g)**2 + ((3.0/2.0)*mass_g*z0_cm*r_cm*re_cm/(r3d_cm**2))**2 + (mass_g*z0e_cm*(1.0-(3.0*z0_cm**2/r3d_cm**2)))**2
    aze_cmss = np.sqrt(aze_cmss) * G / r3d_cm**3

    az_kmsyr = az_cmss * sec_in_yr / 1e5
    aze_kmsyr = aze_cmss * sec_in_yr / 1e5
    ax_kmsyr = ax * asy_to_kms
    ay_kmsyr = ay * asy_to_kms
    axe_kmsyr = axe * asy_to_kms
    aye_kmsyr = aye * asy_to_kms
    a_tot = np.sqrt(az_kmsyr**2 + ax_kmsyr**2 + ay_kmsyr**2)
    ae_tot = np.sqrt((aze_kmsyr*az_kmsyr)**2 + (axe_kmsyr*ax_kmsyr)**2 + (aye_kmsyr*ay_kmsyr)**2)/a_tot

    z0_pc = z0_cm / cm_in_pc
    z0e_pc = z0e_cm / cm_in_pc

    vz_cms = (GM/r3d_cm) - vproj_cms**2
    print "Stars that are unbound:"
    for jj in ii:
        if (vz_cms[jj] < 0.0):
            print names[jj]
    vz_cms = np.sqrt(abs(vz_cms))
    vze_cms = np.sqrt((G*masse_g/(2*r_cm))**2 + (GM*re_cm/(2*r_cm**2))**2 + (vproj_cms*vproje_cms)**2)/vz_cms

    hdr = '%15s  %5s  %6s  %17s  %14s  %17s  %14s  %8s  %8s %8s  %7s  %7s'
    fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %8.3f mas/yr^2  %8.2f sigma  %2d epochs  %5s  %5s  %8.3f  %8.3f'
#    idsig = (np.where((sigmaR < -sigma) & (np.abs(sigmaT) < sigma)))[0] #sig accel
#    idnp = (np.where(((sigmaR > sigma) | (np.abs(sigmaT) > sigma)) & (cnt > 23)))[0] #non physical accel
#    writeLatex = (np.where((sigmaR < -sigma) & (np.abs(sigmaT) < sigma) & (pass_f == 1)))[0]
#    idsig = (np.where(((sigmaR < -22.0) | ((cnt > 24) & (sigmaR < -sigma))) & (np.abs(sigmaT) < sigma)))[0] #sig accel
#    writeLatex = (np.where(((sigmaR < -22.0) | ((cnt > 24) & (sigmaR < -sigma))) & (np.abs(sigmaT) < sigma) & (pass_f == 1)))[0]

    print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)','a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs', 'A/V Fit', 'J/A Fit', 'X PVal', 'Y PVal')
    if len(ii) > 0:
        for i in ii:
            if pass_f[i] == 1:
                fit = 'Acc'
            elif pass_f[i] == 0:
                fit = 'Vel'
            if pass_fj[i] == 1:
                jfit = 'Jrk'
            elif pass_fj[i] == 0:
                jfit = 'Acc'
            print fmt % (names[i], mag[i], r[i], ar[i]*1e3, sigmaR[i], at[i]*1e3, sigmaT[i], cnt[i], fit, jfit, xprob[i], yprob[i])
    hdr = '%15s  %6s  %7s  %10s  %7s  %12s  %7s  %15s  %7s  %17s  %7s '
    fmt = '%15s  %6.3f  %4.3f pc  %8.3f  %6.3f pc  %8.3f  %6.3f km/s  %8.3f  %6.3f km/s/yr  %8.3f  %6.3f km/s/yr'
    print hdr % ('Name','z (mpc)','error','3D r (mpc)','error','v_z (km/s)','error','a_z (km/s/yr)','error','a_tot (km/s/yr)','error')
    for i in ii:
        print fmt % (names[i],z0_pc[i]*1e3,z0e_pc[i]*1e3,r3d_pc[i]*1e3,r3de_pc[i]*1e3,vz_cms[i]*1e-5,vze_cms[i]*1e-5,az_kmsyr[i],aze_kmsyr[i],a_tot[i],ae_tot[i])

    gar = np.array([-1.22383614,-0.28086215,-0.83028028,-1.35719794,-1.79287078,-13.06707682,-2.29810387,-2.12296267,-0.43965733,-0.23296113,-1.30797324,-1.66108454,-0.60035385,-10.55786155,-5.44479274,-8.53332086,-1.10424009,-1.19682068,-0.09531442,-0.13620346,-0.701943,-0.16460509,-0.0824143,-0.13303168])

    gare = np.array([0.1644098,0.05243817,0.01859352,0.03893692,0.23743405,0.36248114,0.08768124,0.38602278,0.038,0.03614995,0.07198101,0.12188026,0.03723603,1.03278454,0.44558917,1.10195083,0.19600005,0.16222211,0.01601008,0.01758228,0.12285136,0.01019788,0.01246004,0.01482086])

    gat = np.array([8.80289986e-02,9.56684468e-02,-1.43758109e-02,-1.54601910e-01,-6.04594385e-01,8.71827579e-01,3.96407105e-01,-4.77876050e-01,1.19751744e-03,-1.05422936e-03,2.30047847e-01,-1.09536121e-01,-1.05883212e-01,-8.37185477e-01,8.93932891e-01,3.40190169e+00,-1.80847498e-01,2.73364718e-01,-2.55178605e-02,-9.41371920e-03,1.39985800e-01,4.99315900e-02,2.06853272e-02,-7.11149823e-03])

    gate = np.array([0.18069513,0.06162996,0.01945393,0.04114693,0.28603533,0.44012273,0.12999281,0.31356865,0.03800065,0.0378576,0.06204998,0.08193667,0.03700958,1.47020996,0.29753433,2.01897231,0.19600187,0.14706671,0.01999194,0.01645795,0.11063268,0.01183234,0.0118637,0.01524279])

    names = [names[ij] for ij in ii]
    ar = ar[ii]
    are = are[ii]
    at = at[ii]
    ate = ate[ii]

    inours = ['S0-54','S0-38','S0-40','S0-49','S0-68','S0-61','S0-30','S0-36','S0-8','S0-28','S0-5','S0-3','S0-4','S1-13','S0-15']

    py.clf()
#    py.errorbar(gar,ar*1e3,xerr=gare,yerr=are*1e3,fmt='.',color='b',ecolor='b')
    gsig = gar/gare
    osig = ar/are
    py.plot(gsig,osig,'bo')
    getline = [-100,100]
    py.plot(getline,getline,color='g')
    for i in range(len(names)):
        for j in range(len(inours)):
            if (names[i] == inours[j]):
#                py.errorbar(gar[i],ar[i]*1e3,xerr=gare[i],yerr=are[i]*1e3,fmt='.',color='r',ecolor='r')
                py.plot(gsig[i],osig[i],'ro')
#    py.axis([-])
    pl1=py.plot([],[],'ro')
    py.axis([-60,0,-100,40])
    legend(pl1,['In Sample'],loc=3,numpoints=1)
    py.xlabel('Gillessen a_r (sigma)')
    py.ylabel('a_rad (sigma)')
    savefig('/u/schappell/our_gillessen_arad.png')

    py.clf()
#    py.errorbar(gat,at*1e3,xerr=gate,yerr=ate*1e3,fmt='.')
    gsig = gat/gate
    osig = at/ate
    py.plot(gsig,osig,'bo')
    getline = [-100,100]
    py.plot(getline,getline,color='g')
    for i in range(len(names)):
        for j in range(len(inours)):
            if (names[i] == inours[j]):
#                py.errorbar(gat[i],at[i]*1e3,xerr=gate[i],yerr=ate[i]*1e3,fmt='.',color='r',ecolor='r')
                py.plot(gsig[i],osig[i],'ro')
    py.axis([-15,15,-15,15])
#    pl2=py.plot([],[],'r.')
    legend(pl1,['In Sample'],loc=3,numpoints=1)
    py.xlabel('Gillessen a_t (sigma)')
    py.ylabel('a_tan (sigma)')
    savefig('/u/schappell/our_gillessen_atan.png')
#    print ' %7s  %7s  %8s  %8s  %8s  %8s  ' %('G Name','O Name','G Sig T','G Sig R','O Sig T','O Sig R')
#    for i in range(len(names)):
    print gar/gare,gat/gate,ar/are,at/ate



def twoProp15():

    imgFile='/u/ghezgroup/data/gc/09maylgs/combo/mag09maylgs_sgra_dim_kp.fits'
    sgra=[624.5,726.3]
    scale = 0.00995
    img = pyfits.getdata(imgFile)
    imgsize = (img.shape)[0]

    # Make axes for images in arcsec
    pixL = np.arange(0,imgsize)
    xL = [-1*(xpos - sgra[0])*scale for xpos in pixL]
    yL = [(ypos - sgra[1])*scale for ypos in pixL]

    # Plots
    fig = py.figure(figsize=(8,8))
    fig.subplots_adjust(left=0.1,right=0.95,top=0.95)
    fig.clf()
    ax = fig.add_subplot(111)

    ax.imshow(np.log10(img+1), aspect='equal', interpolation='bicubic',
              extent=[max(xL), min(xL), min(yL), max(yL)],vmin=2.2,vmax=5,
              origin='lowerleft', cmap=py.cm.gray_r)

#    for ii in range(len(x0)):
#        specT = str(type[ii])
#        p1 = ax.plot(x0[ii],y0[ii],color='g',marker='o',mfc='g',ms=5)
#        if (specT == 'U'):
#            p2 = ax.plot(x0[ii],y0[ii],color='m',marker='D',mfc='m',ms=5)
#        if ((vzt[ii] == 0.0) | (vzt[ii] != vzt[ii])):
#            p3 = ax.plot(x0[ii],y0[ii],color='r',marker='o',mfc='None',mec='r',mew=1,ms=11)
        
    # Set up legend
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=10)

    py.plot([1.5,-.7,-1.15,1.05,1.5],[.6,1.6,.65,-.35,.6],'r')
    py.plot([0.95,-1.25,-1.7,.5,.95],[-.64,.36,-.66,-1.66,-.64],'r')
#    py.plot([-.6,-1.56,-2.52,-1.56,-.6],[1.4,1.82,-.4,-.82,1.4],'r')
    py.plot([1.9,-1.3,-2.1,1.1,1.9],[.8,1.8,-.8,-1.8,.8],'g--')
    py.text(0.0,.7,'North',color='m')
    py.text(-0.2,-.7,'South',color='m')
#    py.text(-1.6,.35,'West',color='m')
    p1=py.plot([],[],'r',label='Kn5')
    p2=py.plot([],[],'g--',label='Kn3')
    ax.legend((p1,p2),['Kn5','Kn3'],loc=3)
    py.plot([0],[0],'k+',ms=7,mew=2)
#    lgdLines = ax.get_legend().get_lines()
#    py.setp(lgdLines, visible=False)

    ax.set_xlabel('RA (arcsec)')
    ax.set_ylabel('Dec (arcsec)')
    ax.axis([2.2,-2.2,-2.,2.])
    fig.savefig('/u/schappell/pointings.png')

    gamma=np.array([-1.0,-0.5,0.0,0.5,1.0,1.5])#,2.0])
    naccelo=np.array([0.0,0.1,0.2,0.4,0.9,2.2])#,5.7])
    uppero=np.array([0.3,0.4,0.8,1.1,2.0,2.9])#,3.8])
    lowero=np.array([0.0,0.1,0.2,0.4,0.9,2.2])#,3.8])

    naccel = naccelo*2.2*8.5
    upper=uppero*sqrt(2.2*8.5)
    lower=lowero*sqrt(2.2*8.5)

    py.clf()
#    new=py.errorbar(gamma,naccel,yerr=[lower,upper],fmt='ro',ecolor='r',capsize=1.5,elinewidth=1.5,label='$N_a$ (mag < 17)')
    new=py.errorbar(gamma,naccel,yerr=[lower,upper],fmt='ro',ecolor='r',label='$N_a$ (mag < 17)')
    old=py.errorbar(gamma,naccelo,yerr=[lowero,uppero],fmt='ko',ecolor='k',label='$N_a$ (mag < 15.5)')
    ours=py.plot([-3,3],[2,2],color='m',label='Observed $N_a$ (<15.5 mag)')
    tdo=fill([0.05-.6,0.05-.6,.05+.29,.05+.29],[-15,145,145,-15],'b',alpha=0.2,label='Velocity Fit')
    bw=fill([7.0/4.0,7.0/4.0,3.0/2.0,3.0/2.0],[-15,145,145,-15],'g',alpha=0.2,label='Bahcall-Wolf')
#    tdo=py.axvspan(0.05-.6,0.06+.29,facecolor='b',alpha=0.5)
#    bw=py.axvspan(1.5,1.75,facecolor='g',alpha=0.5)
    legend(loc=2,numpoints=1)
    py.axis([-1.1,1.8,-5,60])
    py.xlabel('Gamma')
    py.ylabel('Number of Accel. Old Stars')
    py.savefig('/u/schappell/numaccel.png')



def whatAjerk(alnDir='14_06_18/',poly='polyfit_nzj/fit',points='points_nz/',printThese=['S0-2','S0-3','S0-1','S0-5','S0-20','S0-19','S0-8','S0-49','S0-38']):

#    pdb.set_trace()
    fitRoot = root + alnDir + poly
    fitFile = fitRoot + '.accelFormal'
    t0File = fitRoot + '.t0'

    _fit = asciidata.open(fitFile)
    _t0 = asciidata.open(t0File)

    t0x = _t0[1].tonumpy()
    t0y = _t0[2].tonumpy()
    names = _fit[0].tonumpy()
    x0 = _fit[1].tonumpy()
    vx = _fit[2].tonumpy()
    ax = _fit[3].tonumpy()
    jx = _fit[4].tonumpy()
    x0e = _fit[5].tonumpy()
    vxe = _fit[6].tonumpy()
    axe = _fit[7].tonumpy()
    jxe = _fit[8].tonumpy()
    xchi2 = _fit[9].tonumpy()
    xq = _fit[10].tonumpy()

    y0 = _fit[11].tonumpy()
    vy = _fit[12].tonumpy()
    ay = _fit[13].tonumpy()
    jy = _fit[14].tonumpy()
    y0e = _fit[15].tonumpy()
    vye = _fit[16].tonumpy()
    aye = _fit[17].tonumpy()
    jye = _fit[18].tonumpy()
    ychi2 = _fit[19].tonumpy()
    yq = _fit[20].tonumpy()

    #Transfer jerk from x and y to radial and tangential
    r = np.hypot(x0,y0)
    jr = ((jx*x0) + (jy*y0)) / r
    jt = ((jx*y0) - (jy*x0)) / r
    jre = (jxe*x0/r)**2 + (jye*y0/r)**2 + (y0*x0e*jt/r**2)**2 + (x0*y0e*jt/r**2)**2
    jre = np.sqrt(jre)
    jte = (jxe*y0/r)**2 + (jye*x0/r)**2 + (y0*x0e*jr/r**2)**2 + (x0*y0e*jr/r**2)**2
    jte = np.sqrt(jte)

    sigr = jr/jre
    sigt = jt/jte
    sigx = jx/jxe
    sigy = jy/jye

    fmt = '%15s  %15.3f  %15.2f  %15.3f  %15.2f  %10.2f  %7.2f  %7.2f'
    hdr = '%15s  %18s  %14s  %18s  %14s  %14s  %6s'

    def jerkAccel_equ(origin):
        xbh,ybh,zbh,z0,vz,az = origin
        radius = math.sqrt((x_prime-xbh)**2+(y_prime-ybh)**2+(z0-zbh)**2)
        jx_equ = GM_mas_yr*vx_prime/radius**3-3*ax_prime*((-x_prime+xbh)*vx_prime+(-y_prime+ybh)*vy_prime+(-z0+zbh)*vz)/radius**2-jx_prime
        jy_equ = GM_mas_yr*vy_prime/radius**3-3*ay_prime*((-x_prime+xbh)*vx_prime+(-y_prime+ybh)*vy_prime+(-z0+zbh)*vz)/radius**2-jy_prime
#        jz_equ = GM_mas_yr*vz/radius**3 - 3*az*((x_prime - xbh)*vx_prime+(y_prime - ybh)*vy_prime + (z0 - zbh)*vz)/radius**2 -jz
        ax_equ = -GM_mas_yr*(x_prime-xbh)/radius**3-ax_prime
        ay_equ = -GM_mas_yr*(y_prime-xbh)/radius**3-ay_prime
        az_equ = -GM_mas_yr*(z0-zbh)/radius**3-az
        atan_equ=-1.0*(ay_prime*(z0-zbh)-az*(y_prime-ybh)+az*(x_prime-xbh)-ax_prime*(z0-zbh)+ax_prime*(y_prime-ybh)-ay_prime*(x_prime-xbh))
        return (jx_equ, jy_equ, ax_equ, ay_equ, az_equ, atan_equ)

#    for i in range(len(names)):
    if printThese!=[]:
        print hdr % ('Name','j_x (m-as/yr^3)','j_x (sigma)','j_y (m-as/yr^3)','j_y (sigma)','Position (mas)','Time')
        for i in range(len(printThese)):
            for j in range(len(names)):
                if str(names[j]) == str(printThese[i]):
                    print fmt % (names[j],jx[j]*1e6,sigx[j],jy[j]*1e6,sigy[j],x0[j]*1e3,y0[j]*1e3,t0y[j])
                    if str(names[j]) == str('S0-49'):
                        print x0[j],y0[j]
                        print vx[j],vy[j]
                        print ax[j],ay[j]
                        print jx[j],jy[j]
                        print GM_as_yr
 #                   pointsTab = asciidata.open(root + alnDir + points + str(names[j]) + '.points')
 #                   ptime = pointsTab[0].tonumpy()
 #                   px = pointsTab[1].tonumpy()
 #                   py = pointsTab[2].tonumpy()
 #                   pxe = pointsTab[3].tonumpy()
 #                   pye = pointsTab[4].tonumpy()
#        for i in range(len(printThese)):
#            for j in range(len(names)):
#                if str(names[j]) == str(printThese[i]):
#                    GM_mas_yr = GM_as_yr*1e9
#                    x_prime = x0[j]*1e3
#                    y_prime = y0[j]*1e3
#                    vx_prime = vx[j]*1e3
#                    vy_prime = vy[j]*1e3
#                    ax_prime = ax[j]*1e3
#                    ay_prime = ay[j]*1e3
#                    jx_prime = jx[j]*1e3
#                    jy_prime = jy[j]*1e3
#                    print names[j]
#                    solvedArray = opter.fsolve(jerkAccel_equ,[0,0,0,0,0,0])
#                    print 'BH coordinates:'
#                    print solvedArray[0],solvedArray[1],solvedArray[2]
#                    print 'Z position, velocity, accel:'
#                    print solvedArray[3],solvedArray[4],solvedArray[5]
#                    print ''
#                    print GM_as_yr
#                    print x_prime, y_prime,vx_prime,vy_prime,ax_prime,ay_prime,jx_prime,jy_prime
                    
    else:
#        pdb.det_trace()
        idx = np.where((abs(sigx) > 3.0) | (abs(sigy) > 3.0))[0]
        print hdr % ('Name','j_x (m-as/yr^3)','j_x (sigma)','j_y (m-as/yr^3)','j_y (sigma)','Position (mas)','Time')
        for j in idx:
            print fmt % (names[j],jx[j]*1e6,sigx[j],jy[j]*1e6,sigy[j],x0[j]*1e3,y0[j]*1e3,t0y[j])



def plotStarJerk(starName,alnDir='14_06_18/', root_tmp='/g/ghez/align/',
                 align = 'align/align_d_rms_1000_abs_t', poly='polyfit_nz/fit', points='points_nz/', polyj='polyfit_nzj/fit',
                 f_test=True, pvalue=4.0, chainsDir = 'efit/chains_S0-2_newRV2/',outTolatex=False,outToAL=False,
                 starlist='all', plotSigAcc=False, magCut=22, nEpochs=14, verbose = True):
    """
    Plots figures for a given star's jerk fit, if it has one
    """

    outdir = home + 'plots/'

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
    jxchi2 = _fit[9].tonumpy()

    jy0 = _fit[11].tonumpy()
    jvy = _fit[12].tonumpy()
    jay = _fit[13].tonumpy()
    jy = _fit[14].tonumpy()
    jy0e = _fit[15].tonumpy()
    jvye = _fit[16].tonumpy()
    jaye = _fit[17].tonumpy()
    jye = _fit[18].tonumpy()
    jychi2 = _fit[19].tonumpy()

    for i in range(len(jnames)):
        if (str(jnames[i]) == str(starName)):
            t0x = jt0x[i]
            x0 = jx0[i]
            vx = jvx[i]
            ax = jax[i]
            jx = jx[i]
            x0e = jx0e[i]
            vxe = jvxe[i]
            axe = jaxe[i]
            jxe = jxe[i]

            t0y = jt0y[i]
            y0 = jy0[i]
            vy = jvy[i]
            ay = jay[i]
            jy = jy[i]
            y0e = jy0e[i]
            vye = jvye[i]
            aye = jaye[i]
            jye = jye[i]

    pointsTab = asciidata.open(root_tmp + alnDir + points + starName + '.points')
    time = pointsTab[0].tonumpy()
    x = pointsTab[1].tonumpy()
    y = pointsTab[2].tonumpy()
    xerr = pointsTab[3].tonumpy()
    yerr = pointsTab[4].tonumpy()

    dt = time - t0x
    fitLineX = x0 + (vx * dt) + (ax * dt * dt / 2.0) + (jx * dt * dt * dt / 6.0)
    fitSigX = np.sqrt(x0e**2 + (vxe * dt)**2 + (axe * dt * dt / 2.0)**2 + (jxe * dt * dt * dt / 6.0)**2)

    fitLineY = y0 + (vy * dt) + (ay * dt * dt / 2.0) + (jy * dt * dt * dt / 6.0)
    fitSigY = np.sqrt(y0e**2 + (vye * dt)**2 + (aye * dt * dt / 2.0)**2 + (jye * dt * dt * dt / 6.0)**2)

    diffX = x - fitLineX
    diffY = y - fitLineY

    py.clf()
    dateTicLoc = MultipleLocator(3)
    dateTicRng = [1995, 2015]

    maxErr = np.array([xerr, yerr]).max()
    resTicRng = [-3*maxErr, 3*maxErr]

    from matplotlib.ticker import FormatStrFormatter
    fmtX = FormatStrFormatter('%5i')
    fmtY = FormatStrFormatter('%6.2f')

#    idtmp=np.where((time > 2006) & (abs(sigY) < 5))
    idtmp=np.where(time > 1990)
#    xyplim=[min(x[idtmp])-.002,max(x[idtmp])+.002,min(y[idtmp])-.002,max(y[idtmp])+.002]

    paxes = subplot(2,1,1)
    py.plot(time, fitLineX, 'b-')
    py.plot(time, fitLineX + fitSigX, 'b--')
    py.plot(time, fitLineX - fitSigX, 'b--')
    errorbar(time, x, yerr=xerr, fmt='k.')
    rng = axis()
    py.axis(dateTicRng + [rng[2], rng[3]])
    py.xlabel('Date (yrs)')
    py.ylabel('X (arcsec)')
    py.title(starName+' Jerk Fit')
#    title('X chi2_red = %4.2f' % fitx.chi2red)
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)
    
#    paxes = subplot(3, 2, 2)
    paxes = subplot(2,1,2)
    py.plot(time, fitLineY, 'b-')
    py.plot(time, fitLineY + fitSigY, 'b--')
    py.plot(time, fitLineY - fitSigY, 'b--')
    py.errorbar(time, y, yerr=yerr, fmt='k.')
    rng = axis()
    py.axis(dateTicRng + [rng[2], rng[3]])
    py.xlabel('Date (yrs)')
    py.ylabel('Y (arcsec)')
#    title('Y chi2_red = %4.2f' % fity.chi2red)
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)

    savefig(outdir+'plotJerkXY_t_' + starName + '.png')
    py.close(2)

    figure(3)
    clf()
    subplots_adjust(hspace=0.25,left=0.1, right=0.95,top=0.9, bottom=0.1, wspace=0.25)
    paxes=subplot(1,2,1)
    py.errorbar(time,diffX*1e3,yerr=xerr*1e3,fmt='b.')
    py.plot([0,3000],[0,0],'k')
    py.axis(dateTicRng+[-10,10])
    py.xlabel('Date (yrs)')
    py.ylabel('X Residuals (mas)')
    py.title('Velocity Fit')
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

    paxes=subplot(1,2,2)
    py.errorbar(time,diffY*1e3,yerr=yerr*1e3,fmt='b.')
    py.plot([0,3000],[0,0],'k')
    py.axis(dateTicRng+[-10,10])
    py.xlabel('Date (yrs)')
    py.ylabel('X Residuals (mas)')
#    py.title('Acceleration Fit')
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

#    paxes=subplot(2,2,3)
#    py.errorbar(time,diffYv*1e3,yerr=yerr*1e3,fmt='b.')
#    plot([0,3000],[0,0],'k')
#    axis(dateTicRng+[-10,10])
#    xlabel('Date (yrs)')
#    ylabel('Y Residuals (mas)')
#    paxes.xaxis.set_major_formatter(fmtX)
#    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

#    paxes=subplot(2,2,4)
#    errorbar(time,diffYa*1e3,yerr=yerr*1e3,fmt='b.')
#    plot([0,3000],[0,0],'k')
#    axis(dateTicRng+[-10,10])
#    xlabel('Date (yrs)')
#    ylabel('Y Residuals (mas)')
#    title('Acceleration Fit')
#    paxes.xaxis.set_major_formatter(fmtX)
#    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

    savefig(outdir+'plot_JerkResiduals_'+starName+'.png')
    py.close(3)

    py.figure(4)
    py.clf()
    py.plot(fitLineX[idtmp], fitLineY[idtmp], 'b-')
    py.errorbar(x[idtmp], y[idtmp],xerr=xerr[idtmp], yerr=yerr[idtmp], fmt='k.')
    py.title(starName+' Jerk Fit')
#    rng = axis()
#    axis(xyplim)
    py.xlabel('X (arcsec)')
    py.ylabel('Y (arcsec)')
#    title('Y chi2_red = %4.2f' % fity.chi2red)
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.xaxis.set_major_formatter(fmtX)
#    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)
    
#    paxes = subplot(3, 2, 3)
#    plot(time, np.zeros(len(time)), 'b-')
#    plot(time, fitSigX, 'b--')
#    plot(time, -fitSigX, 'b--')
#    errorbar(time, x - fitLineX, yerr=xerr, fmt='k.')
    #axis(dateTicRng + resTicRng)
#    axis([2005,2011,-0.001,0.001])
#    xlabel('Date (yrs)')
#    ylabel('X Residuals (arcsec)')
#    paxes.get_xaxis().set_major_locator(dateTicLoc)

#    paxes = subplot(3, 2, 4)
#    plot(time, np.zeros(len(time)), 'b-')
#    plot(time, fitSigY, 'b--')
#    plot(time, -fitSigY, 'b--')
#    errorbar(time, y - fitLineY, yerr=yerr, fmt='k.')
    #axis(dateTicRng + resTicRng)
#    axis([2005,2011,-0.001,0.001])
#    xlabel('Date (yrs)')
#    ylabel('Y Residuals (arcsec)')
#    paxes.get_xaxis().set_major_locator(dateTicLoc)

#    bins = np.arange(-7, 7, 1)
#    subplot(3, 2, 5)
#    (n, b, p) = hist(sigX, bins)
#    setp(p, 'facecolor', 'k')
#    axis([-5, 5, 0, 20])
#    xlabel('X Residuals (sigma)')
#    ylabel('Number of Epochs')

#    subplot(3, 2, 6)
#    (n, b, p) = hist(sigY, bins)
#    axis([-5, 5, 0, 20])
#    setp(p, 'facecolor', 'k')
#    xlabel('Y Residuals (sigma)')
#    ylabel('Number of Epochs')

#    subplots_adjust(wspace=0.4, hspace=0.4, right=0.95, top=0.95)
    savefig(outdir+'plotJerkY_X_' + starName + '.png')
    py.close(4)



def roughDensity(alnDir='14_06_18/', root_tmp='/g/ghez/align/',
                 align = 'align/align_d_rms_1000_abs_t', updateErr = True,poly='polyfit_nz/fit',
                 points='points_nz/', sigma=3.0, polyj='polyfit_nzj/fit',f_test=True, pvalue=4.0, sig3d=1.0, errCut = 0.1,
                 chainsDir = 'efit/chains_S0-2_newRV2/',binNum=3.0, slope=[-10.0,10.0,0.1],inter=[0.0,10.0,0.1],
                 starlist='all', magCut=15.5, nEpochs=14, radCut = 1.7, verbose = True, guess=[0.8,-2.5,0.5]):
    """
    Plots
    """
#    outdir = root + alnDir + 'plots/'
    outdir = home + 'plots/'

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
    xchi2r = s.getArray('xchi2r')
    ychi2r = s.getArray('ychi2r')
        
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
    jxchi2 = _fit[9].tonumpy()
#    jxq = _fit[10].tonumpy()

    jy0 = _fit[11].tonumpy()
    jvy = _fit[12].tonumpy()
    jay = _fit[13].tonumpy()
    jy = _fit[14].tonumpy()
    jy0e = _fit[15].tonumpy()
    jvye = _fit[16].tonumpy()
    jaye = _fit[17].tonumpy()
    jye = _fit[18].tonumpy()
    jychi2 = _fit[19].tonumpy()
#    jyq = _fit[20].tonumpy()



    ###########
    #
    # F Test -- Accel or Velocity fit?
    #
    ###########
    if f_test == True:
        # Pass in the star names and run the F test
        pass_f,xprob,yprob = run_f_test(names,pvalue,root_tmp,alnDir,align,poly,points,verbose=verbose)
        pass_fj = run_fjerk_test(names,pvalue,root_tmp,alnDir,align,poly,points,verbose=verbose)
    else:
        pass_f = np.ones(len(names)) # check everything
        pass_fj = np.ones(len(names))

    idex = np.where(pass_fj ==1)[0]

    for i in idex:
        if (str(names[i]) == 'S0-16'):
            pass_f[i] == 1

    # S0-16 has such a high inclination that it is being seen as a velocity, force the f test vel/acc to be acc
    r = np.hypot(x0, y0)
    idex = np.where((pass_fj ==1) & (pass_f == 1) & (cnt > nEpochs) & (mag < magCut) & (r < radCut))[0]

    for i in idex:
#        pdb.set_trace()
        jdex = np.where(jnames == str(names[i]))
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
        xchi2r[i] = jxchi2[jdex]/(cnt[i] - 4.0)
        ychi2r[i] = jychi2[jdex]/(cnt[i] - 4.0)
        t0x[i] = jt0x[jdex]
        t0y[i] = jt0y[jdex]
#        ar[i] = jar[jdex]
#        at[i] = jat[jdex]
#        are[i] = jare[jdex]
#        ate[i] = jate[jdex]
#    sigmaR = ar/are
#    sigmaT = at/ate


    # Get accelerations in the radial/tangential direction
    r = np.hypot(x0, y0)

    # Make an epochs cut
    idx = np.where((cnt > nEpochs) & (mag < magCut) & (r < radCut))[0]
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
    xchi2r = xchi2r[idx]
    ychi2r = ychi2r[idx]
    r = r[idx]
    t0x = t0x[idx]
    t0y = t0y[idx]
    names = [names[nn] for nn in idx]

    print 'Number in sample cut:'
    print len(idx)

    ar = ((ax*x0) + (ay*y0)) / r
    at = ((ax*y0) - (ay*x0)) / r
    are =  (axe*x0/r)**2 + (aye*y0/r)**2
    are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    are =  np.sqrt(are)
    ate =  (axe*y0/r)**2 + (aye*x0/r)**2
    ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
    ate =  np.sqrt(ate)

#    writeJLatex = np.where((pass_fj ==1) & (pass_f == 1) & (r < 1.0))[0]
#    pdb.set_trace()

    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, alnDir=alnDir, chainsDir = chainsDir)
#        jx0e,jy0e,jvxe,jvye = nzErr(jx0e, jy0e, jvxe, jvye, jt0x, jt0y, jmag, alnDir=alnDir, chainsDir = chainsDir)

    if updateErr:
        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17.1])
        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5
        

        #Finding the delta to add in quad to the error depending on the radius and mag, whichever one, radius or mag
        #error is larger
        delta = mag*0.0

        for i in range(len(mag)):
            for j in range(len(atBins)-1):
                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
                    delta[i] = deltaArr[j]            

        ate = np.sqrt(ate**2 + delta**2)
        are = np.sqrt(are**2 + delta**2)

    r3dMax =  radCut * dist / au_in_pc

    x0_pc = x0 * dist / au_in_pc
    y0_pc = y0 * dist / au_in_pc
    x0_cm = x0_pc * cm_in_pc
    y0_cm = y0_pc * cm_in_pc
    x0e_pc = x0e * dist / au_in_pc
    y0e_pc = y0e * dist / au_in_pc
    x0e_cm = x0e_pc * cm_in_pc
    y0e_cm = y0e_pc * cm_in_pc
#    r2d = np.sqrt(x0**2 + y0**2) # arcsec
#    r2d_pc = r2d * dist / au_in_pc
    r_cm = r * dist * cm_in_pc / au_in_pc
    re_cm = np.sqrt((x0_cm*x0e_cm/r_cm)**2 + (y0_cm*y0e_cm/r_cm)**2)
    ar_cmss = ar * asy_to_kms * 1e5 / sec_in_yr
    are_cmss = are * asy_to_kms * 1e5 / sec_in_yr
#    vx_cms = vx * asy_to_kms * 1e5
#    vy_cms = vy * asy_to_kms * 1e5
#    vxe_cms = vxe * asy_to_kms * 1e5
#    vye_cms = vye * asy_to_kms * 1e5
#    vproj_cms = np.sqrt(vx_cms**2 + vy_cms**2)
#    vproje_cms = np.sqrt((vx_cms*vxe_cms)**2 + (vy_cms*vye_cms)**2)/vproj_cms
    z0_cm = np.sqrt(abs(((GM * r_cm / ar_cmss)**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)

#    az_cmss = GM * z0_cm / ((x0_cm**2 + y0_cm**2 + z0_cm**2)**(1.5))

    z0e_cm = (abs(G*r_cm/ar_cmss))**(4.0/3.0)*(mass_g)**(-2.0/3.0)*masse_g**2 + ((abs(GM/ar_cmss))**(2.0/3.0)*(abs(r_cm))**(-1.0/3.0)-r_cm)**2*re_cm**2
    z0e_cm += (abs(GM*r_cm))**(4.0/3.0)*(abs(ar_cmss))**(-10.0/3.0)*are_cmss**2
    z0e_cm = np.sqrt(z0e_cm)/(3.0 * z0_cm)

    r3d_cm = np.sqrt(x0_cm**2 + y0_cm**2 + z0_cm**2)
    r3de_cm = np.sqrt(((x0_cm*x0e_cm)**2 + (y0_cm*y0e_cm)**2 + (z0_cm*z0e_cm)**2)/r3d_cm**2)
    r3d_pc = r3d_cm / cm_in_pc
    r3de_pc = r3de_cm / cm_in_pc

#    vz_cms = (GM/r3d_cm) - vproj_cms**2

    #say that if a star is unbound or within 1sig of being unbound, set projected radius as its 3D radius

    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
#    cur.execute('SELECT * FROM unknownSims')
    oldProb = np.zeros(len(names))

    print "Stars that are unbound within 1 sig:"
    for i in range(len(names)):
#        tmpName = str(names[i])
        if ((((GM * (r_cm[i]+re_cm[i]) / ar_cmss[i])**2)**(1.0/3.0) - r_cm[i]**2) <= 0):
            print names[i]
            r3d_cm[i] = r_cm[i]
            r3de_cm[i] = re_cm[i]
            r3d_pc[i] = r_cm[i] / cm_in_pc
            r3de_pc[i] = re_cm[i] / cm_in_pc

#    print 'Stars that are not in database at all:'
    for i in range(len(names)):
        tmpName = str(names[i])
#        if (tmpName == 'S1-38'):
#            pdb.set_trace()
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
#                print tmpName

    print 'Star, Prob Old, r in arcsec:'
    for i in range(len(names)):
        print names[i],oldProb[i],r3d_pc[i]*au_in_pc / dist
    #pdb.set_trace()
    #Have the probs for all stars in decided sample, determine bins

    sig_r3d = np.abs(r3d_pc / r3de_pc)

    fmdex = np.where((sig_r3d >= sig3d) & (r3de_pc < errCut) & (r3d_pc <= r3dMax) & (oldProb > 0.0))[0]
#    maxR = np.max(np.log10(r3d_pc[fmdex]))+0.000001
#    minR = np.min(np.log10(r3d_pc[fmdex]))-0.000001

    maxR = np.max(r3d_pc[fmdex])+0.000001
    minR = np.min(r3d_pc[fmdex])-0.000001
    binSize = (maxR-minR)/binNum
    radLim = np.array([minR+binSize*i for i in range(binNum+1)])
#    radLim = np.array([binSize*i for i in range(binNum+1)])
    number = np.zeros(binNum)
    volume = np.zeros(binNum)
    radius3D = np.zeros(binNum)

    for i in range(binNum):
#        pdb.set_trace()
#        tmpDex = np.where((np.log10(r3d_pc) >= radLim[i]) & (np.log10(r3d_pc) < radLim[i+1]) & (sig_r3d >= sig3d) & (r3de_pc < errCut))[0]
        tmpDex = np.where((oldProb > 0.0) & (r3d_pc >= radLim[i]) & (r3d_pc < radLim[i+1]) & (sig_r3d >= sig3d) & (r3de_pc < errCut))[0]
        tmpSum = oldProb[tmpDex]
        number[i] = tmpSum.sum()
#        volume[i] = 4.0*pi*((10.0**radLim[i+1])**3 - (10.0**radLim[i])**3)/3.0
        volume[i] = 4.0*pi*(radLim[i+1]**3 - radLim[i]**3)/3.0
        radius3D[i] = np.median((r3d_pc[tmpDex]))
        pdb.set_trace()

    densityBin = number/volume #number density in pc^-3
    numerr = np.sqrt(number)
    err_pl = np.log10((numerr+number)/number)
    error_density = numerr/volume

#    pdb.set_trace()
    sval,ival=scLeast2(radius3D,densityBin,error=error_density,slope=slope,inter=inter,guess=guess)

    print number

    plOld = np.where((oldProb > 0.0) & (sig_r3d >= sig3d) & (r3de_pc < errCut) & (r3d_pc <= r3dMax))[0]
    ploz = np.zeros(len(plOld))+10.0**0.05
    
    rho_pl = np.log10(densityBin)
    r_pl = np.log10(radius3D)
    rLeft = np.array([radLim[i] for i in range(binNum)])
    errpatch = np.array([[densityBin[0]-1.0,error_density[1],error_density[2]],error_density])
 #   fitLine = r_pl*sval + ival
    py.clf()
    py.bar(rLeft,densityBin,width=binSize,bottom=1.0,color='w')
    py.errorbar(radius3D,densityBin,yerr=errpatch,fmt='o')
    py.plot((r3d_pc[plOld]),ploz,'D')
 #   py.plot(r_pl,fitLine)
    py.xscale('log')
    py.yscale('log')
    py.xlabel('Log Radius (pc)')
    py.ylabel('Log Density (pc^-3)')
    py.savefig(outdir + 'roughRho.png')

    pdb.set_trace()

    print "log radius and log density"
    print r_pl
    print rho_pl
    print err_pl
    


def protolikeGamma(numz=1000.,rrange=[1.8,100.0,1000.],dataCut = 100.0,projCut=1.7,grange=[-5.0,2.0,1000.],alnDir='14_06_18/', root_tmp='/g/ghez/align/',
                   align = 'align/align_d_rms_1000_abs_t', updateErr = True,poly='polyfit_nz/fit',
                   points='points_nz/', polyj='polyfit_nzj/fit', outTolatex=False,outToAL=False,
                   starlist='all',magCut=15.5,nEpochs=14,likelihood=True):
    #Read in 
    #run histaccel every time or write out to .txt and read from there?    

    names,Rproj_a,ar_a,are_a=histAccel(alnDir=alnDir,root_tmp=root_tmp,align = align,radCut = projCut,poly=poly, points=points,polyj=polyj,
                                       outTolatex=outTolatex,outToAL=outToAL,starlist=starlist, magCut=magCut, nEpochs=nEpochs, likelihood=likelihood)

    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
#    cur.execute('SELECT * FROM unknownSims')
    oldProb = np.zeros(len(names))

#    print "Stars that are unbound within 1 sig:"
#    for i in range(len(names)):
#        tmpName = str(names[i])
#        if ((((GM * (r_cm[i]+re_cm[i]) / ar_cmss[i])**2)**(1.0/3.0) - r_cm[i]**2) <= 0):
#            print names[i]
#            r3d_cm[i] = r_cm[i]
#            r3de_cm[i] = re_cm[i]
#            r3d_pc[i] = r_cm[i] / cm_in_pc
#            r3de_pc[i] = re_cm[i] / cm_in_pc

#    print 'Stars that are not in database at all:'
    for i in range(len(names)):
        tmpName = str(names[i])
#        if (tmpName == 'S1-38'):
#            pdb.set_trace()
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

    drad = (rrange[1]-rrange[0])/rrange[2]
    radCut = np.array([rrange[0] + i*drad for i in range(rrange[2])])
    rcut = radCut * dist * cm_in_au
    Rcut = projCut * dist * cm_in_au
    rDataCut = dataCut * dist * cm_in_au
#    rmin = radMin * dist * cm_in_au

    #test
#    Rproj = rcut/2.0
#    ar = GM*Rproj/(np.sqrt(Rproj**2 + (Rproj/2.0)**2))**3 
#    are = ar/25.0

    dgam = (grange[1]-grange[0])/grange[2]
    gamma = np.array([grange[0] + i*dgam for i in range(grange[2])])
    no2 = np.where(gamma != 2.0)
    gamma = gamma[no2]
    grange[2] = len(gamma)

    lnL = np.zeros([rrange[2],grange[2]])
    print ""
    print 'Stars in likelihood and probOld'

    for j in range(len(Rproj_a)):
        if ((oldProb[j] > 0.0)): #& (ar_a[j] < -1.0*are_a[j])):
            print names[j],oldProb[j]

#            pdb.set_trace()
            Rproj = Rproj_a[j]
            ar = ar_a[j]
            are = are_a[j]
            print ar,are
            prob_star = np.zeros([rrange[2],grange[2]])

            for k in range(rrange[2]):
#                pdb.set_trace()
                rcut_tmp = rcut[k]
                
                zcut = np.sqrt(rcut_tmp**2 - Rproj**2)
                dz = zcut/numz
#                z = np.array([0.0 + i*(zcut/numz) for i in range(numz)])
                z = np.array([(dz/100.0) + i*dz for i in range(numz)])
                
                if (ar < 0.0):
                    a_p = -1.0*GM*Rproj/(np.sqrt(Rproj**2 + z**2))**3
                    prob_ap = np.exp(-1.0*(ar - a_p)**2/(2.0*are**2))
#                norm1 = np.sum(prob_ap) #*dz
                    norm1 = abs(are*math.sqrt(pi/2.0)*(special.erf((ar - min(a_p))/(math.sqrt(2)*are)) - special.erf((ar - max(a_p))/(math.sqrt(2)*are))))
#                    norm1 = abs(are*math.sqrt(pi/2.0)*(special.erf((ar - min(a_p))/(math.sqrt(2)*are)) - special.erf((ar/(math.sqrt(2)*are)))))
                    if (norm1 == 0.0):
                        norm1 = integrate.quad(lambda x: math.exp(-1.0*(ar - x)**2/(2.0*are**2)), min(a_p), max(a_p))
#                        norm1 = integrate.quad(lambda x: math.exp(-1.0*(ar - x)**2/(2.0*are**2)), min(a_p), 0.0)
                        norm1 = abs(norm1[0])
                else:
                    prob_ap = (Rproj**2 + z**2)**(2.5)/(3.0*z*GM*Rproj)
#                    norm1 = np.max(z)
                    norm1 = math.sqrt(rDataCut**2 - Rproj**2)
                prob_a = prob_ap/norm1
                

        #   pdb.set_trace()
#            rho_irho = np.zeros([numz,grange[2]])
#             for i in range(numz):
#                z_tmp = z[i]
#                rho_tmp = (gamma+3.0)*(np.sqrt(Rproj**2 + z_tmp**2))**gamma/(rcut**(gamma+3.0)-rmin**(gamma+3.0))
#                norm2 = np.sum(rho_tmp)*dgam
#                rho_irho[i] = np.abs(rho_tmp/norm2)

            #Sum over z
#            prob_gam = np.array([np.sum(prob_a * rho_irho[:,i]) * dz for i in range(grange[2])])
#            pdb.set_trace()

                for i in range(grange[2]):
#                    pdb.set_trace()
                    g_tmp = gamma[i]
                    rho_tmp = (np.sqrt(Rproj**2 + z**2))**(-1.0*g_tmp)
#                    norm2 = -2.0*zcut*Rcut**2*((zcut/Rcut)**2+1.0)**(g_tmp/2.0)*(zcut**2+Rcut**2)**(g_tmp/-2.0)*special.hyp2f1(0.5,((g_tmp/2.0)-1.0),1.5,-1.0*(zcut/Rcut)**2)/(g_tmp-2.0)
                    norm2a = integrate.quad(lambda x: x**(2.0-g_tmp), 0.0, rcut_tmp)
                    norm2b = integrate.quad(lambda x: x**2.0 * (Rcut**2.0 + x**2)**(g_tmp/-2.0), 0.0, math.sqrt(rcut_tmp**2.0 - Rcut**2.0))
                    norm2 = norm2a[0] - norm2b[0]
#                    test = integrate.dblquad(lambda zprime,Rprime: 2.0*Rprime*(zprime**2+Rprime**2)**(g_tmp/-2.0), 0, Rcut, lambda zprime: 0, 
#                                             lambda zprime: math.sqrt(rcut_tmp**2.0 - zprime**2.0))
#                    test = test[0]
                    prob_z = rho_tmp/norm2

                    prob_rcut_gam = prob_a*prob_z*dz
                    try:
                        lnL[k,i] += oldProb[j] * math.log(np.sum(prob_rcut_gam))
                        prob_star[k,i] = np.sum(prob_rcut_gam)
                    except:
                        pdb.set_trace()
#                if (((str(names[j]) == str('S1-13')) | (str(names[j]) == str('S0-35')) | (str(names[j]) == str('S1-17'))) &
#                    (i==0)):
#                    pdb.set_trace()

#            lnL = lnL + oldProb[j] * np.log(prob_gam)
#            pdb.set_trace()
    print lnL
    pdb.set_trace()






def likeGamma_nw(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
              poly='polyfit_nz/fit',points='points_nz/',R2Dcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,dataCut=100.,
              starlist='all',chainsDir='efit/chains_S0-2_newRV2/',grange=[-14.0,2.0,100.],rrange=[1.71,10.0,100.],numz=100.):

    #Load info on catalogue, x and y positions temp, they are assigned by which t0 is used for that star
    #Won't be using that t0, will define global t0 for all stars
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)

    names = s.getArray('name')
    xtmp = s.getArray('x0')
    ytmp = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

    projRtmp = np.sqrt(xtmp**2 + ytmp**2)

    cut1 = np.where((projRtmp < R2Dcut) & (cnt > nEpochs) & (mag > magCut))[0]
    names = [names[nn] for nn in cut1]
    projRtmp = projRtmp[cut1] * dist / au_in_pc

    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()
    oldProb = np.zeros(len(names))

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

    #Making plot of 2D distribution of old stars only
    py.clf()
    binNum = 3.0
    plOld = np.where(oldProb > 0.0)[0]
    ploz = np.zeros(len(plOld))+10.0**0.05
    maxR = np.max(projRtmp[plOld])+0.00001
    minR = np.min(projRtmp[plOld])
    binSize = (maxR - minR)/binNum
    radLim = np.array([minR+binSize*i for i in range(binNum + 1)])
    number = np.zeros(binNum)
    volume = np.zeros(binNum)
    radius2D = np.zeros(binNum)

#    pdb.set_trace()
    for i in range(binNum):
        tmpDex = np.where((oldProb > 0.0) & (projRtmp >= radLim[i]) & (projRtmp < radLim[i+1]))[0]
        tmpSum = oldProb[tmpDex]
        number[i] = tmpSum.sum()
        volume[i] = pi*(radLim[i+1]**2.0 - radLim[i]**2.0)
        radius2D[i] = np.median(projRtmp[tmpDex])
    densityBin = number/volume
    numerr = np.sqrt(number)
    err_pl = np.log10((numerr+number)/number)
    error_density = numerr/volume
    rho_pl = np.log10(densityBin)
    r_pl = np.log10(radius2D)
    rLeft = np.array([radLim[i] for i in range(binNum)])
    py.bar(rLeft,densityBin,width=binSize,bottom=1.0,color='w')
    py.errorbar(radius2D,densityBin,yerr=error_density,fmt='o')
    py.plot(projRtmp[plOld],ploz,'D')
    py.xscale('log')
    py.yscale('log')
    py.xlabel('Log Projected Radius (pc)')
    py.ylabel('Log Surface Density (pc$^{-2}$)')
    py.savefig('rough2D.png')
    pdb.set_trace()

    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
        ori_x0e = origin_val[18][0]
        ori_y0e = origin_val[19][0]
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]
        t_0 = 2000.0 #hard coded t_0 of sgr*

    drad = (rrange[1]-rrange[0])/rrange[2]
    radCut = np.array([rrange[0] + i*drad for i in range(rrange[2])])
    rcut = radCut * dist * cm_in_au
    dataCut = dataCut * dist * cm_in_au
    R2Dcut = R2Dcut * dist * cm_in_au

    dgam = (grange[1]-grange[0])/grange[2]
    gamma = np.array([grange[0] + i*dgam for i in range(grange[2])])
    no2 = np.where(gamma != 2.0)
    gamma = gamma[no2]
    grange[2] = len(gamma)

    lnL = np.zeros([grange[2],rrange[2]])

    for i in range(len(names)):
        if (oldProb[i] > 0.0):
            tmpName = names[i]
            print tmpName,oldProb[i]

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
            sig2x = np.sum(xerr**2.0)
            sig2y = np.sum(yerr**2.0)
            ws_tx = sig2x * np.sum((time - globalt0)/xerr**2.0)
            ws_t2x = sig2x * np.sum((time - globalt0)**2.0/xerr**2.0)
            ws_ty = sig2y * np.sum((time - globalt0)/yerr**2.0)
            ws_t2y = sig2y * np.sum((time - globalt0)**2.0/yerr**2.0)

            ws_xi = sig2x * np.sum(xp / xerr**2.0)
            ws_yi = sig2y * np.sum(yp / yerr**2.0)
            ws_xi_t = sig2x * np.sum(xp * (time - globalt0) / xerr**2.0)
            ws_yi_t = sig2y * np.sum(yp * (time - globalt0) / yerr**2.0)

            for k in range(rrange[2]):
                rmax = rcut[k]

                x_tilda = (ws_xi * ws_t2x - ws_tx * ws_xi_t) / (ws_t2x - ws_tx**2.0)
                y_tilda = (ws_yi * ws_t2y - ws_tx * ws_yi_t) / (ws_t2y - ws_ty**2.0)
                ratio_tilda = y_tilda/x_tilda
                zmax = math.sqrt(rmax**2.0 - ((x_tilda**2.0 + y_tilda**2.0)/(1.0 - ws_tx * ws_ty * GM / (2.0 * rmax**3.0))**2.0))
                dz = zmax / numz
                zp = np.array([h * dz for h in range(numz)])

                #rmax = rcut[k]

                #rmin = np.min(np.sqrt(xp**2 + yp**2))
                #rtmp = np.array([rmin + h*((rmax-rmin)/numz) for h in range(numz)])
                #accel = -1.0*GM/rtmp**3.0

#                dz = math.sqrt(rmax**2.0 - np.min(xp**2.0 + yp**2.0))/numz
#                ztmp = np.array([h*dz for h in range(numz)])

                sol_x0 = np.zeros(numz)
                sol_y0 = np.zeros(numz)

                ws_x = np.zeros(numz)
                ws_x2 = np.zeros(numz)
                ws_x_t = np.zeros(numz)
                ws_y = np.zeros(numz)
                ws_y2 = np.zeros(numz)
                ws_y_t = np.zeros(numz)

                for q in range(numz):
                    ztmp = zp[q]
                    xfunc = lambda chi: chi*(1.0 - ws_tx**2.0*GM/(2.0*(chi**2.0 + ratio_tilda**2.0*chi**2.0 + ztmp**2.0)**(3.0/2.0))) - x_tilda
                    sol_x0[q] = fsolve(xfunc, x_tilda)
                    sol_y0[q] = sol_x0[q] * ratio_tilda
                    atmp = -1.0*GM/(ztmp**2.0 * sol_x0[q]**2.0 + sol_y0[q]**2.0)**(3.0/2.0)
                    ws_x[q] = sig2x * np.sum((xp + 0.5*(time - globalt0)**2.0*sol_x0[q]*atmp)/xerr**2.0)
                    ws_x2[q] = sig2x * np.sum((xp + 0.5*(time - globalt0)**2.0*sol_x0[q]*atmp)**2.0/xerr**2.0)
                    ws_x_t[q] = sig2x * np.sum((time - globalt0)*(xp + 0.5*(time - globalt0)**2.0*sol_x0[q]*atmp)/xerr**2.0)
                    ws_y[q] = sig2y * np.sum((yp + 0.5*(time - globalt0)**2.0*sol_y0[q]*atmp)/yerr**2.0)
                    ws_y2[q] = sig2y * np.sum((yp + 0.5*(time - globalt0)**2.0*sol_y0[q]*atmp)**2.0/yerr**2.0)
                    ws_y_t[q] = sig2y * np.sum((time - globalt0)*(yp + 0.5*(time - globalt0)**2.0*sol_y0[q]*atmp)/yerr**2.0)

#                    atmp = -1.0*GM/(np.sqrt((ztmp[q])**2.0 + xp**2.0 + yp**2.0))**3.0
#                    ws_x[q] = sig2x * np.sum(xp*(1.0 + 0.5*(time - globalt0)**2.0*atmp)/xerr**2.0)
#                    ws_x2[q] = sig2x * np.sum(xp**2.0*(1.0 + 0.5*(time - globalt0)**2.0*atmp)**2.0/xerr**2.0)
#                    ws_x_t[q] = sig2x * np.sum((time - globalt0)*xp*(1.0 + 0.5*(time - globalt0)**2.0*atmp)/xerr**2.0)
#                    ws_y[q] = sig2y * np.sum(yp*(1.0 + 0.5*(time - globalt0)**2.0*atmp)/yerr**2.0)
#                    ws_y2[q] = sig2y * np.sum(yp**2.0*(1.0 + 0.5*(time - globalt0)**2.0*atmp)**2.0/yerr**2.0)
#                    ws_y_t[q] = sig2y * np.sum((time - globalt0)*yp*(1.0 + 0.5*(time - globalt0)**2.0*atmp)/yerr**2.0)



                likex = (ws_x2 - ws_x**2.0 - ((ws_x_t - ws_x*ws_tx)**2.0/(ws_t2x - ws_tx**2.0)))/(-2.0*sig2x)
                likey = (ws_y2 - ws_y**2.0 - ((ws_y_t - ws_y*ws_ty)**2.0/(ws_t2y - ws_ty**2.0)))/(-2.0*sig2y)
                bandaid = np.max(likex + likey)
                prob_a = np.exp(likex + likey - bandaid)

#                ml_x = ((ws_x * ws_t2x) - (ws_tx * ws_x_t))/(ws_t2x - ws_tx**2.0)
#                ml_y = ((ws_y * ws_t2y) - (ws_ty * ws_y_t))/(ws_t2y - ws_ty**2.0)

                #r3d_tmp = np.sqrt(ml_x**2.0 + ml_y**2.0 + ztmp**2.0)
                r3d_tmp = np.sqrt(sol_x0**2.0 + sol_y0**2.0 + zp**2.0)
                

                for j in range(grange[2]):
                    gtmp = gamma[j]
                    rho_tmp = r3d_tmp**(-1.0*gtmp)
                    norma = integrate.quad(lambda x: x**(2.0-gtmp), 0.0, rmax)
                    normb = integrate.quad(lambda x: x**2.0 * (R2Dcut**2.0 + x**2.0)**(gtmp/-2.0), 0.0, 
                                           math.sqrt(rmax**2.0 - R2Dcut**2.0))
                    norm = norma[0] - normb[0]
                    prob_z = rho_tmp/norm

                    bandaid2 = np.max(np.abs(prob_z))
                    prob_zp = prob_z/bandaid2

                    prob_rcut_gam = prob_a*prob_z*dz

                    #pdb.set_trace()
                    try:
                        lnL[j,k] += oldProb[i] * (math.log(np.sum(prob_rcut_gam)) + bandaid + math.log(bandaid2))
                    except:
                        pdb.set_trace()

    pdb.set_trace()



def likeGamma(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
              poly='polyfit_nz/fit',points='points_nz/',nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,
              starlist='all',chainsDir='efit/chains_S0-2_newRV2/',Rcut = 1.7,grange=[-2.0,1.9,10.],
              maxZ=100.0,maxR=2.0,GLpoly = 10.0):

    #Load info on catalogue, x and y positions temp, they are assigned by which t0 is used for that star
    #Won't be using that t0, will define global t0 for all stars
    s = loadPop(root_tmp=root_tmp,alnDir=alnDir,starlist=starlist,align=align,poly=poly,points=points)

    names = s.getArray('name')
    xff = s.getArray('x0')
    yff = s.getArray('y0')
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    cut1 = np.where((mag < magCut) & (cnt > nEpochs) & ((xff**2.0 + yff**2.0) <= Rcut**2.0))[0]
    names = [names[nn] for nn in cut1]

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

    dgam = (grange[1]-grange[0])/(grange[2] - 1.0)
    gamma = np.array([grange[0] + i*dgam for i in range(grange[2])])
    no2 = np.where(gamma != 2.0)
    gamma = gamma[no2]
    grange[2] = len(gamma)
    lnL = np.zeros(grange[2])
    norm = np.zeros(grange[2])
#    normTest = np.zeros(grange[2])

    maxZ = maxZ * dist * cm_in_au
    maxR = maxR * dist * cm_in_au

    weight_gl,gl,check_gl = GaussLegendreWeights(GLpoly)
    if check_gl != 0:
        print 'GAUSS LEGENDRE PROBLEM'
        print 'FULL STOP'
        stop

    x_gl = gl * maxR
    y_gl = gl * maxR
    z_gl = gl * maxZ/2.0 + maxZ/2.0
    #pdb.set_trace()

#    dz = maxZ / (numz - 1.0)
#    z = np.array([0.0 + i*dz for i in range(numz)])
#    dxy = 2.0*maxR / (numxy - 1.0)
#    x = np.array([-1.0*maxR + i*dxy for i in range(numxy)])
#    y = np.array([-1.0*maxR + i*dxy for i in range(numxy)])

#    rMAX = math.sqrt(maxZ**2.0 + maxR**2.0)
#    start_time = time.time()
#    for i in range(grange[2]):
#        norma = integrate.quad(lambda x: x**(2.0-gamma[i]),0.0,rMAX)
#        normb = integrate.quad(lambda x: x**2.0 * (maxR**2.0 + x**2.0)**(gamma[i]/-2.0), 0.0, maxZ)
#        norm[i] = norma[0] - normb[0]

#    print 'Running time for 1D in minutes:'
#    print (time.time() - start_time)/60.0

#    start_time = time.time()
    for i in range(grange[2]):
        #normblah = integrate.dblquad(lambda zprime,Rprime: Rprime*(zprime**2.0 + Rprime**2.0)**(gamma[i]/-2.0), 0.0, maxR,
                                  #lambda Rprime: 0.0, lambda Rprime: maxZ)
        normt = integrate.tplquad(lambda zprime,yprime,xprime: (xprime**2.0 + yprime**2.0 + zprime**2.0)**(gamma[i]/-2.0),
                                  -1.0*maxR, maxR, lambda xprime: -1.0*maxR, lambda xprime: maxR, lambda xprime, yprime: 
                                  0.0, lambda xprime, yprime: maxZ)
        norm[i] = normt[0]
        #pdb.set_trace()
#    print 'Running time for 2D in minutes:'
#    print (time.time() - start_time)/60.0

#    pdb.set_trace()

    print ''
    print ''
    print ''
    print 'Stars and prob of being old:'

    for i in range(len(names)):
        #Cycle through stars and consider those w/ non-zero prob of being late-type
        if (oldProb[i] > 0.0):
            tmpName = names[i]
            print tmpName,oldProb[i]

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

            #Cycle through x0, y0, and z
            #accel_lnL = np.array([])
            #rho_prob = np.zeros([numxy,numxy,numz,grange[2]])
            starlike = np.zeros(grange[2])
            for n in range(grange[2]):
                gtmp = gamma[n]
                accel_lnL = np.array([])
                rho_prob = np.array([])
                for j in range(len(x_gl)):
                    xtmp = x_gl[j]
                    xweight = weight_gl[j]
                    for m in range(len(z_gl)):
                        ztmp = z_gl[m]
                        zweight = weight_gl[m]
                        for k in range(len(y_gl)):
                            #if (math.sqrt(xtmp**2.0 + y[k]**2.0) < maxR):
                            ytmp = y_gl[k]
                            yweight = weight_gl[k]
                            atmp = -0.5*GM / (xtmp**2.0 + ytmp**2.0 + ztmp**2.0)**(3.0/2.0)

                            xtilda = (xp - xtmp - xtmp*atmp*(time - globalt0)**2.0)
                            ytilda = (yp - ytmp - ytmp*atmp*(time - globalt0)**2.0)

                            axarray = np.zeros(2)
                            ayarray = np.zeros(2)
                            axarray[0] = np.sum(xtilda*(time - globalt0)/xerr**2.0)
                            axarray[1] = np.sum(xtilda*(time - globalt0)**3.0/xerr**2.0)
                            ayarray[0] = np.sum(ytilda*(time - globalt0)/yerr**2.0)
                            ayarray[1] = np.sum(ytilda*(time - globalt0)**3.0/yerr**2.0)

                            x_af1 = np.dot(axarray,inv_xf)
                            y_af1 = np.dot(ayarray,inv_yf)

                            accel_tmp = (np.sum(xtilda**2.0/xerr**2.0) - np.dot(x_af1,axarray) + 
                                         np.sum(ytilda**2.0/yerr**2.0) - np.dot(y_af1,ayarray))/-2.0

                            #for n in range(grange[2]):
                                #gtmp = gamma[n]
                            rho_tmp = xweight*yweight*zweight*(xtmp**2.0 + ytmp**2.0 + ztmp**2.0)**(gtmp/-2.0) / norm[n]
                                #starlike[n] += math.exp(accel_lnL[j,k,m] - star_constant) * rho_prob
                            accel_lnL = np.append(accel_lnL, accel_tmp)
                            rho_prob = np.append(rho_prob, rho_tmp)

                            #cycle through gamma, add to likelihood for individual star
                            #for n in range(grange[2]):
                                #gtmp = gamma[n]
                                #rho = (xtmp**2.0 + ytmp**2.0 + ztmp**2.0)**(gtmp/-2.0)
                                #rho_prob[j,k,m,n] = rho/norm[n]
                            #Can't have 4D array with number of points needed

            #star_constant = np.max(accel_lnL)
            #Cycle through again for gamma
            #for j in range(len(x)):
                #xtmp = x[j]
                #for m in range(len(z)):
                    #ztmp = z[m]
                    #for k in range(len(y)):
                        #if (math.sqrt(xtmp**2.0 + y[k]**2.0) < maxR):
                            #ytmp = y[k]
                            #for n in range(grange[2]):
                                #gtmp = gamma[n]
                                #rho_prob = (xtmp**2.0 + ytmp**2.0 + ztmp**2.0)**(gtmp/-2.0) / norm[n]
                                #starlike[n] += math.exp(accel_lnL[j,k,m] - star_constant) * rho_prob

                if (n == 0):
                    star_constant = np.max(accel_lnL)
                starlike[n] = np.sum(np.exp(accel_lnL - star_constant)*rho_prob)

            #starlike = np.array([np.sum(np.exp(accel_lnL - star_constant)*rho_prob[:,:,:,n]) for n in range(grange[2])])
            lnL += oldProb[i]*(np.log(starlike) + star_constant)
    pdb.set_trace()


##################################################################
# Recursive generation of the Legendre polynomial of order n
def Legendre(n,x):
	x=array(x)
	if (n==0):
		return x*0+1.0
	elif (n==1):
		return x
	else:
		return ((2.0*n-1.0)*x*Legendre(n-1,x)-(n-1)*Legendre(n-2,x))/n
 
##################################################################
# Derivative of the Legendre polynomials
def DLegendre(n,x):
	x=array(x)
	if (n==0):
		return x*0
	elif (n==1):
		return x*0+1.0
	else:
		return (n/(x**2-1.0))*(x*Legendre(n,x)-Legendre(n-1,x))
##################################################################
# Roots of the polynomial obtained using Newton-Raphson method
def LegendreRoots(polyorder,tolerance=1e-20):
	if polyorder<2:
		err=1 # bad polyorder no roots can be found
	else:
		roots=[]
		# The polynomials are alternately even and odd functions. So we evaluate only half the number of roots. 
		for i in range(1,int(polyorder)/2 +1):
			x=cos(pi*(i-0.25)/(polyorder+0.5))
			error=10*tolerance
		        iters=0
		        while (error>tolerance) and (iters<1000):
		                dx=-Legendre(polyorder,x)/DLegendre(polyorder,x)
		                x=x+dx
		                iters=iters+1
		                error=abs(dx)
			roots.append(x)
		# Use symmetry to get the other roots
		roots=array(roots)
		if polyorder%2==0:
			roots=concatenate( (-1.0*roots, roots[::-1]) )
		else:
			roots=concatenate( (-1.0*roots, [0.0], roots[::-1]) )
		err=0 # successfully determined roots
	return [roots, err]
##################################################################
# Weight coefficients
def GaussLegendreWeights(polyorder):
	W=[]
	[xis,err]=LegendreRoots(polyorder)
	if err==0:
		W=2.0/( (1.0-xis**2)*(DLegendre(polyorder,xis)**2) )
		err=0
	else:
		err=1 # could not determine roots - so no weights
	return [W, xis, err]



def getContourLevels(probDist,percLevels=[.997,.95,.6827]):
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
    cdf = cdf/float(max(cdf))
    
    # Determine point at which we reach 68% level
    percents = np.array(percLevels)
    levels = np.zeros(len(percents), dtype=float)
    for ii in range(len(levels)):
        # Get the index of the pixel at which the CDF
        # reaches this percentage (the first one found)
        idx = (np.where(cdf < percents[ii]))[0]
        
        # Now get the level of that pixel
        levels[ii] = pixSort[idx[-1]]
    return levels


def plotModelParam(flag='for_posterNOW',numbins=25,Cplot=True):
    posterior = np.loadtxt('/u/schappell/pmnOld/'+flag+'_.txt')#'_post_equal_weights.dat')
    priors = np.loadtxt('/u/schappell/pmnOld/'+flag+'_priors.txt')
    gamma = posterior[:,2]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha = posterior[:,3]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta = posterior[:,4]*(priors[2,1] - priors[2,0]) + priors[2,0]
    #tmp = math.sqrt(2.0) * priors[4,1] *special.erfinv(2.0*posterior[:,3] - 1.0) + (priors[4,0]/(math.sqrt(2.0)*priors[4,1]))
    #rbreak = math.sqrt(2.0) * priors[4,1] * tmp + priors[4,0]
    #rbreak = np.exp(posterior[:,5]*(math.log(priors[3,1]) - math.log(priors[3,0])) + math.log(priors[3,0]))
    rbreak = posterior[:,5]*(priors[3,1] - priors[3,0]) + priors[3,0]
    if (Cplot==True):
        cfact = posterior[:,6]
    weights = posterior[:,0]
        #print np.sum(weights)
    #else:
     #   weights = posterior[:,4]
    #weights -= np.min(weights)

    #1D posteriors
    py.clf()
    hist,bins,junk=py.hist(gamma,bins=4*numbins,normed=1,weights=weights)
    py.xlabel('$\gamma$ (Inner Slope)')
    py.ylabel('Weighted Posterior')
    print 'Gamma --------------------------'
    print 'Mean: '+str(np.average(gamma,weights=weights))+' +/- '+str(math.sqrt(np.average((gamma-np.average(gamma,weights=weights))**2,weights=weights)))
    print 'Median: '+str(np.median(gamma))
    tmpdex = np.argmax(hist)
    print 'Max bin: '+str(bins[tmpdex])+' to '+str(bins[tmpdex+1])
    levels = getContourLevels(hist)#,percLevels=[0.999999998027,0.999999426697,0.999936657516])
    tmpdex=np.where(hist > levels[2])[0]
    ax = gca()
    ylim = ax.get_ylim()
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    print '68% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[1])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    print '95% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[0])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    print '99.7% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    print ' '
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+flag+'_gamma_post.png')
    py.clf()

    hist,bins,junk=py.hist(alpha,bins=4*numbins,normed=1,weights=weights)
    py.xlabel(r'$\alpha$ (Outer Slope)')
    py.ylabel('Weighted Posterior')
    ax = gca()
    ylim = ax.get_ylim()
    print 'Alpha --------------------------'
    print 'Mean: '+str(np.average(alpha,weights=weights))+' +/- '+str(math.sqrt(np.average((alpha-np.average(alpha,weights=weights))**2,weights=weights)))
    print 'Median: '+str(np.median(alpha))
    tmpdex = np.argmax(hist)
    print 'Max bin: '+str(bins[tmpdex])+' to '+str(bins[tmpdex+1])
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    print '68% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[1])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    print '95% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[0])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    print '99.7% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    print ' '
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+flag+'_alpha_post.png')
    py.clf()

    hist,bins,junk=py.hist(delta,bins=4*numbins,normed=1,weights=weights)
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel('Weighted Posterior')
    ax = gca()
    ylim = ax.get_ylim()
    print 'Delta --------------------------'
    print 'Mean: '+str(np.average(delta,weights=weights))+' +/- '+str(math.sqrt(np.average((delta-np.average(delta,weights=weights))**2,weights=weights)))
    print 'Median: '+str(np.median(delta))
    tmpdex = np.argmax(hist)
    print 'Max bin: '+str(bins[tmpdex])+' to '+str(bins[tmpdex+1])
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    print '68% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[1])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    print '95% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[0])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    print '99.7% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    print ' '
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+flag+'_delta_post.png')
    py.clf()

    hist,bins,junk=py.hist(rbreak,bins=4*numbins,normed=1,weights=weights)
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel('Weighted Posterior')
    ax = gca()
    ylim = ax.get_ylim()
    print 'r_break --------------------------'
    print 'Mean: '+str(np.average(rbreak,weights=weights))+' +/- '+str(math.sqrt(np.average((rbreak-np.average(rbreak,weights=weights))**2,weights=weights)))
    print 'Median: '+str(np.median(rbreak))
    tmpdex = np.argmax(hist)
    print 'Max bin: '+str(bins[tmpdex])+' to '+str(bins[tmpdex+1])
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
    print '68% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[1])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='k')
    print '95% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    tmpdex=np.where(hist > levels[0])[0]
    py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
    print '99.7% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
    print ' '
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+flag+'_rbreak_post.png')
    py.clf()

    if (Cplot==True):
        hist,bins,junk=py.hist(cfact,bins=4*numbins,normed=1,weights=weights)
        py.xlabel('C Factor')
        py.ylabel('Weighted Posterior')
        ax = gca()
        ylim = ax.get_ylim()
        print 'C Factor --------------------------'
        print 'Mean: '+str(np.average(cfact,weights=weights))+' +/- '+str(math.sqrt(np.average((cfact-np.average(cfact,weights=weights))**2,weights=weights)))
        print 'Median: '+str(np.median(cfact))
        tmpdex = np.argmax(hist)
        print 'Max bin: '+str(bins[tmpdex])+' to '+str(bins[tmpdex+1])
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])[0]
        py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
        py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='limegreen')
        print '68% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
        tmpdex=np.where(hist > levels[1])[0]
        py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='k')
        py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='k')
        print '95% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
        tmpdex=np.where(hist > levels[0])[0]
        py.plot([np.min(bins[tmpdex]),np.min(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
        py.plot([np.max(bins[tmpdex]),np.max(bins[tmpdex])],ylim,linewidth=3.0,color='grey')
        print '99.7% interval : '+str(np.min(bins[tmpdex]))+' to '+str(np.max(bins[tmpdex+1]))
        print ' '
        py.ylim(ylim)
        py.savefig('/u/schappell/plots/'+flag+'_cfactor_post.png')
        py.clf()


    #2D contours
    py.clf()
    hist,ybins,xbins = np.histogram2d(gamma,alpha,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
    py.xlabel(r'$\alpha$ (Outer Slope)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/'+flag+'_gamma_alpha.png')
    print 'GAMMA AND ALPHA --------------------------'
    tmpdex = unravel_index(hist.argmax(),hist.shape)
    print 'GAMMA'
    print 'Max bin: '+str(np.median(Y[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[1])
    print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    print 'ALPHA'
    print 'Max bin: '+str(np.median(X[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    print ' '

    hist,ybins,xbins = np.histogram2d(gamma,delta,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/'+flag+'_gamma_delta.png')
    print 'GAMMA AND DELTA --------------------------'
    tmpdex = unravel_index(hist.argmax(),hist.shape)
    print 'GAMMA'
    print 'Max bin: '+str(np.median(Y[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    print 'DELTA'
    print 'Max bin: '+str(np.median(X[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    print ' '

    hist,ybins,xbins = np.histogram2d(gamma,rbreak,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/'+flag+'_gamma_rbreak.png')
    print 'GAMMA AND r_BREAK --------------------------'
    tmpdex = unravel_index(hist.argmax(),hist.shape)
    print 'GAMMA'
    print 'Max bin: '+str(np.median(Y[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    print 'r_BREAK'
    print 'Max bin: '+str(np.median(X[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    print ' '


    if (Cplot==True):
        hist,ybins,xbins = np.histogram2d(gamma,cfact,bins=numbins,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
        #        py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.ylabel('$\gamma$ (Inner Slope)')
        py.savefig('/u/schappell/plots/'+flag+'_gamma_cfactor.png')
        print 'GAMMA AND C FACTOR --------------------------'
        tmpdex = unravel_index(hist.argmax(),hist.shape)
        print 'GAMMA'
        print 'Max bin: '+str(np.median(Y[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        print 'C FACTOR'
        print 'Max bin: '+str(np.median(X[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        print ' '

        hist,ybins,xbins = np.histogram2d(alpha,cfact,bins=numbins,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.ylabel(r'$\alpha$ (Outer Slope)')
        py.savefig('/u/schappell/plots/'+flag+'_alpha_cfactor.png')
        print 'ALPHA AND C FACTOR --------------------------'
        tmpdex = unravel_index(hist.argmax(),hist.shape)
        print 'ALPHA'
        print 'Max bin: '+str(np.median(Y[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        print 'C FACTOR'
        print 'Max bin: '+str(np.median(X[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        print ' '

        hist,ybins,xbins = np.histogram2d(delta,cfact,bins=numbins,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.ylabel(r'$\delta$ (Sharpness)')
        py.savefig('/u/schappell/plots/'+flag+'_delta_cfactor.png')
        print 'DELTA AND C FACTOR --------------------------'
        tmpdex = unravel_index(hist.argmax(),hist.shape)
        print 'DELTA'
        print 'Max bin: '+str(np.median(Y[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        print 'C FACTOR'
        print 'Max bin: '+str(np.median(X[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        print ' '

        hist,ybins,xbins = np.histogram2d(rbreak,cfact,bins=numbins,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.ylabel(r'$r_{break}$ (prc)')
        py.savefig('/u/schappell/plots/'+flag+'_rbreak_cfactor.png')
        print 'r_BREAK AND C FACTOR --------------------------'
        tmpdex = unravel_index(hist.argmax(),hist.shape)
        print 'r_BREAK'
        print 'Max bin: '+str(np.median(Y[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
        print 'C FACTOR'
        print 'Max bin: '+str(np.median(X[tmpdex]))
        levels = getContourLevels(hist)
        tmpdex=np.where(hist > levels[2])
        print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[1]) 
        print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        tmpdex=np.where(hist > levels[0])
        print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
        print ' '

    hist,ybins,xbins = np.histogram2d(alpha,delta,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel(r'$\alpha$ (Outer Slope)')
    py.savefig('/u/schappell/plots/'+flag+'_alpha_delta.png')
    print 'ALPHA AND DELTA --------------------------'
    tmpdex = unravel_index(hist.argmax(),hist.shape)
    print 'ALPHA'
    print 'Max bin: '+str(np.median(Y[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    print 'DELTA'
    print 'Max bin: '+str(np.median(X[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    print ' '

    hist,ybins,xbins = np.histogram2d(alpha,rbreak,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel(r'$\alpha$ (Outer Slope)')
    py.savefig('/u/schappell/plots/'+flag+'_alpha_rbreak.png')
    print 'ALPHA AND r_BREAK --------------------------'
    tmpdex = unravel_index(hist.argmax(),hist.shape)
    print 'ALPHA'
    print 'Max bin: '+str(np.median(Y[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    print 'r_BREAK'
    print 'Max bin: '+str(np.median(X[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    print ' '

    hist,ybins,xbins = np.histogram2d(delta,rbreak,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel(r'$\delta$ (Sharpness)')
    py.savefig('/u/schappell/plots/'+flag+'_delta_rbreak.png')
    print 'DELTA AND r_BREAK --------------------------'
    tmpdex = unravel_index(hist.argmax(),hist.shape)
    print 'DELTA'
    print 'Max bin: '+str(np.median(Y[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(Y[tmpdex]))+' to '+str(np.max(Y[tmpdex]))
    print 'r_BREAK'
    print 'Max bin: '+str(np.median(X[tmpdex]))
    levels = getContourLevels(hist)
    tmpdex=np.where(hist > levels[2])
    print '68% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[1]) 
    print '95% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    tmpdex=np.where(hist > levels[0])
    print '99.7% interval : '+str(np.min(X[tmpdex]))+' to '+str(np.max(X[tmpdex]))
    print ' '


    #Make pyramid plot
    py.clf()
    #left,bottom,width,height,space = 0.12,0.12,0.2,0.2,0.02
    Xarray = np.zeros((10,numbins,numbins))
    Yarray = np.zeros((10,numbins,numbins))
    histarray = np.zeros((10,numbins,numbins))
    ij = 0
    for i in range(3):
        if (i==0):
            yarray = alpha * 1.0
        elif (i==1):
#            yarray = cfact * 1.0
#        elif (i==2):
            yarray = rbreak * 1.0
        else:
            yarray = delta * 1.0
        for j in range(3-i):
            if (j==0):
                xarray = gamma * 1.0
            elif (j==1):
                xarray = delta * 1.0
            elif (j==2):
                xarray = rbreak * 1.0
#            else:
#                xarray = cfact * 1.0
            hist,ybins,xbins = np.histogram2d(yarray,xarray,bins=numbins,weights=weights)
            x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
            y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
            X,Y=np.meshgrid(x,y)
            Xarray[ij,:,:] = X
            Yarray[ij,:,:] = Y
            histarray[ij,:,:] = hist
            ij += 1

    ij = 0
    xticks = np.array([1.0,2.0,0.5,0.1])
    yticks = np.array([2.0,0.5,2.0])
    py.figure(1,figsize=(8.0,8.0))
    for i in range(3):
        for j in range(3-i):
            levels = getContourLevels(histarray[ij,:,:])
            py.axes([0.085+(j*0.3),0.085+(i*0.3),0.27,0.27])
            py.contour(Xarray[ij,:,:],Yarray[ij,:,:],histarray[ij,:,:],levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid'))
            ij += 1
                #pdb.set_trace()
            ax = gca()
            ax.xaxis.set_major_locator(ticker.MultipleLocator(xticks[j]))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(yticks[i]))
            if (i > 0):
                ax.axes.xaxis.set_ticklabels([])
            if (j > 0):
                ax.axes.yaxis.set_ticklabels([])
            if (i == 0):
                if (j == 0):
                    py.xlabel(r'$\gamma$ (Inner Slope)')
                if (j == 1):
                    py.xlabel(r'$\delta$ (Sharpness)')
                if (j == 2):
                    py.xlabel(r'$r_{break}$ (pc)')
            #    if (j == 3):
             #       py.xlabel('C Factor')
            if (j == 0):
                if (i==0):
                    py.ylabel(r'$\alpha$ (Outer Slope)')
                elif (i==1):
            #        py.ylabel('C Factor')
             #   if (i==2):
                    py.ylabel(r'$r_{break}$ (pc)')
                else:
                    py.ylabel(r'$\delta$ (Sharpness)')
                #pdb.set_trace()
    #pdb.set_trace()
    py.savefig('/u/schappell/plots/'+flag+'_pyramidPOSTER.eps',format='eps',dpi=1000)


    #Make pyramid plot with 1D posteriors
    py.clf()
    #left,bottom,width,height,space = 0.12,0.12,0.2,0.2,0.02
    Xarray = np.zeros((10,numbins,numbins))
    Yarray = np.zeros((10,numbins,numbins))
    histarray = np.zeros((10,numbins,numbins))
    ij = 0
    for i in range(4):
        if (i==0):
            yarray = alpha * 1.0
        elif (i==1):
            yarray = cfact * 1.0
        elif (i==2):
            yarray = rbreak * 1.0
        else:
            yarray = delta * 1.0
        for j in range(4-i):
            if (j==0):
                xarray = gamma * 1.0
            elif (j==1):
                xarray = delta * 1.0
            elif (j==2):
                xarray = rbreak * 1.0
            else:
                xarray = cfact * 1.0
            hist,ybins,xbins = np.histogram2d(yarray,xarray,bins=numbins,weights=weights)
            x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
            y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])
            X,Y=np.meshgrid(x,y)
            Xarray[ij,:,:] = X
            Yarray[ij,:,:] = Y
            histarray[ij,:,:] = hist
            ij += 1

    ij = 0
    xticks = np.array([2.0,4.0,0.75,0.3,4.0])#,0.1])
    yticks = np.array([4.0,0.3,1.0,4.0])
    py.figure(1,figsize=(8.0,8.0))
    for i in range(5):
        for j in range(5-i):
            py.axes([0.07+(j*0.175),0.084+(i*0.18),0.165,0.165])
            if ((i+j) < 4):
                levels = getContourLevels(histarray[ij,:,:])
                py.contour(Xarray[ij,:,:],Yarray[ij,:,:],histarray[ij,:,:],levels=levels, 
                           colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid'))

                ij += 1
                #pdb.set_trace()
                ax = gca()
                ax.xaxis.set_major_locator(ticker.MultipleLocator(xticks[j]))
                ax.yaxis.set_major_locator(ticker.MultipleLocator(yticks[i]))
                if (i > 0):
                    ax.axes.xaxis.set_ticklabels([])
                if (j > 0):
                    ax.axes.yaxis.set_ticklabels([])
                if (i == 0):
                    if (j == 0):
                        py.xlabel(r'$\gamma$ (Inner Slope)')
                        gxlim = ax.get_xlim()
                    elif (j == 1):
                        py.xlabel(r'$\delta$ (Sharpness)')
                        dxlim = ax.get_xlim()
                    elif (j == 2):
                        py.xlabel(r'$r_{break}$ (pc)')
                        rxlim = ax.get_xlim()
                    elif (j == 3):
                        py.xlabel('C Factor')
                        cxlim = ax.get_xlim()
                        xtmp = ax.xaxis.get_major_ticks()
                        xtmp[0].label1.set_visible(False)
                if (j == 0):
                    if (i==0):
                        py.ylabel(r'$\alpha$ (Outer Slope)')
                    elif (i==1):
                        py.ylabel('C Factor')
                        ytmp = ax.yaxis.get_major_ticks()
                        ytmp[0].label1.set_visible(False)
                    elif (i==2):
                        py.ylabel(r'$r_{break}$ (pc)')
                    elif (i==3):
                        py.ylabel(r'$\delta$ (Sharpness)')

            else:
                ax = gca()
                yticks1 = np.array([0.1,10.0,1.0,0.06,1.0])
                #pdb.set_trace()
                ax.xaxis.set_major_locator(ticker.MultipleLocator(xticks[j]))
                ax.yaxis.set_major_locator(ticker.MultipleLocator(yticks1[i]))
                if (j > 0):
                    ax.axes.yaxis.tick_right()
                    ax.yaxis.set_label_position('right')
                if (j < 4):
                    ax.axes.xaxis.set_ticklabels([])
                if (j == 0):
                    hist,bins,junk=py.hist(gamma,weights=weights,normed=1,bins=100,histtype='step',color='b')
                    py.xlim(gxlim)
                    py.ylabel('Posterior')
                elif (j == 1):
                    hist,bins,junk=py.hist(delta,weights=weights,normed=1,bins=100,histtype='step',color='b')
                    py.xlim(dxlim)
                elif (j == 2):
                    hist,bins,junk=py.hist(rbreak,weights=weights,normed=1,bins=100,histtype='step',color='b')
                    py.xlim(rxlim)
                elif (j == 3):
                    hist,bins,junk=py.hist(cfact,weights=weights,normed=1,bins=100,histtype='step',color='b')
                    py.xlim(cxlim)
                elif (j == 4):
                    hist,bins,junk=py.hist(alpha,weights=weights,normed=1,bins=100,histtype='step',color='b')
                    py.xlabel(r'$\alpha$ (Outer Slope)')
                    ytmp = ax.yaxis.get_major_ticks()
                    ytmp[0].label1.set_visible(False)
                    xtmp = ax.yaxis.get_major_ticks()
                    xtmp[0].label1.set_visible(False)


                #pdb.set_trace()
    #pdb.set_trace()
    py.savefig('/u/schappell/plots/'+flag+'_pyramidALL.eps',format='eps',dpi=1000)





def plotGammaAlpha(file='/u/schappell/pmnOld/for_posterNOW_post_equal_weights.dat',numbins=15):
    posterior = np.loadtxt(file)
    gamma = posterior[:,0]
    alpha = posterior[:,1]
    #delta = posterior[:,2]
    #rbreak = posterior[:,3]
    weights = posterior[:,4]

    #weights *= 1e308
    weights -= np.min(weights)

    hist,ybins,xbins = np.histogram2d(gamma,alpha,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])

    levels = getContourLevels(hist)

    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'))
    #py.xlim([1,6])
    #py.ylim([-1.75,1.5])
    py.xlabel(r'$\alpha$ (Outer Slope)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/gamma_alpha_pmn.png')



def plotGammaDelta(file='/u/schappell/pmnOld/for_posterNOW_post_equal_weights.dat',numbins=15):
    posterior = np.loadtxt(file)
    gamma = posterior[:,0]
    delta = posterior[:,2]
    weights = posterior[:,4]

    #weights *= 1e308
    weights -= np.min(weights)

    hist,ybins,xbins = np.histogram2d(gamma,delta,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])

    levels = getContourLevels(hist)

    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'))
    #py.xlim([1,6])
    #py.ylim([-1.75,1.5])
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/gamma_delta_pmn.png')


def plotGamma_rbreak(file='/u/schappell/pmnOld/for_posterNOW_post_equal_weights.dat',numbins=15):
    posterior = np.loadtxt(file)
    gamma = posterior[:,0]
    #alpha = posterior[:,1]
    #delta = posterior[:,2]
    rbreak = posterior[:,3]
    weights = posterior[:,4]

    #weights *= 1e308
    weights -= np.min(weights)

    hist,ybins,xbins = np.histogram2d(gamma,rbreak,bins=numbins,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numbins)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numbins)])

    levels = getContourLevels(hist)

    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'))
    #py.xlim([1,6])
    #py.ylim([-1.75,1.5])
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/gamma_rbreak_pmn.png')



def multiStarPlot(starlist,fitlist=[-1],alnDir='14_06_18/',align='align/align_d_rms_1000_abs_t',tag='multipleStars',
                  poly='polyfit_nz/fit',polyj='polyfit_nzj/fit',points='points_nz/',oneColor=True):

    rootDir = '/g/ghez/align/'

    #plot multiple stars on the same
    #oneColor plots all the points for all the stars in black
    #fitlist, list of what order fit for each star, 0 for vel, 1 for accel, 2 for jerk
    #leave as [-1] for no fit overplot

    s = starset.StarSet(rootDir + alnDir + align)
    s.loadPolyfit(rootDir + alnDir + poly, accel=0, arcsec=0)
    s.loadPolyfit(rootDir + alnDir + poly, accel=1, arcsec=0)

    names = s.getArray('name')

    #Get fit info
    if (fitlist[0] != -1):

        #Get jerk info
        fitRoot = rootDir + alnDir + polyj
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

        jy0 = _fit[11].tonumpy()
        jvy = _fit[12].tonumpy()
        jay = _fit[13].tonumpy()
        jy = _fit[14].tonumpy()
        jy0e = _fit[15].tonumpy()
        jvye = _fit[16].tonumpy()
        jaye = _fit[17].tonumpy()
        jye = _fit[18].tonumpy()


    py.clf()
    for i in range(len(starlist)):
    
        ii = names.index(starlist[i])
        star = s.stars[ii]        

        pointsTab = asciidata.open(rootDir + alnDir + points + starlist[i] + '.points')

        time = pointsTab[0].tonumpy()
        x = pointsTab[1].tonumpy()
        y = pointsTab[2].tonumpy()
        xerr = pointsTab[3].tonumpy()
        yerr = pointsTab[4].tonumpy()

        #if ((starlist[i]== 'S0-24') | (starlist[i] == 'S0-59')):
        #    print 'skip S0-24'
        if (oneColor == True):
            py.errorbar(x,y,xerr=xerr,yerr=yerr,fmt='k.')
        else:
            py.errorbar(x,y,xerr=xerr,yerr=yerr)

        if (fitlist[0] == -1):
            continue
        else:
            if (fitlist[i] == 0):
                print starlist[i] + ' fit with velocity'
                fitxv = star.fitXv
                fityv = star.fitYv
                dtv = time - fitxv.t0
                fitLineXv = fitxv.p + (fitxv.v * dtv)
                fitLineYv = fityv.p + (fityv.v * dtv)
                py.plot(fitLineXv,fitLineYv,'b-')

            elif (fitlist[i] == 1):
                print starlist[i] + ' fit with acceleration'
                fitxa = star.fitXa
                fitya = star.fitYa
                dta = time - fitxa.t0
                fitLineXa = fitxa.p + (fitxa.v * dta) + (fitxa.a * dta * dta / 2.0)
                fitLineYa = fitya.p + (fitya.v * dta) + (fitya.a * dta * dta / 2.0)
                py.plot(fitLineXa,fitLineYa,'b-')

            elif (fitlist[i] == 2):
                ij = np.where(jnames == starlist[i])[0]
                print starlist[i] + ' fit with jerk'
                dtj = time - jt0x[ij]
                fitLineXj = jx0[ij] + (jvx[ij] * dtj) + (jax[ij] * dtj * dtj / 2.0) + (jx[ij] * dtj * dtj * dtj / 6.0)
                fitLineYj = jy0[ij] + (jvy[ij] * dtj) + (jay[ij] * dtj * dtj / 2.0) + (jy[ij] * dtj * dtj * dtj / 6.0)
                py.plot(fitLineXj,fitLineYj,'b-')

    #py.xlim([-.5,.5])
    #py.ylim([-.5,.5])
    py.xlim([.53,.68])
    py.ylim([0.05,.2])
    py.xlabel('X (arcsec)')
    py.ylabel('Y (arcsec)')

    py.savefig('/u/schappell/plots/'+tag+'.png')



def gammaHist_do_over(getOnesig=False):
    tdo = np.loadtxt('/u/schappell/Downloads/gamma.txt')
    sch = np.loadtxt('/u/schappell/pmnOld/for_posterNOW_post_equal_weights.dat')

    tdo_gamma = tdo[:,0]
    tmp_hist = tdo[:,1]

    sch_hist,bins,patches=py.hist(sch[:,0],bins=35,normed=True)
    tdo_hist = np.zeros(len(sch_hist))

    for i in range(len(sch_hist)):
        tmp = np.where((tdo_gamma >= bins[i]) & (tdo_gamma < bins[i+1]))
        tdo_hist[i] = np.sum(tmp_hist[tmp])

    tdo_hist *= np.sum(sch_hist)/np.sum(tdo_hist)

    bin_width = bins[1]-bins[0]
    py.clf()
    py.bar(bins[:-1],tdo_hist,width=bin_width)
    if (getOnesig==True):
        one_sig= [14,15,16,17,18,19,20]
        py.bar(bins[one_sig],tdo_hist[one_sig],width=bin_width,color='gold')
    py.bar(bins[:-1],sch_hist,width=bin_width,color='grey')
    if (getOnesig==True):
        one_sig=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
        py.bar(bins[one_sig],sch_hist[one_sig],width=bin_width,color='limegreen')
    py.ylabel('Normalized Posterior')
    py.xlabel('$\gamma$')
    py.savefig('/u/schappell/gamm_hist_tdo_over.png')
    print bins



def plotMockData(gamma=-1.0,alpha=4.0,delta=4.0,r_break=0.5,R2d_cut=1.7):

    tmp = np.loadtxt('/u/schappell/code/c/stars_mn.dat')
    pOld = tmp[:,3]
    tmp = np.loadtxt('/u/schappell/code/c/stars_r.dat')
    radius = tmp[:,0]
    radius /= cm_in_pc
    R2d_cut *= dist / au_in_pc
    #bins = np.array([0.01,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    bins = np.array([0.01,0.025,0.05,0.075,0.1,0.2,0.4,0.6,0.8,1.0])
    numBins = len(bins) - 1

    py.clf()
    hist,bins,tmp=py.hist(radius,bins=bins,weights=pOld)
    py.xlabel('Radius (pc)')
    py.ylabel('Number')
    py.title('Gamma: '+str(gamma)+' Alpha: '+str(alpha)+' Delta: '+str(delta)+
             ' R_break:'+str(r_break)+' pc')
    py.savefig('/u/schappell/plots/starHist_mock_g'+str(gamma)+'_a'+str(alpha)+'_d'+str(delta)+
               '_br'+str(r_break)+'_pc.png')
    bsize = np.array([bins[i+1] - bins[i] for i in range(numBins)])
   # volume = np.array([bins[i+1]**3 - bins[i]**3 for i in range(numBins)])*4.0*pi/3.0
    intVol = np.zeros(numBins+1)
    for i in range(numBins+1):
        if (bins[i] < R2d_cut):
            intVol[i] = 4.0*pi*bins[i]**3/3.0
            print 'Using sphere for below R2d cut'
        else:
            intVol[i] = pi*R2d_cut**2*bins[i]

    volume = np.array([intVol[i+1] - intVol[i] for i in range(numBins)])

    density = hist/volume
    r3D = np.array([(bins[i]+bins[i+1])/2.0 for i in range(numBins)])
    rLeft = bins[:numBins]

    error_density = (np.sqrt(hist)/volume)
    for i in range(numBins):
        if (error_density[i] >= density[i]):
            error_density[i] = density[i] - 100.0

    #pdb.set_trace()
    py.clf()
    py.bar(rLeft,density,width=bsize,bottom=100.0,color='w')
    py.errorbar(r3D,density,yerr=error_density,fmt='o')
    py.xscale('log')
    py.yscale('log')
    py.xlabel('Log Radius (pc)')
    #py.xlim([0.04,1.0])
    py.ylabel('Log Density (pc^-3)')
    py.title('Gamma: '+str(gamma)+' Alpha: '+str(alpha)+' Delta: '+str(delta)+
             ' R_break:'+str(r_break)+' pc')
    py.savefig('/u/schappell/plots/roughRho_mock_g'+str(gamma)+'_a'+str(alpha)+'_d'+str(delta)+
               '_br'+str(r_break)+'_pc.png')



def projectionDensity(numBins=70,nonRadial=1,gamma=-1.0,alpha=4.0,delta=9.0,
                      rbreak=0.5,numR=1000,maxZ=1.0,cfact=0.31,gcowsNum=4,label='',label_mock='_mock'):

    cfact /= (1.0 - cfact)
    gcows_stars = np.loadtxt('/u/schappell/code/c/stars_mn'+label+'.dat')
    maser_stars = np.loadtxt('/u/schappell/code/c/maser_mn'+label+'.dat')

    gcowsHist, gcows_bins = np.histogram(gcows_stars[:,0]/cm_in_au/dist,bins=gcowsNum)

    if (nonRadial>0):
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,1]
        gcows_cols = gcows[:,0]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        gcows_surface = np.zeros(len(gcowsHist))
        for i in range(len(gcowsHist)):
            gcowsdex = np.where((gcows_R > gcows_bins[i]) & (gcows_R <= gcows_bins[i+1]))[0]
            gcows_surface[i] = plate_scale**2 * len(gcowsdex)
    else:
        gcows_surface = np.array([gcows_bins[i+1]**2 - gcows_bins[i]**2 for i in range(len(gcows_bins)-1)])*pi
    gcows_middle = np.array([(gcows_bins[i+1]+gcows_bins[i])/2.0 for i in range(len(gcows_bins)-1)])

    maserHist,maser_bins = np.histogram(maser_stars[:,0]/cm_in_au/dist,bins=gcows_bins)
    maser_surface = np.array([maser_bins[i+1]**2 - maser_bins[i]**2 for i in range(len(maser_bins)-1)])*pi
    maser_middle = np.array([(maser_bins[i+1]+maser_bins[i])/2.0 for i in range(len(maser_bins)-1)])

    pdb.set_trace()
    py.clf()
    norm = cfact * np.sum(maserHist/maser_surface)*(maser_bins[1]-maser_bins[0])
    norm += np.sum(gcowsHist/gcows_surface)*(gcows_bins[1]-gcows_bins[0])

    modR_one,modD_one=projectBPL(gamma=gamma,alpha=alpha,delta=delta,rbreak=rbreak,min_R=np.min(gcows_bins)*dist/au_in_pc,
                                 max_R=np.max(maser_bins)*dist/au_in_pc,numR=numR,maxZ=maxZ,norm=norm)

    #if (envelope==True):
     #   modArray = np.zeros((16,len(modD_one)))
      #  tmpDex = 0
        #gval = np.array([(gminmax[1]-gminmax[0])*i/9.0 for i in range(10)])+gminmax[0]
        #aval = np.array([(aminmax[1]-aminmax[0])*i/9.0 for i in range(10)])+aminmax[0]
        #dval = np.array([(dminmax[1]-dminmax[0])*i/9.0 for i in range(10)])+dminmax[0]
        #rval = np.array([(rminmax[1]-rminmax[0])*i/9.0 for i in range(10)])+rminmax[0]
       # for i in range(len(gminmax)):
        #    for j in range(len(aminmax)):
         #       for k in range(len(dminmax)):
          #          for l in range(len(rminmax)):
           #             modR,modD=projectBPL(gamma=gminmax[i],alpha=aminmax[j],delta=dminmax[k],rbreak=rminmax[l],
            #                                 min_R=np.min(gcows_bins)*dist/au_in_pc,max_R=np.max(maser_bins)*dist/au_in_pc,
             #                                numR=numR,maxZ=maxZ,norm=norm)
              #          modArray[tmpDex] = modD
               #         tmpDex += 1
                        
                   # py.plot(modR*au_in_pc/dist,modD*dist/au_in_pc,color='gray')
        #for i in range(16):
         #   for j in range(16):
          #      if (i > j):
           #         py.fill_between(modR*au_in_pc/dist,modArray[i]*dist/au_in_pc,modArray[j]*dist/au_in_pc,color='gray')


    py.errorbar(gcows_middle,gcowsHist/gcows_surface,yerr=np.sqrt(gcowsHist)/gcows_surface,fmt='ok')
    py.errorbar(maser_middle,cfact*maserHist/maser_surface,yerr=np.sqrt(maserHist)/maser_surface,fmt='ob')

    py.plot(modR_one*au_in_pc/dist,modD_one*dist/au_in_pc,'g')

    py.xscale('log')
    py.yscale('log')
    py.xlim([0.2,20])
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/surfaceDensity_real.png')

    mockData_maser = np.loadtxt('/u/schappell/code/c/maser_mn'+label_mock+'.dat')
    mockData_gcows = np.loadtxt('/u/schappell/code/c/stars_mn'+label_mock+'.dat')

    gcowsHist, gcows_bins = np.histogram(mockData_gcows[:,0]/cm_in_au/dist,bins=gcowsNum)

    if (nonRadial>0):
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,1]
        gcows_cols = gcows[:,0]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        gcows_surface = np.zeros(len(gcowsHist))
        py.close()
        for i in range(len(gcowsHist)):
            gcowsdex = np.where((gcows_R > gcows_bins[i]) & (gcows_R <= gcows_bins[i+1]))[0]
            gcows_surface[i] = plate_scale**2 * len(gcowsdex)
            py.plot(gcows_cols[gcowsdex],gcows_rows[gcowsdex],'.')

        mockxy = np.loadtxt('/u/schappell/code/c/xy_pixels'+label+'.dat')
        py.plot(mockxy[:,0],mockxy[:,1],'.')
        py.xlabel('Pixels')
        py.ylabel('Pixels')
        py.savefig('/u/schappell/plots/gcows_area.png')
        #pdb.set_trace()
    else:
        gcows_surface = np.array([gcows_bins[i+1]**2 - gcows_bins[i]**2 for i in range(len(gcows_bins)-1)])*pi
    gcows_middle = np.array([(gcows_bins[i+1]+gcows_bins[i])/2.0 for i in range(len(gcows_bins)-1)])

    maserHist,maser_bins = np.histogram(mockData_maser[:,0]/cm_in_au/dist,bins=numBins)
    maser_surface = np.array([maser_bins[i+1]**2 - maser_bins[i]**2 for i in range(len(maser_bins)-1)])*pi
    maser_middle = np.array([(maser_bins[i+1]+maser_bins[i])/2.0 for i in range(len(maser_bins)-1)])

    py.clf()
    py.errorbar(gcows_middle,gcowsHist/gcows_surface,yerr=np.sqrt(gcowsHist)/gcows_surface,fmt='ok')
    py.errorbar(maser_middle,cfact*maserHist/maser_surface,yerr=np.sqrt(maserHist)/maser_surface,fmt='ob')

    norm = cfact * np.sum(maserHist/maser_surface)*(maser_bins[1]-maser_bins[0])
    norm += np.sum(gcowsHist/gcows_surface)*(gcows_bins[1]-gcows_bins[0])
    modR,modD=projectBPL(gamma=gamma,alpha=alpha,delta=delta,rbreak=rbreak,min_R=np.min(gcows_bins)*dist/au_in_pc,
                         max_R=np.max(maser_bins)*dist/au_in_pc,numR=numR,maxZ=maxZ,norm=norm)
    
    py.plot(modR*au_in_pc/dist,modD*dist/au_in_pc,'g')

    py.xscale('log')
    py.yscale('log')
    #py.ylim([2,20])
    #py.xlim([0.01,20])
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/surfaceDensity_mock.png')
    py.clf()




def projectRhoEnvelope(numBins=70,nonRadial=1,numR=50,maxZ=1.0,gcowsNum=4,datLabel='',datDir='/u/schappell/code/c/',
                       mnlabel='',mndir='/u/schappell/pmnOld/',R2d_cut=5.0,cfact=0.27):

    posterior = mndir+mnlabel+'post_equal_weights.dat'
    priors = mndir+mnlabel+'priors.txt'

    gcows_stars = np.loadtxt(datDir+'stars_mn'+datLabel+'.dat')

    maser_stars = np.loadtxt(datDir+'maser_mn'+datLabel+'.dat')
    maser_R2d = np.sort(maser_stars[:,0]/cm_in_au/dist)
    binSize = len(maser_R2d) / numBins
    maser_bins = np.zeros(numBins+1)
    for i in range(numBins):
        maser_bins[i] = maser_R2d[i*binSize]
    maser_bins[numBins] = np.max(maser_R2d)
    maserHist,maser_bins = np.histogram(maser_stars[:,0]/cm_in_au/dist,bins=maser_bins)
    
    #maserHist *= cfact / ( 1.0 - cfact)

    maser_surface = np.array([maser_bins[i+1]**2 - maser_bins[i]**2 for i in range(len(maser_bins)-1)])*pi
    maser_middle = np.array([(maser_bins[i+1]+maser_bins[i])/2.0 for i in range(len(maser_bins)-1)])

    gcows_R2d = np.sort(gcows_stars[:,0]/cm_in_au/dist)
    binSize = len(gcows_R2d) / gcowsNum
    gcows_bins = np.zeros(gcowsNum+1)
    for i in range(gcowsNum):
        gcows_bins[i] = gcows_R2d[i*binSize]
    gcows_bins[gcowsNum] = np.max(gcows_R2d)
    gcowsHist, gcows_bins = np.histogram(gcows_stars[:,0]/cm_in_au/dist,bins=gcows_bins,weights=gcows_stars[:,3])

    if (nonRadial>0):
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,1]
        gcows_cols = gcows[:,0]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        gcows_surface = np.zeros(len(gcowsHist))
        for i in range(len(gcowsHist)):
            gcowsdex = np.where((gcows_R > gcows_bins[i]) & (gcows_R <= gcows_bins[i+1]))[0]
            gcows_surface[i] = plate_scale**2 * len(gcowsdex)
    else:
        gcows_surface = np.array([gcows_bins[i+1]**2 - gcows_bins[i]**2 for i in range(len(gcows_bins)-1)])*pi
    gcows_middle = np.array([(gcows_bins[i+1]+gcows_bins[i])/2.0 for i in range(len(gcows_bins)-1)])

    maser_diff = np.array([(maser_bins[i+1] - maser_bins[i]) for i in range(len(maserHist))])
    gcows_diff = np.array([(gcows_bins[i+1] - gcows_bins[i]) for i in range(len(gcowsHist))])
    norm = np.sum(maserHist*maser_diff/maser_surface)
    norm += np.sum(gcowsHist*gcows_diff/gcows_surface)

    modR,modD,modLow,modHigh = projBPL_median_sigma(posterior,priors,maxZ=maxZ,min_R=np.min(gcows_bins)*dist/au_in_pc,
                                                    max_R=np.max(maser_bins)*dist/au_in_pc,numR=numR,R2d_cut=R2d_cut,norm=norm)

    py.clf()

    gcDex = np.where(modR*au_in_pc/dist < 0.98*R2d_cut)[0]
    py.fill_between(modR[gcDex]*au_in_pc/dist,modLow[gcDex]*dist/au_in_pc,modHigh[gcDex]*dist/au_in_pc,color='gray')

    scDex = np.where(modR*au_in_pc/dist >= R2d_cut)[0]
    py.fill_between(modR[scDex]*au_in_pc/dist,modLow[scDex]*dist/au_in_pc*cfact/(1.0-cfact),
                    modHigh[scDex]*dist/au_in_pc*cfact/(1.0-cfact),color='gray')

    maserHist *= cfact / (1.0 - cfact)

    py.errorbar(gcows_middle,gcowsHist/gcows_surface,yerr=np.sqrt(gcowsHist)/gcows_surface,fmt='ok')
    py.errorbar(maser_middle,maserHist/maser_surface,yerr=np.sqrt(maserHist)/maser_surface,fmt='ob')

    py.plot(modR[gcDex]*au_in_pc/dist,modD[gcDex]*dist/au_in_pc,'g')
    py.plot(modR[scDex]*au_in_pc/dist,modD[scDex]*dist/au_in_pc*cfact/(1.0-cfact),'g')

    py.xscale('log')
    py.yscale('log')
    py.xlim([0.2,20])
    py.ylim([0.5,20])
    py.xlabel('Radius (arcsec)')
    py.ylabel("Surface Density K' > 15.5 mag (arcsec$^{-2}$)")
    py.savefig('/u/schappell/plots/surfaceDensity_oneSigma.png')

    #delta_mass = modD * modR * (modR[1]-modR[0])
    #delta_mass_low = modLow * modR * (modR[1] - modR[0])
    #delta_mass_high = modHigh * modR * (modR[1] - modR[0])
    #pdb.set_trace()


def projectRhoEnvelopeSingle(numBins=4,nonRadial=1,numR=500,maxZ=1.0,datLabel='',datDir='/u/schappell/code/c/',
                             mnlabel='',mndir='/u/schappell/pmnOld/',R2d_cut=5.0):

    posterior = mndir+mnlabel+'post_equal_weights.dat'
    priors = mndir+mnlabel+'priors.txt'

    gcows_stars = np.loadtxt(datDir+'stars_mn'+datLabel+'.dat')

    gcowsHist, gcows_bins = np.histogram(gcows_stars[:,0]/cm_in_au/dist,bins=numBins)

    if (nonRadial>0):
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,1]
        gcows_cols = gcows[:,0]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        gcows_surface = np.zeros(len(gcowsHist))
        for i in range(len(gcowsHist)):
            gcowsdex = np.where((gcows_R > gcows_bins[i]) & (gcows_R <= gcows_bins[i+1]))[0]
            gcows_surface[i] = plate_scale**2 * len(gcowsdex)
    else:
        gcows_surface = np.array([gcows_bins[i+1]**2 - gcows_bins[i]**2 for i in range(len(gcows_bins)-1)])*pi
    gcows_middle = np.array([(gcows_bins[i+1]+gcows_bins[i])/2.0 for i in range(len(gcows_bins)-1)])

    norm = np.sum(gcowsHist/gcows_surface)*(gcows_bins[1]-gcows_bins[0])

    R2d_cut *= dist / au_in_pc

    posterior = np.loadtxt(posterior)
    priors = np.loadtxt(priors)
    gamma_array = posterior[:,0]*(priors[1] - priors[0]) + priors[0]

    dR = R2d_cut/numR
    R2d = np.array([i*dR+dR for i in range(numR)])
    median_rho = np.zeros(len(R2d))
    lower_rho = np.zeros(len(R2d))
    upper_rho = np.zeros(len(R2d))
    rho = np.zeros((len(gamma_array),len(R2d)))

    for i in range(len(R2d)):
        #tmp_rho = np.zeros(len(gamma_array))
        #print 'Start R2d: '+str(R2d[i])+' pc, '+str(i)
        for j in range(len(gamma_array)):
            rho[j,i] = maxZ * (R2d[i]**2 + maxZ**2)**(gamma_array[j]/-2.0)*((maxZ**2/R2d[i]**2)+1.0)**(gamma_array[j]/2.0)
            rho[j,i] *= scipy.special.hyp2f1(0.5,gamma_array[j]/2.0,1.5,-1.0*maxZ**2/R2d[i]**2)

    for j in range(len(gamma_array)):
        rho[j,:] /= np.sum(rho[j,:])
    for i in range(len(R2d)):
        tmp_median = np.median(rho[:,i])
        median_rho[i] = np.median(rho[:,i])

        tmp_percent = 0.0
        delta_rho = tmp_median * 0.001
        upper_sigma = tmp_median * 1.0
        #print 'Found median, now find lower and upper sigma'
        while ((tmp_percent < 0.341) & (upper_sigma < np.max(rho[:,i]))):
            upper_sigma += delta_rho
            tmpdex = np.where((rho[:,i] >= tmp_median) & (rho[:,i] <= upper_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        tmp_percent = 0.0
        lower_sigma = tmp_median * 1.0
        while ((tmp_percent < 0.341) & (lower_sigma > np.min(rho[:,i]))):
            lower_sigma -= delta_rho
            tmpdex = np.where((rho[:,i] <= tmp_median) & (rho[:,i] >= lower_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        #pdb.set_trace()
        #lower_sigma = min(tmp_rho)
        #upper_sigma = max(tmp_rho)

        lower_rho[i] = lower_sigma * 1.0
        upper_rho[i] = upper_sigma * 1.0

    norm = norm / np.sum(median_rho*dR)
    median_rho *= norm
    lower_rho *= norm
    upper_rho *= norm

    py.clf()

    py.fill_between(R2d*au_in_pc/dist,lower_rho*dist/au_in_pc,upper_rho*dist/au_in_pc,color='gray')

    py.errorbar(gcows_middle,gcowsHist/gcows_surface,yerr=np.sqrt(gcowsHist)/gcows_surface,fmt='ok')

    py.plot(R2d*au_in_pc/dist,median_rho*dist/au_in_pc,'g')

    py.xscale('log')
    py.yscale('log')
    py.xlim([0.2,5])
    py.ylim([0.4,10])
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/surfaceDensity_SINGLE_oneSigma.png')



def only_projDen(numBins=70,alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,nonRadial=1,
                 align='align/align_d_rms_1000_abs_t',poly='polyfit_nz/fit',points='points_nz/',
                 starlist='all',chainsDir='efit/chains/',innerCut=5.0,Rcut=5.0,nEpochs=14.,magCut=15.5,
                 globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,outerCut=10.0,ak_correct=True,gcowsNum=8,
                 maserDate='13jullgs1',numR=1000,maxZ=1.0,onlySchodel=False,magCutMaser=17.75,schodel=True):

    if (onlySchodel==False):
        names,r2d,ar,are,mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                         Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                         chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)

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


    #Index of stars used for multinest
        pdex = np.where(oldProb > 0.0)[0]

        r2d_accel = r2d[pdex] / cm_in_au / dist #r2d from accelerations in arcsec
        oldProb_accel = oldProb[pdex]
        mag_accel = mag[pdex]

    R2d,oldProb,mag = maserStars_out(maserDate=maserDate,innerCut=innerCut,outerCut=outerCut,magCut=magCutMaser,schodel=schodel)

    if (onlySchodel==False):

        gcowsHist, gcows_bins = np.histogram(r2d_accel,bins=gcowsNum)

        if (nonRadial>0):
            gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
            gcows_rows = gcows[:,1]
            gcows_cols = gcows[:,0]
            plate_scale = 0.00995
            gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
            gcows_surface = np.zeros(len(gcowsHist))
            for i in range(len(gcowsHist)):
                gcowsdex = np.where((gcows_R > gcows_bins[i]) & (gcows_R <= gcows_bins[i+1]))[0]
                gcows_surface[i] = plate_scale**2 * len(gcowsdex)
        else:
            gcows_surface = np.array([gcows_bins[i+1]**2 - gcows_bins[i]**2 for i in range(len(gcows_bins)-1)])*pi
        gcows_middle = np.array([(gcows_bins[i+1]+gcows_bins[i])/2.0 for i in range(len(gcows_bins)-1)])

    maserHist,maser_bins = np.histogram(R2d,bins=numBins)
    maser_surface = np.array([maser_bins[i+1]**2 - maser_bins[i]**2 for i in range(len(maser_bins)-1)])*pi
    maser_middle = np.array([(maser_bins[i+1]+maser_bins[i])/2.0 for i in range(len(maser_bins)-1)])

    py.clf()
    if (onlySchodel==False):
        py.errorbar(gcows_middle,gcowsHist/gcows_surface,yerr=np.sqrt(gcowsHist)/gcows_surface,fmt='ok',label='GCOWS')
    py.errorbar(maser_middle,maserHist/maser_surface,yerr=np.sqrt(maserHist)/maser_surface,fmt='ob',label='Schodel')

    py.xscale('log')
    py.yscale('log')
    if (onlySchodel == False):
        py.xlim([0.2,20])
        py.legend()
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/surfaceDensity_only_real.png')

    py.clf()
    bins = np.array([6.0+0.5*i for i in range(29)])
    if (onlySchodel==False):
        tmpDex = np.where(mag < magCut)[0]
        weights = mag * 0.0 + float(len(mag_accel))/float(len(tmpDex))
        py.hist(mag,bins=bins,weights=weights,label='Schodel',histtype='step',linewidth=2.0)
        py.hist(mag_accel,bins=bins,label='GCOWS',histtype='step',linewidth=2.0)
        py.legend(loc=2)
        py.ylabel('Number stars (Normalized to GCOWS < 15.5 mag)')

    py.xlabel("K' (mag)")
    py.savefig('/u/schappell/plots/KLF_hist.png')
        




def combinedMNhist(gamma,flag):
    histInfo = np.loadtxt('/u/schappell/pmnOld/'+flag+'_hist_bins.dat')
    middle = histInfo[:,0]
    number = histInfo[:,1]

    levels = getContourLevels(number)

    py.clf()
    py.plot(middle,number,'b')
    py.plot([gamma,gamma],[0,100],'k--')
    py.plot([-3,2],[levels[0],levels[0]],color='limegreen')
    py.plot([-3,2],[levels[1],levels[1]],color='gray')
    py.plot([-3,2],[levels[2],levels[2]],color='k')

    py.savefig('/u/schappell/plots/gammaHist_'+flag+'.png')



def compareMocksets(g1=-1.0,g2=-0.5,g3=-0.1,g4=1.0,a=6.0,d=2.0,br=0.5,numS=27,Rcut=1.7):
    stats1 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g1)+'_a'+str(a)+'_d'+
                        str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                        str(Rcut)+'arcsec_meanG_maxLG_evid.dat')
    stats2 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g2)+'_a'+str(a)+'_d'+
                        str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                        str(Rcut)+'arcsec_meanG_maxLG_evid.dat')
    stats3 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g3)+'_a'+str(a)+'_d'+
                        str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                        str(Rcut)+'arcsec_meanG_maxLG_evid.dat')
    stats4 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g4)+'_a'+str(a)+'_d'+
                        str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                        str(Rcut)+'arcsec_meanG_maxLG_evid.dat')

    info1 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g1)+'_a'+str(a)+'_d'+
                       str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                       str(Rcut)+'arcsec_hist_gamma_r2d_3d_logar_bins.dat')
    info2 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g2)+'_a'+str(a)+'_d'+
                       str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                       str(Rcut)+'arcsec_hist_gamma_r2d_3d_logar_bins.dat')
    info3 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g3)+'_a'+str(a)+'_d'+
                       str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                       str(Rcut)+'arcsec_hist_gamma_r2d_3d_logar_bins.dat')
    info4 = np.loadtxt('/u/schappell/pmnOld/broken_g'+str(g4)+'_a'+str(a)+'_d'+
                       str(d)+'_br'+str(br)+'_nStar'+str(numS)+'_Rcut'+
                       str(Rcut)+'arcsec_hist_gamma_r2d_3d_logar_bins.dat')

    label1 = 'Gamma: '+str(g1)
    label2 = 'Gamma: '+str(g2)
    label3 = 'Gamma: '+str(g3)
    label4 = 'Gamma: '+str(g4)
    py.clf()
    a,b,c=py.hist(stats1[:,2],bins=20,label=label1,histtype='step',normed=1,linewidth=1.5)
    py.hist(stats2[:,2],bins=b,label=label2,histtype='step',normed=1,linewidth=1.5)
    py.hist(stats3[:,2],bins=b,label=label3,histtype='step',normed=1,linewidth=1.5)
    py.hist(stats4[:,2],bins=b,label=label4,histtype='step',normed=1,linewidth=1.5)
    py.legend()
    py.xlabel('Evidence')
    py.ylabel('Normalized Frequency')
    py.savefig('/u/schappell/plots/compareMockSets_evidence.png')

    py.clf()
    py.plot(info1[:,0],info1[:,1]/np.sum(info1[:,1]),linewidth=1.5,label=label1)
    py.plot(info2[:,0],info2[:,1]/np.sum(info2[:,1]),linewidth=1.5,label=label2)
    py.plot(info3[:,0],info3[:,1]/np.sum(info3[:,1]),linewidth=1.5,label=label3)
    py.plot(info4[:,0],info4[:,1]/np.sum(info4[:,1]),linewidth=1.5,label=label4)
    py.legend(loc=2)
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Combined Normalized Posterior')
    py.savefig('/u/schappell/plots/compareMockSets_gamma.png')

    py.clf()
    a,b,c=py.hist(info1[:,2]/(dist*cm_in_au),weights=info1[:,3]/100,label=label1,histtype='step',linewidth=1.5)
    py.hist((info2[:,2]/(dist*cm_in_au)),weights=(info2[:,3]/100),bins=b,label=label2,histtype='step',linewidth=1.5)
    py.hist((info3[:,2]/(dist*cm_in_au)),weights=(info3[:,3]/100),bins=b,label=label3,histtype='step',linewidth=1.5)
    py.hist((info4[:,2]/(dist*cm_in_au)),weights=(info4[:,3]/100),bins=b,label=label4,histtype='step',linewidth=1.5)
   # py.errorbar(info1[:,2]/(dist*cm_in_au),info1[:,3]/100,yerr=np.sqrt(info1[:,3]/100),label=label1,fmt='o')
   # py.errorbar(info2[:,2]/(dist*cm_in_au),info2[:,3]/100,yerr=np.sqrt(info2[:,3]/100),label=label2,fmt='o')
    py.legend(loc=2)
    py.xlabel('R2d (arcsec)')
    py.ylabel('Number Stars')
    py.savefig('/u/schappell/plots/compareMockSets_R2d.png')

    py.clf()
    #b = np.array([0.0,0.01,0.025,0.05,0.075,0.1,0.2,0.4,0.6,0.8,1.0])
    b = 10.0**(np.arange(-2,0.0,0.2))
    py.hist(info1[:,4]/cm_in_pc,weights=info1[:,5]/100,bins=b,label=label1,histtype='step',linewidth=1.5)
    py.hist(info2[:,4]/cm_in_pc,weights=(info2[:,5]/100),bins=b,label=label2,histtype='step',linewidth=1.5)
    py.hist(info3[:,4]/cm_in_pc,weights=(info3[:,5]/100),bins=b,label=label3,histtype='step',linewidth=1.5)
    py.hist(info4[:,4]/cm_in_pc,weights=(info4[:,5]/100),bins=b,label=label4,histtype='step',linewidth=1.5)
    #py.plot([0.1,0.1],[0,12],'k--')
    #py.errorbar(info1[:,4]/cm_in_pc,info1[:,5]/100,yerr=np.sqrt(info1[:,5]/100),label=label1,fmt='o')
    #py.errorbar(info2[:,4]/cm_in_pc,info2[:,5]/100,yerr=np.sqrt(info2[:,5]/100),label=label2,fmt='o')
    py.legend(loc=2)
    py.xscale('log')
    py.xlim([0.02,1.0])
    py.ylim([0,10])
    py.xlabel('r 3D (pc)')
    py.ylabel('Number Stars')
    py.savefig('/u/schappell/plots/compareMockSets_r3d.png')

    py.clf()
  #  py.errorbar(info1[:,6],info1[:,7]/numS1,yerr=np.sqrt(info1[:,7]/numS1),label=label1,fmt='o')
  #  py.errorbar(info2[:,6],info2[:,7]/numS2,yerr=np.sqrt(info2[:,7]/numS2),label=label2,fmt='o')
  #  py.legend()
  #  py.xlabel('a_r')
  #  py.ylabel('Number')
  #  py.savefig('/u/schappell/plots/compareMockSets_ar.png')
  #  py.clf()




def compareTwoPost(file1='data_mthread_FINAL_free_maser_cuts_5.0_15.0arcsec',label1='W/ Accel',
                   file2='data_mthread_NOACCEL_free_maser_cuts_5.0_15.0arcsec',label2='W/o Accel',numBins=50):

    posterior1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_.txt')
    priors1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_priors.txt')
    gamma1 = posterior1[:,2]*(priors1[0,1] - priors1[0,0]) + priors1[0,0]
    alpha1 = posterior1[:,3]*(priors1[1,1] - priors1[1,0]) + priors1[1,0]
    delta1 = posterior1[:,4]*(priors1[2,1] - priors1[2,0]) + priors1[2,0]
    rbreak1 = posterior1[:,5]*(priors1[3,1] - priors1[3,0]) + priors1[3,0]
    cfact1 = posterior1[:,6]
    weights1 = posterior1[:,0]

    posterior2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_.txt')
    priors2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_priors.txt')
    gamma2 = posterior2[:,2]*(priors2[0,1] - priors2[0,0]) + priors2[0,0]
    alpha2 = posterior2[:,3]*(priors2[1,1] - priors2[1,0]) + priors2[1,0]
    delta2 = posterior2[:,4]*(priors2[2,1] - priors2[2,0]) + priors2[2,0]
    rbreak2 = posterior2[:,5]*(priors2[3,1] - priors2[3,0]) + priors2[3,0]
    cfact2 = posterior2[:,6]
    weights2 = posterior2[:,0]

    py.clf()
    py.hist(gamma1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1)
    py.hist(gamma2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed')
    py.legend(loc=2)
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_gamma.png')

    py.clf()
    py.hist(alpha1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1)
    py.hist(alpha2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed')
    py.legend(loc=4)
    py.xlabel(r'$\alpha$ (Outer Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_alpha.png')

    py.clf()
    py.hist(delta1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1)
    py.hist(delta2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed')
    py.legend()
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_delta.png')

    py.clf()
    py.hist(rbreak1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1)
    py.hist(rbreak2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed')
    py.legend(loc=2)
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_rbreak.png')

    py.clf()
    py.hist(cfact1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1)
    py.hist(cfact2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed')
    py.legend()
    py.xlabel(r'C Factor')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_cfact.png')



def compareTD(file='data_mthread_FINAL_free_maser_cuts_5.0_15.0arcsec',numBins=50):

    posterior1 = np.loadtxt('/u/schappell/pmnOld/'+file+'_post_equal_weights.dat')
    priors1 = np.loadtxt('/u/schappell/pmnOld/'+file+'_priors.txt')
    gamma1 = posterior1[:,0]*(priors1[0,1] - priors1[0,0]) + priors1[0,0]
    weights1 = posterior1[:,5]

    file2=np.loadtxt('/u/schappell/Downloads/gamma.txt')
    gamma2 = file2[:,0]
    hist = file2[:,1]
    tmpsum = 0
    for i in range(len(hist)-1):
        if ((gamma2[i] != gamma2[i+1]) & (hist[i] == hist[i+1])):
            tmpsum += hist[i] * (gamma2[i+1] - gamma2[i])

    py.clf()
    py.hist(gamma1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label='W/ accel')
    py.plot(gamma2,hist/tmpsum,'--',linewidth=2,label='Do 2013')
    py.legend(loc=2)
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTDO_gamma.png')
    pdb.set_trace()




def ar_r2d_sample_limits(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                         poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
                         Rcut=1.7,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
                         file1='data_oldC_oldDen_free_maser_cuts_5.0_15.0arcsec',label1='W/ Accel',
                         file2='data_no_accel_gcows_free_maser_cuts_5.0_15.0arcsec',label2='W/o Accel',numBins=50):



    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
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

    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

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

    #pdb.set_trace()
    for pp in range(len(r2d)):
        ar_cmss = -ar[pp] * 1.0
        are_cmss = are[pp] * 1.0
        r_cm = r2d[pp] * 1.0
        #pdb.set_trace()
        if (ar_cmss + (3.0*are_cmss) <= (GM/r_cm**2)):
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print (str(names[pp])+' has significant ar, '+str((ar_cmss+(3.0*are_cmss))*sec_in_yr/asy_to_kms/1e5)+
                   ' and lower a_z, '+str(az_cmss*sec_in_yr/asy_to_kms/1e5))


    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    #pdex = np.where(oldProb == 1.0)[0]
    num_stars = len(pdex)
    #tmprint = 0
    #print 'Stars used:'
    #for p in pdex:
     #   print tmprint+5,names[p],oldProb[p]
      #  tmprint += 1


    #Save info to file for c to read
    savenames = [names[pp] for pp in pdex]
    r2d = r2d[pdex] / cm_in_pc
    ar = ar[pdex]*sec_in_yr/1.0e5
    are = are[pdex]*sec_in_yr/1.0e5
    oldProb = oldProb[pdex]

    r2dline = np.arange(0.01*cm_in_pc,0.07*cm_in_pc,1e13)
    a0line = (GM/r2dline**2 * sec_in_yr / 1.0e5)

    reported = ['S0-26','S0-7','S0-4','S0-40','S0-1','S0-3','S0-19','S0-2','S0-16','S0-8','S0-20',
                'S0-103','S0-28','S0-27','S0-52','S0-61','S0-30','S0-38','S0-17', 'S0-68','S0-70',
                'S1-2','S1-3','S0-36','S1-12','S1-13','S1-3','S0-15','irs16C','irs16SW','S1-12',
                'S1-14','S1-3','S1-12']

    py.clf()
    py.subplot(211)

    for pp in range(len(r2d)):
        if (-ar[pp] - (0.5*are[pp]) > 0.0):
            new_source=1
            r_cm = r2d[pp] * dist * cm_in_au
            ar_cmss = ar[pp] * asy_to_kms * 1e5 / sec_in_yr
            are_cmss = are[pp] * asy_to_kms * 1e5 / sec_in_yr
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print str(savenames[pp])+' has significant ar, '+str(-ar[pp]-(3.0*are[pp]))
            #p1=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='b*',ms=14,capsize=4,elinewidth=2)
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    p2=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='go',ms=8,capsize=4,elinewidth=2)
                    new_source=0
            if (new_source==1):
                p1=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='b*',ms=14,capsize=4,elinewidth=2)

        else:
            print str(savenames[pp])+' has upperlimit, '+str(pp)
            upperlimit = -ar[pp]+3.0*are[pp]
            if (upperlimit <= 0.0):
                upperlimit = 0.0
            py.plot(r2d[pp],upperlimit,'k_',ms=10)
            py.arrow(r2d[pp],upperlimit,0,-1,hold=True,color='black',width=0.0003,
                     head_length=1,linewidth=0.0005,label='Upper Limits')
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    print 'have problem'
                    pdb.set_trace()

    py.plot(r2dline/cm_in_pc,a0line,'--k')
    py.legend([p1[0],p2[0]],['New Sources','Prev. Reported'],numpoints=1)
    py.text(0.017, 14, r'$|$a$|_{max} = \frac{GM}{R^2}$',fontsize=20)
    #py.yscale('log')
    py.xlim(0.01,0.07)
    py.ylim(-2,20)
    py.xlabel('Projected Radius (pc)')
    py.ylabel(r'$|$a$_R|$ (km/s/yr)')
    #py.title('Acceleration Detections & Upper Limits')

    posterior1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_post_equal_weights.dat')
    priors1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_priors.txt')
    gamma1 = posterior1[:,0]*(priors1[0,1] - priors1[0,0]) + priors1[0,0]
    alpha1 = posterior1[:,1]*(priors1[1,1] - priors1[1,0]) + priors1[1,0]
    delta1 = posterior1[:,2]*(priors1[2,1] - priors1[2,0]) + priors1[2,0]
    rbreak1 = posterior1[:,3]*(priors1[3,1] - priors1[3,0]) + priors1[3,0]
    cfact1 = posterior1[:,4]
    weights1 = posterior1[:,5]

    posterior2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_post_equal_weights.dat')
    priors2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_priors.txt')
    gamma2 = posterior2[:,0]*(priors2[0,1] - priors2[0,0]) + priors2[0,0]
    alpha2 = posterior2[:,1]*(priors2[1,1] - priors2[1,0]) + priors2[1,0]
    delta2 = posterior2[:,2]*(priors2[2,1] - priors2[2,0]) + priors2[2,0]
    rbreak2 = posterior2[:,3]*(priors2[3,1] - priors2[3,0]) + priors2[3,0]
    cfact2 = posterior2[:,4]
    weights2 = posterior2[:,5]

    py.subplot(212)
    py.hist(gamma1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1,color='k')
    py.hist(gamma2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed',color='k')
    py.legend(loc=2)
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/ar_r2d_sample_reported_limits_compare_gamma.png')
    py.clf()
    #pdb.set_trace()



def projectBPL(gamma=-1.0,alpha=4.0,delta=9.0,rbreak=0.5,min_R=0.01,max_R=0.75,
               numR=1000,maxZ=1.0,norm=1.0):
    #For given broken power law (alpha, gamma, delta, and r_break), return projected
    #density as a function of R_2D, between min R and max R given
    #numR is the number of R_2D positions desired
    #min and max R and max Z ALL in pc
    #RETURNS: arrays of R_2D and corresponding projected densities, each with length of numR

    dR = (max_R - min_R)/(numR - 1)
    R2d = np.array([min_R + i*dR for i in range(numR)])

    projDen = np.zeros(numR)

    for i in range(numR):
        projDen[i] = projectBPL_atR2D(R2d[i],gamma=gamma,alpha=alpha,delta=delta,rbreak=rbreak,maxZ=maxZ)
        #tmpInt = integrate.quad(lambda x: (math.sqrt(R2d[i]**2+x**2)/rbreak)**(-1.0*gamma)*
         #                       (1.0+(math.sqrt(R2d[i]**2+x**2)/rbreak)**delta)**((gamma-alpha)/delta),
          #                      -1.0*maxZ, maxZ)
        #projDen[i] = tmpInt[0]

    projDen *= norm / np.sum(projDen*dR)
    return R2d,projDen



def projectBPL_atR2D(R2d,gamma=-1.0,alpha=4.0,delta=9.0,rbreak=0.5,maxZ=1.0):
    #For given broken power law, finds projected surface density at given position R2D
    #R2d in pc
    #Rho_0 set to 1

    tmpInt = integrate.quad(lambda x: (math.sqrt(R2d**2+x**2)/rbreak)**(-1.0*gamma)*
                            (1.0+(math.sqrt(R2d**2+x**2)/rbreak)**delta)**((gamma-alpha)/delta),
                            -1.0*maxZ, maxZ)
    return tmpInt[0]



def projBPL_median_sigma(posteriorFile,priorsFile,maxZ=1.0,min_R=0.01,max_R=0.75,numR=1000,R2d_cut=5.0,norm=1.0):

    R2d_cut *= dist / au_in_pc

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]
    cfact_array = posterior[:,4]

    dR = (max_R - min_R)/(numR - 1)
    R2d = np.array([min_R + i*dR for i in range(numR)])
    median_rho = np.zeros(len(R2d))
    lower_rho = np.zeros(len(R2d))
    upper_rho = np.zeros(len(R2d))
    rho = np.zeros((len(gamma_array),len(R2d)))

    for i in range(len(R2d)):
        #tmp_rho = np.zeros(len(gamma_array))
        print 'Start R2d: '+str(R2d[i])+' pc, '+str(i)
        for j in range(len(gamma_array)):
            rho[j,i] = projectBPL_atR2D(R2d[i],gamma=gamma_array[j],alpha=alpha_array[j],
                                        delta=delta_array[j],rbreak=rbreak_array[j],maxZ=maxZ)
            if (R2d[i] < R2d_cut):
                rho[j,i] *= cfact_array[j]
            else:
                rho[j,i] *= (1.0 - cfact_array[j])
    for j in range(len(gamma_array)):
        rho[j,:] /= np.sum(rho[j,:])
    for i in range(len(R2d)):
        tmp_median = np.median(rho[:,i])
        median_rho[i] = np.median(rho[:,i])

        #tmpHist,tmpbins = np.histogram(tmp_rho,bins=50)
        #levels = getContourLevels(tmpHist)
        #tmpdex = np.where(tmpHist >= levels[0])[0]
        #upper_rho[i] = tmpbins[max(tmpdex)+1]
        #lower_rho[i] = tmpbins[min(tmpdex)]

        #pdb.set_trace()

        tmp_percent = 0.0
        delta_rho = tmp_median * 0.001
        upper_sigma = tmp_median * 1.0
        #print 'Found median, now find lower and upper sigma'
        while ((tmp_percent < 0.341) & (upper_sigma < np.max(rho[:,i]))):
            upper_sigma += delta_rho
            tmpdex = np.where((rho[:,i] >= tmp_median) & (rho[:,i] <= upper_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        tmp_percent = 0.0
        lower_sigma = tmp_median * 1.0
        while ((tmp_percent < 0.341) & (lower_sigma > np.min(rho[:,i]))):
            lower_sigma -= delta_rho
            tmpdex = np.where((rho[:,i] <= tmp_median) & (rho[:,i] >= lower_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        #pdb.set_trace()
        #lower_sigma = min(tmp_rho)
        #upper_sigma = max(tmp_rho)

        lower_rho[i] = lower_sigma * 1.0
        upper_rho[i] = upper_sigma * 1.0

    norm = norm / np.sum(median_rho*dR)
    median_rho *= norm
    lower_rho *= norm
    upper_rho *= norm
    pdb.set_trace()

    return R2d,median_rho,lower_rho,upper_rho






def PBPL_lsq(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
             poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
             Rcut=5.0,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
             max_r=5.0,schodel=False,maserDate='13jullgs1',ak_correct=True,isval=[-2.0,2.0,10],
             osval=[0.0,10.0,10],sval=[0.1,10.0,10],bval=[0.1,15.0,10],numBins=50,nonRadial=1,onlySchodel=False):

    if (schodel==True):
        ak_correct=True
    else:
        ak_correct=False #if using schodel's data, which is ext corrected, have to correct our data
        #otherwise, using stars from maser field, don't ext correct anything


    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    if (onlySchodel==False):
        names, r2d, ar, are = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                        Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                        chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)


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
    #ar = ar[pdex]
    #are = are[pdex]
        oldProb = oldProb[pdex]
    #np.savetxt('/u/schappell/code/c/stars_mn.dat',np.transpose([r2d,ar,are,oldProb]),delimiter=' ')

    #Rcut *= dist / au_in_pc #Rcut now in pc
    #pdb.set_trace()

    #mnRoot = '/u/schappell/pmnOld/data'+label
    #if ((situation == 2) | (situation == 4)):
    #    mnRoot += '_fixed_a'+str(mnAlpha)+'_d'+str(mnDelta)+'_br'+str(mnBreak)+'pc'
    #else:
    #    mnRoot += '_free'
    #if (situation > 2):
    #    mnRoot += '_maser_cuts_'+str(innerCut)+'_'+str(outerCut)+'arcsec'
    #else:
    #    mnRoot += '_cut_'+str(Rcut)+'arcsec'

    #mnRoot += '_'

    #Rcut *= dist / au_in_pc
    #if (situation > 2):
    #maserStars(maserDate=maserDate,innerCut=innerCut,outerCut=outerCut,magCut=magCut,schodel=schodel)
    #outerCut *= dist / au_in_pc
    #innerCut *= dist / au_in_pc

        maserInfo = np.loadtxt('/u/schappell/code/c/maser_mn.dat')

        r2d_m = maserInfo[:,0]
        oldProb_m = maserInfo[:,1]
        r2d /= dist * cm_in_au
        r2d_m /= dist * cm_in_au
        hist_m,bins_m = np.histogram(r2d_m,bins=numBins,weights=oldProb_m)
        sigma_m = np.sqrt(hist_m)
        middle_m = np.array([(bins_m[i+1] +bins_m[i])/2.0 for i in range(len(hist_m))])
        surface_area_m = np.array([bins_m[i+1]**2 - bins_m[i]**2 for i in range(len(hist_m))])*pi
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,0]
        gcows_cols = gcows[:,1]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        surface_area = np.zeros(len(hist))
        for i in range(len(hist)):
            gcowsdex = np.where((gcows_R > bins[i]) & (gcows_R <= bins[i+1]))[0]
            surface_area[i] = plate_scale**2 * len(gcowsdex)

        hist_m /= surface_area_m
        sigma_m /= surface_area_m

        hist,bins = np.histogram(r2d,bins=10,weights=oldProb)
        sigma = np.sqrt(hist)
        middle = np.array([(bins[i+1] +bins[i])/2.0 for i in range(len(hist))])

    else:
        schodelInfo = np.loadtxt('/u/schappell/code/c/onlySchodel_mn.dat')
        r2d = schodelInfo[:,0] / dist / cm_in_au
        oldProb = schodelInfo[:,1]
        hist,bins = np.histogram(r2d,bins=numBins,weights=oldProb)
        sigma = np.sqrt(hist)
        middle = np.array([(bins[i+1] +bins[i])/2.0 for i in range(len(hist))])
        surface_area = np.array([bins[i+1]**2 - bins[i]**2 for i in range(len(hist))])*pi

    hist /= surface_area
    sigma /= surface_area

    if (onlySchodel==False):
        middle = np.append(middle,middle_m)
        hist = np.append(hist,hist_m)
        sigma = np.append(sigma,sigma_m)

    sq_diff = np.zeros((isval[2],osval[2],sval[2],bval[2]))
    inner_slope = np.array([isval[0]+((isval[1]-isval[0])*i/(isval[2]-1)) for i in range(isval[2])])
    outer_slope = np.array([osval[0]+((osval[1]-osval[0])*i/(osval[2]-1)) for i in range(osval[2])])
    sharpness = np.array([sval[0]+((sval[1]-sval[0])*i/(sval[2]-1)) for i in range(sval[2])])
    break_R = np.array([bval[0]+((bval[1]-bval[0])*i/(bval[2]-1)) for i in range(bval[2])])

    for isi in range(isval[2]):
        is_tmp = inner_slope[isi]
        for osi in range(osval[2]):
            os_tmp = outer_slope[osi]
            for si in range(sval[2]):
                s_tmp = sharpness[si]
                for bi in range(bval[2]):
                    b_tmp = break_R[bi]
                    model = (middle/b_tmp)**(-1.0*is_tmp) * (1.0 + (middle/b_tmp)**s_tmp)**((is_tmp-os_tmp)/s_tmp)

                    model *= np.sum(hist) / np.sum(model)
                    #if ((np.sum(model))!=(np.sum(model))):
                     #   pdb.set_trace()
                    sq_diff[isi,osi,si,bi] = np.sum((hist - model)**2/sigma**2)

    py.clf()
    pdb.set_trace()






def maserStars_out(maserDate='13jullgs1',innerCut=5.0,outerCut=15.0,magCut=15.5,schodel=False,lmagCut=0.0,
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

    return R2d[mdex],oldProb[mdex],mag[mdex]




def GCOWS_schodel_overlap(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                          poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',Rcut=5.0,
                          nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0):

    ak_correct = True
    nonRadial = 1

    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                         Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                         chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)


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
 

    #Save info to file for c to read
    r2d = r2d[pdex] / dist / cm_in_au
#    ar = ar[pdex]
 #   are = are[pdex]
    oldProb = oldProb[pdex]



    #now schodel stars
    R2d_schodel = np.array([])
    mag_schodel = np.array([])
    x_schodel = np.array([])
    y_schodel = np.array([])

    ext_scale = 0.02705 #arcsec/pixel
    ext_center = [808,611]
    extMap = pyfits.getdata('/u/schappell/Downloads/AKs_fg6.fits')
    
    dbfile = '/g/ghez/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

    cur.execute('SELECT r2d,k,x,y FROM schoedel2010')
    for row in cur:
        Xext = int(round(row[2]/ext_scale))+ext_center[0]
        Yext = int(round(row[3]/ext_scale))+ext_center[0]
        if ((Xext < 0) | (Yext < 0)):
            print 'Something is wrong, star is calculated as being off extinction map'
            return
        else:
            if ((Xext < 1600) & (Yext < 1600)):
                R2d_schodel = np.append(R2d_schodel,row[0])
                mag_schodel = np.append(mag_schodel,row[1] + 2.7 - extMap[Xext,Yext])
                x_schodel = np.append(x_schodel,row[2])
                y_schodel = np.append(y_schodel,row[3])

    in_gcows = np.zeros(len(mag_schodel))
    plate_scale = 0.00995
    gcows_zero = np.array([1500.0,1500.0])
    gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')
    gcows_dex = np.where(gcows > 0)
        #saves x and y positions, in pixels
    for ii in range(len(mag_schodel)):
            #pdb.set_trace()
        try:
            x_pixels = int(round((-1.0*x_schodel[ii] / plate_scale) + gcows_zero[0]))
            y_pixels = int(round((y_schodel[ii] / plate_scale) + gcows_zero[1]))
            in_gcows[ii] = gcows[y_pixels,x_pixels]
        except:
            continue
        if ((x_pixels < 0.0) | (y_pixels < 0.0)):
                in_gcows[ii] = 0
    idx = np.where((in_gcows==1) & (R2d_schodel < Rcut))[0]

    R2d_schodel = R2d_schodel[idx]
    mag_schodel = mag_schodel[idx]
    x_schodel = x_schodel[idx]
    y_schodel = y_schodel[idx]

    gcows_dex = np.where(gcows > 0)

    schodelHist,bins = np.histogram(R2d_schodel)
    gcowsHist,gcows_bins = np.histogram(r2d,bins=bins)

    gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
    gcows_rows = gcows_dex[0]
    gcows_cols = gcows_dex[1]
    gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
    gcows_surface = np.zeros(len(gcowsHist))
    for i in range(len(schodelHist)):
        gcowsdex = np.where((gcows_R > bins[i]) & (gcows_R <= bins[i+1]))[0]
        gcows_surface[i] = plate_scale**2 * len(gcowsdex)

    schodel_middle = np.array([(bins[i+1]+bins[i])/2.0 for i in range(len(bins)-1)])
    dex = range(len(gcowsHist))

    py.clf()
    py.errorbar(schodel_middle[dex],gcowsHist/gcows_surface[dex],yerr=np.sqrt(gcowsHist)/gcows_surface[dex],fmt='ok')
    py.errorbar(schodel_middle,schodelHist/gcows_surface,yerr=np.sqrt(schodelHist)/gcows_surface,fmt='ob')
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/overlap_surface_density.png')

    py.clf()
    bins = np.array([6.0+0.5*i for i in range(29)])
    py.hist(mag_schodel,bins=bins,label='Schodel',histtype='step',linewidth=2.0)
    py.hist(mag,bins=bins,label='GCOWS',histtype='step',linewidth=2.0)
    py.legend(loc=2)
    py.xlabel("K' (mag)")
    py.ylabel('Number stars')
    py.savefig('/u/schappell/plots/overlap_KLF.png')





def IAU16_proc(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
               poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
               Rcut=1.7,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
               file1='data_oldC_oldDen_free_maser_cuts_5.0_15.0arcsec',label1='W/ Accel',
               file2='data_no_accel_gcows_free_maser_cuts_5.0_15.0arcsec',label2='W/o Accel',numBins=50):

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
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

    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

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

    #pdb.set_trace()
    for pp in range(len(r2d)):
        ar_cmss = -ar[pp] * 1.0
        are_cmss = are[pp] * 1.0
        r_cm = r2d[pp] * 1.0
        #pdb.set_trace()
        if (ar_cmss + (3.0*are_cmss) <= (GM/r_cm**2)):
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print (str(names[pp])+' has significant ar, '+str((ar_cmss+(3.0*are_cmss))*sec_in_yr/asy_to_kms/1e5)+
                   ' and lower a_z, '+str(az_cmss*sec_in_yr/asy_to_kms/1e5))


    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    #pdex = np.where(oldProb == 1.0)[0]
    num_stars = len(pdex)
    #tmprint = 0
    #print 'Stars used:'
    #for p in pdex:
     #   print tmprint+5,names[p],oldProb[p]
      #  tmprint += 1


    #Save info to file for c to read
    savenames = [names[pp] for pp in pdex]
    r2d = r2d[pdex] / cm_in_pc
    ar = ar[pdex]*sec_in_yr/1.0e5
    are = are[pdex]*sec_in_yr/1.0e5
    oldProb = oldProb[pdex]

    r2dline = np.arange(0.01*cm_in_pc,0.07*cm_in_pc,1e13)
    a0line = (GM/r2dline**2 * sec_in_yr / 1.0e5)

    reported = ['S0-26','S0-7','S0-4','S0-40','S0-1','S0-3','S0-19','S0-2','S0-16','S0-8','S0-20',
                'S0-103','S0-28','S0-27','S0-52','S0-61','S0-30','S0-38','S0-17', 'S0-68','S0-70',
                'S1-2','S1-3','S0-36','S1-12','S1-13','S1-3','S0-15','irs16C','irs16SW','S1-12',
                'S1-14','S1-3','S1-12']

    py.clf()
    py.subplot(211)
    for pp in range(len(r2d)):
        if (-ar[pp] - (1.0*are[pp]) > 0.0):
            r_cm = r2d[pp] * dist * cm_in_au
            ar_cmss = ar[pp] * asy_to_kms * 1e5 / sec_in_yr
            are_cmss = are[pp] * asy_to_kms * 1e5 / sec_in_yr
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print str(savenames[pp])+' has significant ar, '+str(-ar[pp]-(3.0*are[pp]))
            p1=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='b*',ms=14,capsize=4,elinewidth=2)
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    p2=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='gd',ms=14,capsize=4,elinewidth=2)

        else:
            print str(savenames[pp])+' has upperlimit, '+str(pp)
            upperlimit = -ar[pp]+3.0*are[pp]
            if (upperlimit <= 0.0):
                upperlimit = 0.0
            py.plot(r2d[pp],upperlimit,'k_',ms=10)
            py.arrow(r2d[pp],upperlimit,0,-1,hold=True,color='black',width=0.0003,
                     head_length=1,linewidth=0.0005,label='Upper Limits')
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    print 'have problem'
                    pdb.set_trace()

    py.plot(r2dline/cm_in_pc,a0line,'--k')
    py.legend([p1[0],p2[0]],['New Sources','Prev. Reported'],numpoints=1,loc=2)
    py.text(0.017, 14, r'$|$a$|_{max} = \frac{GM}{R^2}$',fontsize=20)
    #py.yscale('log')
    py.xlim(0.01,0.07)
    py.ylim(-2,20)
    py.xlabel('Projected Radius (pc)')
    py.ylabel(r'$|$a$_R|$ (km/s/yr)')
    #py.title('Acceleration Detections & Upper Limits')



    posterior1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_post_equal_weights.dat')
    priors1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_priors.txt')
    gamma1 = posterior1[:,0]*(priors1[0,1] - priors1[0,0]) + priors1[0,0]
    alpha1 = posterior1[:,1]*(priors1[1,1] - priors1[1,0]) + priors1[1,0]
    delta1 = posterior1[:,2]*(priors1[2,1] - priors1[2,0]) + priors1[2,0]
    rbreak1 = posterior1[:,3]*(priors1[3,1] - priors1[3,0]) + priors1[3,0]
    cfact1 = posterior1[:,4]
    weights1 = posterior1[:,5]

    posterior2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_post_equal_weights.dat')
    priors2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_priors.txt')
    gamma2 = posterior2[:,0]*(priors2[0,1] - priors2[0,0]) + priors2[0,0]
    alpha2 = posterior2[:,1]*(priors2[1,1] - priors2[1,0]) + priors2[1,0]
    delta2 = posterior2[:,2]*(priors2[2,1] - priors2[2,0]) + priors2[2,0]
    rbreak2 = posterior2[:,3]*(priors2[3,1] - priors2[3,0]) + priors2[3,0]
    cfact2 = posterior2[:,4]
    weights2 = posterior2[:,5]

    py.subplot(212)
    py.hist(gamma1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1,color='k')
    py.hist(gamma2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed',color='k')
    py.legend()
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_gamma_AND_ar_vs_r2d.png')






def mass_in_radius(posteriorFile,priorsFile,maxZ=5.0,min_R=0.0001,max_R=1.0,numR=1000,mass_max=1e6,m_radius=0.01):

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]

    dR = (max_R - min_R)/(numR - 1)
    R2d = np.array([min_R + i*dR for i in range(numR)])

    rho = np.zeros((len(gamma_array),len(R2d)))
    mass_out = np.zeros(len(gamma_array))

    for i in range(len(R2d)):
        print 'Start R2d: '+str(R2d[i])+' pc, '+str(i)
        for j in range(len(gamma_array)):
            rho[j,i] = projectBPL_atR2D(R2d[i],gamma=gamma_array[j],alpha=alpha_array[j],
                                        delta=delta_array[j],rbreak=rbreak_array[j],maxZ=maxZ)
            rho[j,i] *= R2d[i] * dR




def projBPL_median_sigma(posteriorFile,priorsFile,maxZ=1.0,min_R=0.01,max_R=0.75,numR=1000,R2d_cut=5.0,norm=1.0):

    R2d_cut *= dist / au_in_pc

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]
    cfact_array = posterior[:,4]

    dR = (max_R - min_R)/(numR - 1)
    R2d = np.array([min_R + i*dR for i in range(numR)])
    median_rho = np.zeros(len(R2d))
    lower_rho = np.zeros(len(R2d))
    upper_rho = np.zeros(len(R2d))
    rho = np.zeros((len(gamma_array),len(R2d)))

    for i in range(len(R2d)):
        #tmp_rho = np.zeros(len(gamma_array))
        print 'Start R2d: '+str(R2d[i])+' pc, '+str(i)
        for j in range(len(gamma_array)):
            rho[j,i] = projectBPL_atR2D(R2d[i],gamma=gamma_array[j],alpha=alpha_array[j],
                                        delta=delta_array[j],rbreak=rbreak_array[j],maxZ=maxZ)
            if (R2d[i] < R2d_cut):
                rho[j,i] *= cfact_array[j]
            else:
                rho[j,i] *= (1.0 - cfact_array[j])
    for j in range(len(gamma_array)):
        rho[j,:] /= np.sum(rho[j,:])
    for i in range(len(R2d)):
        tmp_median = np.median(rho[:,i])
        median_rho[i] = np.median(rho[:,i])

        #tmpHist,tmpbins = np.histogram(tmp_rho,bins=50)
        #levels = getContourLevels(tmpHist)
        #tmpdex = np.where(tmpHist >= levels[0])[0]
        #upper_rho[i] = tmpbins[max(tmpdex)+1]
        #lower_rho[i] = tmpbins[min(tmpdex)]

        #pdb.set_trace()

        tmp_percent = 0.0
        delta_rho = tmp_median * 0.001
        upper_sigma = tmp_median * 1.0
        #print 'Found median, now find lower and upper sigma'
        while ((tmp_percent < 0.341) & (upper_sigma < np.max(rho[:,i]))):
            upper_sigma += delta_rho
            tmpdex = np.where((rho[:,i] >= tmp_median) & (rho[:,i] <= upper_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        tmp_percent = 0.0
        lower_sigma = tmp_median * 1.0
        while ((tmp_percent < 0.341) & (lower_sigma > np.min(rho[:,i]))):
            lower_sigma -= delta_rho
            tmpdex = np.where((rho[:,i] <= tmp_median) & (rho[:,i] >= lower_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        #pdb.set_trace()
        #lower_sigma = min(tmp_rho)
        #upper_sigma = max(tmp_rho)

        lower_rho[i] = lower_sigma * 1.0
        upper_rho[i] = upper_sigma * 1.0

    norm = norm / np.sum(median_rho*dR)
    median_rho *= norm
    lower_rho *= norm
    upper_rho *= norm
    pdb.set_trace()

    return R2d,median_rho,lower_rho,upper_rho






def PBPL_lsq(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
             poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
             Rcut=5.0,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
             max_r=5.0,schodel=False,maserDate='13jullgs1',ak_correct=True,isval=[-2.0,2.0,10],
             osval=[0.0,10.0,10],sval=[0.1,10.0,10],bval=[0.1,15.0,10],numBins=50,nonRadial=1,onlySchodel=False):

    if (schodel==True):
        ak_correct=True
    else:
        ak_correct=False #if using schodel's data, which is ext corrected, have to correct our data
        #otherwise, using stars from maser field, don't ext correct anything


    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    if (onlySchodel==False):
        names, r2d, ar, are = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                        Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                        chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)


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
    #ar = ar[pdex]
    #are = are[pdex]
        oldProb = oldProb[pdex]
    #np.savetxt('/u/schappell/code/c/stars_mn.dat',np.transpose([r2d,ar,are,oldProb]),delimiter=' ')

    #Rcut *= dist / au_in_pc #Rcut now in pc
    #pdb.set_trace()

    #mnRoot = '/u/schappell/pmnOld/data'+label
    #if ((situation == 2) | (situation == 4)):
    #    mnRoot += '_fixed_a'+str(mnAlpha)+'_d'+str(mnDelta)+'_br'+str(mnBreak)+'pc'
    #else:
    #    mnRoot += '_free'
    #if (situation > 2):
    #    mnRoot += '_maser_cuts_'+str(innerCut)+'_'+str(outerCut)+'arcsec'
    #else:
    #    mnRoot += '_cut_'+str(Rcut)+'arcsec'

    #mnRoot += '_'

    #Rcut *= dist / au_in_pc
    #if (situation > 2):
    #maserStars(maserDate=maserDate,innerCut=innerCut,outerCut=outerCut,magCut=magCut,schodel=schodel)
    #outerCut *= dist / au_in_pc
    #innerCut *= dist / au_in_pc

        maserInfo = np.loadtxt('/u/schappell/code/c/maser_mn.dat')

        r2d_m = maserInfo[:,0]
        oldProb_m = maserInfo[:,1]
        r2d /= dist * cm_in_au
        r2d_m /= dist * cm_in_au
        hist_m,bins_m = np.histogram(r2d_m,bins=numBins,weights=oldProb_m)
        sigma_m = np.sqrt(hist_m)
        middle_m = np.array([(bins_m[i+1] +bins_m[i])/2.0 for i in range(len(hist_m))])
        surface_area_m = np.array([bins_m[i+1]**2 - bins_m[i]**2 for i in range(len(hist_m))])*pi
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,0]
        gcows_cols = gcows[:,1]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        surface_area = np.zeros(len(hist))
        for i in range(len(hist)):
            gcowsdex = np.where((gcows_R > bins[i]) & (gcows_R <= bins[i+1]))[0]
            surface_area[i] = plate_scale**2 * len(gcowsdex)

        hist_m /= surface_area_m
        sigma_m /= surface_area_m

        hist,bins = np.histogram(r2d,bins=10,weights=oldProb)
        sigma = np.sqrt(hist)
        middle = np.array([(bins[i+1] +bins[i])/2.0 for i in range(len(hist))])

    else:
        schodelInfo = np.loadtxt('/u/schappell/code/c/onlySchodel_mn.dat')
        r2d = schodelInfo[:,0] / dist / cm_in_au
        oldProb = schodelInfo[:,1]
        hist,bins = np.histogram(r2d,bins=numBins,weights=oldProb)
        sigma = np.sqrt(hist)
        middle = np.array([(bins[i+1] +bins[i])/2.0 for i in range(len(hist))])
        surface_area = np.array([bins[i+1]**2 - bins[i]**2 for i in range(len(hist))])*pi

    hist /= surface_area
    sigma /= surface_area

    if (onlySchodel==False):
        middle = np.append(middle,middle_m)
        hist = np.append(hist,hist_m)
        sigma = np.append(sigma,sigma_m)

    sq_diff = np.zeros((isval[2],osval[2],sval[2],bval[2]))
    inner_slope = np.array([isval[0]+((isval[1]-isval[0])*i/(isval[2]-1)) for i in range(isval[2])])
    outer_slope = np.array([osval[0]+((osval[1]-osval[0])*i/(osval[2]-1)) for i in range(osval[2])])
    sharpness = np.array([sval[0]+((sval[1]-sval[0])*i/(sval[2]-1)) for i in range(sval[2])])
    break_R = np.array([bval[0]+((bval[1]-bval[0])*i/(bval[2]-1)) for i in range(bval[2])])

    for isi in range(isval[2]):
        is_tmp = inner_slope[isi]
        for osi in range(osval[2]):
            os_tmp = outer_slope[osi]
            for si in range(sval[2]):
                s_tmp = sharpness[si]
                for bi in range(bval[2]):
                    b_tmp = break_R[bi]
                    model = (middle/b_tmp)**(-1.0*is_tmp) * (1.0 + (middle/b_tmp)**s_tmp)**((is_tmp-os_tmp)/s_tmp)

                    model *= np.sum(hist) / np.sum(model)
                    #if ((np.sum(model))!=(np.sum(model))):
                     #   pdb.set_trace()
                    sq_diff[isi,osi,si,bi] = np.sum((hist - model)**2/sigma**2)

    py.clf()
    pdb.set_trace()






def maserStars_out(maserDate='13jullgs1',innerCut=5.0,outerCut=15.0,magCut=15.5,schodel=False,lmagCut=0.0,
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

    return R2d[mdex],oldProb[mdex],mag[mdex]




def GCOWS_schodel_overlap(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                          poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',Rcut=5.0,
                          nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0):

    ak_correct = True
    nonRadial = 1

    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                         Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                         chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)


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
 

    #Save info to file for c to read
    r2d = r2d[pdex] / dist / cm_in_au
#    ar = ar[pdex]
 #   are = are[pdex]
    oldProb = oldProb[pdex]



    #now schodel stars
    R2d_schodel = np.array([])
    mag_schodel = np.array([])
    x_schodel = np.array([])
    y_schodel = np.array([])

    ext_scale = 0.02705 #arcsec/pixel
    ext_center = [808,611]
    extMap = pyfits.getdata('/u/schappell/Downloads/AKs_fg6.fits')
    
    dbfile = '/g/ghez/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

    cur.execute('SELECT r2d,k,x,y FROM schoedel2010')
    for row in cur:
        Xext = int(round(row[2]/ext_scale))+ext_center[0]
        Yext = int(round(row[3]/ext_scale))+ext_center[0]
        if ((Xext < 0) | (Yext < 0)):
            print 'Something is wrong, star is calculated as being off extinction map'
            return
        else:
            if ((Xext < 1600) & (Yext < 1600)):
                R2d_schodel = np.append(R2d_schodel,row[0])
                mag_schodel = np.append(mag_schodel,row[1] + 2.7 - extMap[Xext,Yext])
                x_schodel = np.append(x_schodel,row[2])
                y_schodel = np.append(y_schodel,row[3])

    in_gcows = np.zeros(len(mag_schodel))
    plate_scale = 0.00995
    gcows_zero = np.array([1500.0,1500.0])
    gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')
    gcows_dex = np.where(gcows > 0)
        #saves x and y positions, in pixels
    for ii in range(len(mag_schodel)):
            #pdb.set_trace()
        try:
            x_pixels = int(round((-1.0*x_schodel[ii] / plate_scale) + gcows_zero[0]))
            y_pixels = int(round((y_schodel[ii] / plate_scale) + gcows_zero[1]))
            in_gcows[ii] = gcows[y_pixels,x_pixels]
        except:
            continue
        if ((x_pixels < 0.0) | (y_pixels < 0.0)):
                in_gcows[ii] = 0
    idx = np.where((in_gcows==1) & (R2d_schodel < Rcut))[0]

    R2d_schodel = R2d_schodel[idx]
    mag_schodel = mag_schodel[idx]
    x_schodel = x_schodel[idx]
    y_schodel = y_schodel[idx]

    gcows_dex = np.where(gcows > 0)

    schodelHist,bins = np.histogram(R2d_schodel)
    gcowsHist,gcows_bins = np.histogram(r2d,bins=bins)

    gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
    gcows_rows = gcows_dex[0]
    gcows_cols = gcows_dex[1]
    gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
    gcows_surface = np.zeros(len(gcowsHist))
    for i in range(len(schodelHist)):
        gcowsdex = np.where((gcows_R > bins[i]) & (gcows_R <= bins[i+1]))[0]
        gcows_surface[i] = plate_scale**2 * len(gcowsdex)

    schodel_middle = np.array([(bins[i+1]+bins[i])/2.0 for i in range(len(bins)-1)])
    dex = range(len(gcowsHist))

    py.clf()
    py.errorbar(schodel_middle[dex],gcowsHist/gcows_surface[dex],yerr=np.sqrt(gcowsHist)/gcows_surface[dex],fmt='ok')
    py.errorbar(schodel_middle,schodelHist/gcows_surface,yerr=np.sqrt(schodelHist)/gcows_surface,fmt='ob')
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/overlap_surface_density.png')

    py.clf()
    bins = np.array([6.0+0.5*i for i in range(29)])
    py.hist(mag_schodel,bins=bins,label='Schodel',histtype='step',linewidth=2.0)
    py.hist(mag,bins=bins,label='GCOWS',histtype='step',linewidth=2.0)
    py.legend(loc=2)
    py.xlabel("K' (mag)")
    py.ylabel('Number stars')
    py.savefig('/u/schappell/plots/overlap_KLF.png')





def IAU16_proc(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
               poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
               Rcut=1.7,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
               file1='data_oldC_oldDen_free_maser_cuts_5.0_15.0arcsec',label1='W/ Accel',
               file2='data_no_accel_gcows_free_maser_cuts_5.0_15.0arcsec',label2='W/o Accel',numBins=50):

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
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

    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

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

    #pdb.set_trace()
    for pp in range(len(r2d)):
        ar_cmss = -ar[pp] * 1.0
        are_cmss = are[pp] * 1.0
        r_cm = r2d[pp] * 1.0
        #pdb.set_trace()
        if (ar_cmss + (3.0*are_cmss) <= (GM/r_cm**2)):
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print (str(names[pp])+' has significant ar, '+str((ar_cmss+(3.0*are_cmss))*sec_in_yr/asy_to_kms/1e5)+
                   ' and lower a_z, '+str(az_cmss*sec_in_yr/asy_to_kms/1e5))


    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    #pdex = np.where(oldProb == 1.0)[0]
    num_stars = len(pdex)
    #tmprint = 0
    #print 'Stars used:'
    #for p in pdex:
     #   print tmprint+5,names[p],oldProb[p]
      #  tmprint += 1


    #Save info to file for c to read
    savenames = [names[pp] for pp in pdex]
    r2d = r2d[pdex] / cm_in_pc
    ar = ar[pdex]*sec_in_yr/1.0e5
    are = are[pdex]*sec_in_yr/1.0e5
    oldProb = oldProb[pdex]

    r2dline = np.arange(0.01*cm_in_pc,0.07*cm_in_pc,1e13)
    a0line = (GM/r2dline**2 * sec_in_yr / 1.0e5)

    reported = ['S0-26','S0-7','S0-4','S0-40','S0-1','S0-3','S0-19','S0-2','S0-16','S0-8','S0-20',
                'S0-103','S0-28','S0-27','S0-52','S0-61','S0-30','S0-38','S0-17', 'S0-68','S0-70',
                'S1-2','S1-3','S0-36','S1-12','S1-13','S1-3','S0-15','irs16C','irs16SW','S1-12',
                'S1-14','S1-3','S1-12']

    py.clf()
    py.subplot(211)
    for pp in range(len(r2d)):
        if (-ar[pp] - (1.0*are[pp]) > 0.0):
            r_cm = r2d[pp] * dist * cm_in_au
            ar_cmss = ar[pp] * asy_to_kms * 1e5 / sec_in_yr
            are_cmss = are[pp] * asy_to_kms * 1e5 / sec_in_yr
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print str(savenames[pp])+' has significant ar, '+str(-ar[pp]-(3.0*are[pp]))
            p1=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='b*',ms=14,capsize=4,elinewidth=2)
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    p2=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='gd',ms=14,capsize=4,elinewidth=2)

        else:
            print str(savenames[pp])+' has upperlimit, '+str(pp)
            upperlimit = -ar[pp]+3.0*are[pp]
            if (upperlimit <= 0.0):
                upperlimit = 0.0
            py.plot(r2d[pp],upperlimit,'k_',ms=10)
            py.arrow(r2d[pp],upperlimit,0,-1,hold=True,color='black',width=0.0003,
                     head_length=1,linewidth=0.0005,label='Upper Limits')
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    print 'have problem'
                    pdb.set_trace()

    py.plot(r2dline/cm_in_pc,a0line,'--k')
    py.legend([p1[0],p2[0]],['New Sources','Prev. Reported'],numpoints=1,loc=2)
    py.text(0.017, 14, r'$|$a$|_{max} = \frac{GM}{R^2}$',fontsize=20)
    #py.yscale('log')
    py.xlim(0.01,0.07)
    py.ylim(-2,20)
    py.xlabel('Projected Radius (pc)')
    py.ylabel(r'$|$a$_R|$ (km/s/yr)')
    #py.title('Acceleration Detections & Upper Limits')



    posterior1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_post_equal_weights.dat')
    priors1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_priors.txt')
    gamma1 = posterior1[:,0]*(priors1[0,1] - priors1[0,0]) + priors1[0,0]
    alpha1 = posterior1[:,1]*(priors1[1,1] - priors1[1,0]) + priors1[1,0]
    delta1 = posterior1[:,2]*(priors1[2,1] - priors1[2,0]) + priors1[2,0]
    rbreak1 = posterior1[:,3]*(priors1[3,1] - priors1[3,0]) + priors1[3,0]
    cfact1 = posterior1[:,4]
    weights1 = posterior1[:,5]

    posterior2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_post_equal_weights.dat')
    priors2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_priors.txt')
    gamma2 = posterior2[:,0]*(priors2[0,1] - priors2[0,0]) + priors2[0,0]
    alpha2 = posterior2[:,1]*(priors2[1,1] - priors2[1,0]) + priors2[1,0]
    delta2 = posterior2[:,2]*(priors2[2,1] - priors2[2,0]) + priors2[2,0]
    rbreak2 = posterior2[:,3]*(priors2[3,1] - priors2[3,0]) + priors2[3,0]
    cfact2 = posterior2[:,4]
    weights2 = posterior2[:,5]

    py.subplot(212)
    py.hist(gamma1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1,color='k')
    py.hist(gamma2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed',color='k')
    py.legend()
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_gamma_AND_ar_vs_r2d.png')






def mass_in_radius(posteriorFile,priorsFile,maxZ=5.0,min_R=0.0001,max_R=1.0,numR=1000,mass_max=1e6,m_radius=0.01):

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]

    dR = (max_R - min_R)/(numR - 1)
    R2d = np.array([min_R + i*dR for i in range(numR)])

    rho = np.zeros((len(gamma_array),len(R2d)))
    mass_out = np.zeros(len(gamma_array))

    for i in range(len(R2d)):
        print 'Start R2d: '+str(R2d[i])+' pc, '+str(i)
        for j in range(len(gamma_array)):
            rho[j,i] = projectBPL_atR2D(R2d[i],gamma=gamma_array[j],alpha=alpha_array[j],
                                        delta=delta_array[j],rbreak=rbreak_array[j],maxZ=maxZ)
            rho[j,i] *= R2d[i] * dR




def projBPL_median_sigma(posteriorFile,priorsFile,maxZ=1.0,min_R=0.01,max_R=0.75,numR=1000,R2d_cut=5.0,norm=1.0):

    R2d_cut *= dist / au_in_pc

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]
    cfact_array = posterior[:,4]

    dR = (max_R - min_R)/(numR - 1)
    R2d = np.array([min_R + i*dR for i in range(numR)])
    median_rho = np.zeros(len(R2d))
    lower_rho = np.zeros(len(R2d))
    upper_rho = np.zeros(len(R2d))
    rho = np.zeros((len(gamma_array),len(R2d)))

    for i in range(len(R2d)):
        #tmp_rho = np.zeros(len(gamma_array))
        print 'Start R2d: '+str(R2d[i])+' pc, '+str(i)
        for j in range(len(gamma_array)):
            rho[j,i] = projectBPL_atR2D(R2d[i],gamma=gamma_array[j],alpha=alpha_array[j],
                                        delta=delta_array[j],rbreak=rbreak_array[j],maxZ=maxZ)
            if (R2d[i] < R2d_cut):
                rho[j,i] *= cfact_array[j]
            else:
                rho[j,i] *= (1.0 - cfact_array[j])
    for j in range(len(gamma_array)):
        rho[j,:] /= np.sum(rho[j,:])
    for i in range(len(R2d)):
        tmp_median = np.median(rho[:,i])
        median_rho[i] = np.median(rho[:,i])

        #tmpHist,tmpbins = np.histogram(tmp_rho,bins=50)
        #levels = getContourLevels(tmpHist)
        #tmpdex = np.where(tmpHist >= levels[0])[0]
        #upper_rho[i] = tmpbins[max(tmpdex)+1]
        #lower_rho[i] = tmpbins[min(tmpdex)]

        #pdb.set_trace()

        tmp_percent = 0.0
        delta_rho = tmp_median * 0.001
        upper_sigma = tmp_median * 1.0
        #print 'Found median, now find lower and upper sigma'
        while ((tmp_percent < 0.341) & (upper_sigma < np.max(rho[:,i]))):
            upper_sigma += delta_rho
            tmpdex = np.where((rho[:,i] >= tmp_median) & (rho[:,i] <= upper_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        tmp_percent = 0.0
        lower_sigma = tmp_median * 1.0
        while ((tmp_percent < 0.341) & (lower_sigma > np.min(rho[:,i]))):
            lower_sigma -= delta_rho
            tmpdex = np.where((rho[:,i] <= tmp_median) & (rho[:,i] >= lower_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        #pdb.set_trace()
        #lower_sigma = min(tmp_rho)
        #upper_sigma = max(tmp_rho)

        lower_rho[i] = lower_sigma * 1.0
        upper_rho[i] = upper_sigma * 1.0

    norm = norm / np.sum(median_rho*dR)
    median_rho *= norm
    lower_rho *= norm
    upper_rho *= norm
    pdb.set_trace()

    return R2d,median_rho,lower_rho,upper_rho






def PBPL_lsq(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
             poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
             Rcut=5.0,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
             max_r=5.0,schodel=False,maserDate='13jullgs1',ak_correct=True,isval=[-2.0,2.0,10],
             osval=[0.0,10.0,10],sval=[0.1,10.0,10],bval=[0.1,15.0,10],numBins=50,nonRadial=1,onlySchodel=False):

    if (schodel==True):
        ak_correct=True
    else:
        ak_correct=False #if using schodel's data, which is ext corrected, have to correct our data
        #otherwise, using stars from maser field, don't ext correct anything


    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    if (onlySchodel==False):
        names, r2d, ar, are = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                        Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                        chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)


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
    #ar = ar[pdex]
    #are = are[pdex]
        oldProb = oldProb[pdex]
    #np.savetxt('/u/schappell/code/c/stars_mn.dat',np.transpose([r2d,ar,are,oldProb]),delimiter=' ')

    #Rcut *= dist / au_in_pc #Rcut now in pc
    #pdb.set_trace()

    #mnRoot = '/u/schappell/pmnOld/data'+label
    #if ((situation == 2) | (situation == 4)):
    #    mnRoot += '_fixed_a'+str(mnAlpha)+'_d'+str(mnDelta)+'_br'+str(mnBreak)+'pc'
    #else:
    #    mnRoot += '_free'
    #if (situation > 2):
    #    mnRoot += '_maser_cuts_'+str(innerCut)+'_'+str(outerCut)+'arcsec'
    #else:
    #    mnRoot += '_cut_'+str(Rcut)+'arcsec'

    #mnRoot += '_'

    #Rcut *= dist / au_in_pc
    #if (situation > 2):
    #maserStars(maserDate=maserDate,innerCut=innerCut,outerCut=outerCut,magCut=magCut,schodel=schodel)
    #outerCut *= dist / au_in_pc
    #innerCut *= dist / au_in_pc

        maserInfo = np.loadtxt('/u/schappell/code/c/maser_mn.dat')

        r2d_m = maserInfo[:,0]
        oldProb_m = maserInfo[:,1]
        r2d /= dist * cm_in_au
        r2d_m /= dist * cm_in_au
        hist_m,bins_m = np.histogram(r2d_m,bins=numBins,weights=oldProb_m)
        sigma_m = np.sqrt(hist_m)
        middle_m = np.array([(bins_m[i+1] +bins_m[i])/2.0 for i in range(len(hist_m))])
        surface_area_m = np.array([bins_m[i+1]**2 - bins_m[i]**2 for i in range(len(hist_m))])*pi
        gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
        gcows_rows = gcows[:,0]
        gcows_cols = gcows[:,1]
        plate_scale = 0.00995
        gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
        surface_area = np.zeros(len(hist))
        for i in range(len(hist)):
            gcowsdex = np.where((gcows_R > bins[i]) & (gcows_R <= bins[i+1]))[0]
            surface_area[i] = plate_scale**2 * len(gcowsdex)

        hist_m /= surface_area_m
        sigma_m /= surface_area_m

        hist,bins = np.histogram(r2d,bins=10,weights=oldProb)
        sigma = np.sqrt(hist)
        middle = np.array([(bins[i+1] +bins[i])/2.0 for i in range(len(hist))])

    else:
        schodelInfo = np.loadtxt('/u/schappell/code/c/onlySchodel_mn.dat')
        r2d = schodelInfo[:,0] / dist / cm_in_au
        oldProb = schodelInfo[:,1]
        hist,bins = np.histogram(r2d,bins=numBins,weights=oldProb)
        sigma = np.sqrt(hist)
        middle = np.array([(bins[i+1] +bins[i])/2.0 for i in range(len(hist))])
        surface_area = np.array([bins[i+1]**2 - bins[i]**2 for i in range(len(hist))])*pi

    hist /= surface_area
    sigma /= surface_area

    if (onlySchodel==False):
        middle = np.append(middle,middle_m)
        hist = np.append(hist,hist_m)
        sigma = np.append(sigma,sigma_m)

    sq_diff = np.zeros((isval[2],osval[2],sval[2],bval[2]))
    inner_slope = np.array([isval[0]+((isval[1]-isval[0])*i/(isval[2]-1)) for i in range(isval[2])])
    outer_slope = np.array([osval[0]+((osval[1]-osval[0])*i/(osval[2]-1)) for i in range(osval[2])])
    sharpness = np.array([sval[0]+((sval[1]-sval[0])*i/(sval[2]-1)) for i in range(sval[2])])
    break_R = np.array([bval[0]+((bval[1]-bval[0])*i/(bval[2]-1)) for i in range(bval[2])])

    for isi in range(isval[2]):
        is_tmp = inner_slope[isi]
        for osi in range(osval[2]):
            os_tmp = outer_slope[osi]
            for si in range(sval[2]):
                s_tmp = sharpness[si]
                for bi in range(bval[2]):
                    b_tmp = break_R[bi]
                    model = (middle/b_tmp)**(-1.0*is_tmp) * (1.0 + (middle/b_tmp)**s_tmp)**((is_tmp-os_tmp)/s_tmp)

                    model *= np.sum(hist) / np.sum(model)
                    #if ((np.sum(model))!=(np.sum(model))):
                     #   pdb.set_trace()
                    sq_diff[isi,osi,si,bi] = np.sum((hist - model)**2/sigma**2)

    py.clf()
    pdb.set_trace()






def maserStars_out(maserDate='13jullgs1',innerCut=5.0,outerCut=15.0,magCut=15.5,schodel=False,lmagCut=0.0,
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

    return R2d[mdex],oldProb[mdex],mag[mdex]




def GCOWS_schodel_overlap(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
                          poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',Rcut=5.0,
                          nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0):

    ak_correct = True
    nonRadial = 1

    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
                                         Rcut=Rcut,points=points,polyj=polyj,f_test=True,pvalue=pvalue,ak_correct=ak_correct,
                                         chainsDir=chainsDir,starlist='all',magCut=magCut,nEpochs=nEpochs,nonRadial=nonRadial)


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
 

    #Save info to file for c to read
    r2d = r2d[pdex] / dist / cm_in_au
#    ar = ar[pdex]
 #   are = are[pdex]
    oldProb = oldProb[pdex]



    #now schodel stars
    R2d_schodel = np.array([])
    mag_schodel = np.array([])
    x_schodel = np.array([])
    y_schodel = np.array([])

    ext_scale = 0.02705 #arcsec/pixel
    ext_center = [808,611]
    extMap = pyfits.getdata('/u/schappell/Downloads/AKs_fg6.fits')
    
    dbfile = '/g/ghez/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

    cur.execute('SELECT r2d,k,x,y FROM schoedel2010')
    for row in cur:
        Xext = int(round(row[2]/ext_scale))+ext_center[0]
        Yext = int(round(row[3]/ext_scale))+ext_center[0]
        if ((Xext < 0) | (Yext < 0)):
            print 'Something is wrong, star is calculated as being off extinction map'
            return
        else:
            if ((Xext < 1600) & (Yext < 1600)):
                R2d_schodel = np.append(R2d_schodel,row[0])
                mag_schodel = np.append(mag_schodel,row[1] + 2.7 - extMap[Xext,Yext])
                x_schodel = np.append(x_schodel,row[2])
                y_schodel = np.append(y_schodel,row[3])

    in_gcows = np.zeros(len(mag_schodel))
    plate_scale = 0.00995
    gcows_zero = np.array([1500.0,1500.0])
    gcows = pyfits.getdata('/u/schappell/Downloads/NIRC2 radial mask/nirc2_gcows_2010_all_mask.fits')
    gcows_dex = np.where(gcows > 0)
        #saves x and y positions, in pixels
    for ii in range(len(mag_schodel)):
            #pdb.set_trace()
        try:
            x_pixels = int(round((-1.0*x_schodel[ii] / plate_scale) + gcows_zero[0]))
            y_pixels = int(round((y_schodel[ii] / plate_scale) + gcows_zero[1]))
            in_gcows[ii] = gcows[y_pixels,x_pixels]
        except:
            continue
        if ((x_pixels < 0.0) | (y_pixels < 0.0)):
                in_gcows[ii] = 0
    idx = np.where((in_gcows==1) & (R2d_schodel < Rcut))[0]

    R2d_schodel = R2d_schodel[idx]
    mag_schodel = mag_schodel[idx]
    x_schodel = x_schodel[idx]
    y_schodel = y_schodel[idx]

    gcows_dex = np.where(gcows > 0)

    schodelHist,bins = np.histogram(R2d_schodel)
    gcowsHist,gcows_bins = np.histogram(r2d,bins=bins)

    gcows = np.loadtxt('/u/schappell/code/c/gcows_field.dat')
    gcows_rows = gcows_dex[0]
    gcows_cols = gcows_dex[1]
    gcows_R = np.sqrt((gcows_rows-1500.0)**2 + (gcows_cols-1500.0)**2)*plate_scale
    gcows_surface = np.zeros(len(gcowsHist))
    for i in range(len(schodelHist)):
        gcowsdex = np.where((gcows_R > bins[i]) & (gcows_R <= bins[i+1]))[0]
        gcows_surface[i] = plate_scale**2 * len(gcowsdex)

    schodel_middle = np.array([(bins[i+1]+bins[i])/2.0 for i in range(len(bins)-1)])
    dex = range(len(gcowsHist))

    py.clf()
    py.errorbar(schodel_middle[dex],gcowsHist/gcows_surface[dex],yerr=np.sqrt(gcowsHist)/gcows_surface[dex],fmt='ok')
    py.errorbar(schodel_middle,schodelHist/gcows_surface,yerr=np.sqrt(schodelHist)/gcows_surface,fmt='ob')
    py.xlabel('Radius (arcsec)')
    py.ylabel('Surface Density (arcsec$^{-2}$)')
    py.savefig('/u/schappell/plots/overlap_surface_density.png')

    py.clf()
    bins = np.array([6.0+0.5*i for i in range(29)])
    py.hist(mag_schodel,bins=bins,label='Schodel',histtype='step',linewidth=2.0)
    py.hist(mag,bins=bins,label='GCOWS',histtype='step',linewidth=2.0)
    py.legend(loc=2)
    py.xlabel("K' (mag)")
    py.ylabel('Number stars')
    py.savefig('/u/schappell/plots/overlap_KLF.png')





def IAU16_proc(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
               poly='polyfit_nz/fit',points='points_nz/',starlist='all',chainsDir='efit/chains/',
               Rcut=1.7,nEpochs=14.,magCut=15.5,globalt0=2006.0,polyj='polyfit_nzj/fit',pvalue=4.0,
               file1='data_oldC_oldDen_free_maser_cuts_5.0_15.0arcsec',label1='W/ Accel',
               file2='data_no_accel_gcows_free_maser_cuts_5.0_15.0arcsec',label2='W/o Accel',numBins=50):

    #Make cuts in the sample of stars used, in mag, #Epochs, and 2D radius
    #all in cgs units
    names, r2d, ar, are, mag = accelInfo(alnDir=alnDir,root_tmp=root_tmp,updateErr=updateErr,align=align,poly=poly,
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

    connection = sqlite.connect(dbfile)
    # Create a cursor object
    cur = connection.cursor()

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

    #pdb.set_trace()
    for pp in range(len(r2d)):
        ar_cmss = -ar[pp] * 1.0
        are_cmss = are[pp] * 1.0
        r_cm = r2d[pp] * 1.0
        #pdb.set_trace()
        if (ar_cmss + (3.0*are_cmss) <= (GM/r_cm**2)):
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print (str(names[pp])+' has significant ar, '+str((ar_cmss+(3.0*are_cmss))*sec_in_yr/asy_to_kms/1e5)+
                   ' and lower a_z, '+str(az_cmss*sec_in_yr/asy_to_kms/1e5))


    #Index of stars used for multinest
    pdex = np.where(oldProb > 0.0)[0]
    #pdex = np.where(oldProb == 1.0)[0]
    num_stars = len(pdex)
    #tmprint = 0
    #print 'Stars used:'
    #for p in pdex:
     #   print tmprint+5,names[p],oldProb[p]
      #  tmprint += 1


    #Save info to file for c to read
    savenames = [names[pp] for pp in pdex]
    r2d = r2d[pdex] / cm_in_pc
    ar = ar[pdex]*sec_in_yr/1.0e5
    are = are[pdex]*sec_in_yr/1.0e5
    oldProb = oldProb[pdex]

    r2dline = np.arange(0.01*cm_in_pc,0.07*cm_in_pc,1e13)
    a0line = (GM/r2dline**2 * sec_in_yr / 1.0e5)

    reported = ['S0-26','S0-7','S0-4','S0-40','S0-1','S0-3','S0-19','S0-2','S0-16','S0-8','S0-20',
                'S0-103','S0-28','S0-27','S0-52','S0-61','S0-30','S0-38','S0-17', 'S0-68','S0-70',
                'S1-2','S1-3','S0-36','S1-12','S1-13','S1-3','S0-15','irs16C','irs16SW','S1-12',
                'S1-14','S1-3','S1-12']

    py.clf()
    py.subplot(211)
    for pp in range(len(r2d)):
        if (-ar[pp] - (1.0*are[pp]) > 0.0):
            r_cm = r2d[pp] * dist * cm_in_au
            ar_cmss = ar[pp] * asy_to_kms * 1e5 / sec_in_yr
            are_cmss = are[pp] * asy_to_kms * 1e5 / sec_in_yr
            z0_cm = np.sqrt(abs(((GM * r_cm / (ar_cmss+3.0*are_cmss))**2)**(1.0/3.0) - r_cm**2)) #abs value (no sign)
            az_cmss = GM * z0_cm / ((r_cm**2 + z0_cm**2)**(1.5)) #gives lower limit of a_z
            print str(savenames[pp])+' has significant ar, '+str(-ar[pp]-(3.0*are[pp]))
            p1=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='b*',ms=14,capsize=4,elinewidth=2)
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    p2=py.errorbar(r2d[pp],-ar[pp],yerr=are[pp],fmt='gd',ms=14,capsize=4,elinewidth=2)

        else:
            print str(savenames[pp])+' has upperlimit, '+str(pp)
            upperlimit = -ar[pp]+3.0*are[pp]
            if (upperlimit <= 0.0):
                upperlimit = 0.0
            py.plot(r2d[pp],upperlimit,'k_',ms=10)
            py.arrow(r2d[pp],upperlimit,0,-1,hold=True,color='black',width=0.0003,
                     head_length=1,linewidth=0.0005,label='Upper Limits')
            for ii in range(len(reported)):
                if (savenames[pp] == reported[ii]):
                    print 'have problem'
                    pdb.set_trace()

    py.plot(r2dline/cm_in_pc,a0line,'--k')
    py.legend([p1[0],p2[0]],['New Sources','Prev. Reported'],numpoints=1,loc=2)
    py.text(0.017, 14, r'$|$a$|_{max} = \frac{GM}{R^2}$',fontsize=20)
    #py.yscale('log')
    py.xlim(0.01,0.07)
    py.ylim(-2,20)
    py.xlabel('Projected Radius (pc)')
    py.ylabel(r'$|$a$_R|$ (km/s/yr)')
    #py.title('Acceleration Detections & Upper Limits')



    posterior1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_post_equal_weights.dat')
    priors1 = np.loadtxt('/u/schappell/pmnOld/'+file1+'_priors.txt')
    gamma1 = posterior1[:,0]*(priors1[0,1] - priors1[0,0]) + priors1[0,0]
    alpha1 = posterior1[:,1]*(priors1[1,1] - priors1[1,0]) + priors1[1,0]
    delta1 = posterior1[:,2]*(priors1[2,1] - priors1[2,0]) + priors1[2,0]
    rbreak1 = posterior1[:,3]*(priors1[3,1] - priors1[3,0]) + priors1[3,0]
    cfact1 = posterior1[:,4]
    weights1 = posterior1[:,5]

    posterior2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_post_equal_weights.dat')
    priors2 = np.loadtxt('/u/schappell/pmnOld/'+file2+'_priors.txt')
    gamma2 = posterior2[:,0]*(priors2[0,1] - priors2[0,0]) + priors2[0,0]
    alpha2 = posterior2[:,1]*(priors2[1,1] - priors2[1,0]) + priors2[1,0]
    delta2 = posterior2[:,2]*(priors2[2,1] - priors2[2,0]) + priors2[2,0]
    rbreak2 = posterior2[:,3]*(priors2[3,1] - priors2[3,0]) + priors2[3,0]
    cfact2 = posterior2[:,4]
    weights2 = posterior2[:,5]

    py.subplot(212)
    py.hist(gamma1,bins=numBins,normed=1,weights=weights1,histtype='step',linewidth=2,label=label1,color='k')
    py.hist(gamma2,bins=numBins,normed=1,weights=weights2,histtype='step',linewidth=2,label=label2,ls='dashed',color='k')
    py.legend()
    py.xlabel(r'$\gamma$ (Inner Slope)')
    py.ylabel('Normalized Posterior')
    py.savefig('/u/schappell/plots/compareTwoPosts_gamma_AND_ar_vs_r2d.png')






def mass_in_radius(posteriorFile,priorsFile,max_r=1.0,mass_max=1e6,m_radius=0.01,percLevels=[0.6827, .95, .997]):

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]


    rho = np.zeros(len(gamma_array))

    for i in range(len(gamma_array)):
        massNorm = integrate.quad(lambda x: x**2 * x**(-1.0*gamma_array[i])*
                                  (1.0+(x/rbreak_array[i])**delta_array[i])**((gamma_array[i]-alpha_array[i])/delta_array[i]),
                                  0.0, max_r)
        mass_atr = integrate.quad(lambda x: x**2 * x**(-1.0*gamma_array[i])*
                                  (1.0+(x/rbreak_array[i])**delta_array[i])**((gamma_array[i]-alpha_array[i])/delta_array[i]),
                                  0.0, m_radius)
        rho[i] = mass_atr[0] * mass_max / massNorm[0]

    
    np.savetxt('/u/schappell/mass_'+str(m_radius)+'pc_chains.dat',np.transpose([rho]),delimiter=' ')
    tmpHist,tmpbins = np.histogram(rho,bins=1000)
    levels = getContourLevels(tmpHist,percLevels=percLevels)
    for i in range(len(levels)):
        tmpdex = np.where(tmpHist >= levels[i])[0]
        print 'Confidence level: '+str(percLevels[i])
        print 'Mass range: '+str(tmpbins[min(tmpdex)])+' to '+str(tmpbins[max(tmpdex)+1])+' solar masses'





def mass_loss(posteriorFile,priorsFile,max_r=1.0,mass_max=1e6,gamma2=1.75,percLevels=[0.6827, .95, .997]):

    mass2 = 4.0 * pi * max_r**(-1.0*gamma2 + 3.0) / (-1.0*gamma2 + 3.0)
    rho2_prime = max_r**(-1.0*gamma2)

    posterior = np.loadtxt(posteriorFile)
    priors = np.loadtxt(priorsFile)
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]

    mass_diff = np.zeros(len(gamma_array))

    for i in range(len(gamma_array)):
        mass1 = integrate.quad(lambda x: x**2 * x**(-1.0*gamma_array[i])*
                               (1.0+(x/rbreak_array[i])**delta_array[i])**((gamma_array[i]-alpha_array[i])/delta_array[i]),
                               0.0, max_r)
        mass1 = 4.0*pi*mass1[0]
        rho1 = mass_max / mass1
        rho1_prime = max_r**(-1.0*gamma_array[i])*(1.0+(max_r/rbreak_array[i])**delta_array[i])**((gamma_array[i]-alpha_array[i])/delta_array[i])
        rho2 = rho1 * rho1_prime / rho2_prime
        mass_diff[i] = mass2*rho2 - mass_max

    
    np.savetxt('/u/schappell/massDiff_'+str(max_r)+'pc_gamma'+str(gamma2)+'_chains.dat',np.transpose([mass_diff]),delimiter=' ')
    tmpHist,tmpbins = np.histogram(mass_diff,bins=1000)
    levels = getContourLevels(tmpHist,percLevels=percLevels)
    for i in range(len(levels)):
        tmpdex = np.where(tmpHist >= levels[i])[0]
        print 'Confidence level: '+str(percLevels[i])
        print 'Mass diff range: '+str(tmpbins[min(tmpdex)])+' to '+str(tmpbins[max(tmpdex)+1])+' solar masses'





def envelope_3d(mnlabel='data_mthread_FINAL_free_maser_cuts_5.0_15.0arcsec',min_r=0.001,max_r=2.0,numr=100,norm=1e6,norm_r=1.0):

    posterior = np.loadtxt('/u/schappell/pmnOld/'+mnlabel+'_post_equal_weights.dat')
    priors = np.loadtxt('/u/schappell/pmnOld/'+mnlabel+'_priors.txt')
    gamma_array = posterior[:,0]*(priors[0,1] - priors[0,0]) + priors[0,0]
    alpha_array = posterior[:,1]*(priors[1,1] - priors[1,0]) + priors[1,0]
    delta_array = posterior[:,2]*(priors[2,1] - priors[2,0]) + priors[2,0]
    rbreak_array = posterior[:,3]*(priors[3,1] - priors[3,0]) + priors[3,0]
    cfact_array = posterior[:,4]

    dr = (max_r - min_r)/(numr - 1)
    r3d = np.array([min_r + i*dr for i in range(numr)])
    normdex = int((norm_r - min_r) / dr)
    median_rho = np.zeros(len(r3d))
    lower_rho = np.zeros(len(r3d))
    upper_rho = np.zeros(len(r3d))
    rho = np.zeros((len(gamma_array),len(r3d)))

    for i in range(len(r3d)):
        for j in range(len(gamma_array)):
            rho[j,i] = r3d[i]**(-1.0*gamma_array[j])*(1.0+(r3d[i]/rbreak_array[j])**delta_array[j])**((gamma_array[j]-alpha_array[j])/delta_array[j])
    for j in range(len(gamma_array)):
        #pdb.set_trace()
        rho[j,:] *= norm / np.sum(r3d[0:normdex]**2 * rho[j,0:normdex] * dr) / 4.0 / pi
        #rho[j,:] /= np.sum(rho[j,:])
    for i in range(len(r3d)):
        tmp_median = np.median(rho[:,i])
        median_rho[i] = np.median(rho[:,i])

        #tmpHist,tmpbins = np.histogram(tmp_rho,bins=50)
        #levels = getContourLevels(tmpHist)
        #tmpdex = np.where(tmpHist >= levels[0])[0]
        #upper_rho[i] = tmpbins[max(tmpdex)+1]
        #lower_rho[i] = tmpbins[min(tmpdex)]

        #pdb.set_trace()

        tmp_percent = 0.0
        delta_rho = tmp_median * 0.001
        upper_sigma = tmp_median * 1.0
        #print 'Found median, now find lower and upper sigma'
        while ((tmp_percent < 0.341) & (upper_sigma < np.max(rho[:,i]))):
            upper_sigma += delta_rho
            tmpdex = np.where((rho[:,i] >= tmp_median) & (rho[:,i] <= upper_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        tmp_percent = 0.0
        lower_sigma = tmp_median * 1.0
        while ((tmp_percent < 0.341) & (lower_sigma > np.min(rho[:,i]))):
            lower_sigma -= delta_rho
            tmpdex = np.where((rho[:,i] <= tmp_median) & (rho[:,i] >= lower_sigma))[0]
            tmp_percent = float(len(tmpdex))/float(len(rho[:,i]))

        #pdb.set_trace()
        #lower_sigma = min(tmp_rho)
        #upper_sigma = max(tmp_rho)

        lower_rho[i] = lower_sigma * 1.0
        upper_rho[i] = upper_sigma * 1.0

    #norm = norm / np.sum(median_rho*dR)
    #median_rho *= norm
    #lower_rho *= norm
    #upper_rho *= norm
    bw_cusp_low = r3d**(-1.5)
    bw_cusp_high = r3d**(-1.75)
    bw_cusp_low *= median_rho[normdex] / bw_cusp_low[normdex]
    bw_cusp_high *= median_rho[normdex] / bw_cusp_high[normdex]

    gamma0 = 3.0 * norm / norm_r**3 / 4.0 / pi

    #pdb.set_trace()

    py.clf()
    py.fill_between(r3d,bw_cusp_low,bw_cusp_high,color='g')
    py.plot([0.1,max_r],[gamma0,gamma0],color='r')
    py.fill_between(r3d,lower_rho,upper_rho,color='grey')
    py.plot(r3d,median_rho,color='k')
    py.xscale('log')
    py.yscale('log')
    py.xlim([0.1,max_r])
    py.ylim([1e3,1e7])
    py.xlabel('r (pc)')
    py.ylabel(r'Density (solar mass/pc$^3$)')
    py.savefig('/u/schappell/plots/density_3d_cusp.png')




def realCompareMock(real='_mthread_FINAL',mock='tmp0.5_2.5_3_0.5_C0.3_MT',label='',binNum=40):

    realGC = np.loadtxt('/u/schappell/code/c/stars_mn'+real+'.dat')
    realSC = np.loadtxt('/u/schappell/code/c/maser_mn'+real+'.dat')

    mockGC = np.loadtxt('/u/schappell/code/c/stars_mn'+mock+'.dat')
    mockSC = np.loadtxt('/u/schappell/code/c/maser_mn'+mock+'.dat')
    mockr = np.loadtxt('/u/schappell/code/c/maser_mn'+mock+'.dat') / cm_in_pc

    real_R2d = np.append(realGC[:,0],realSC[:,0]) / (cm_in_au * dist)
    real_pOld = np.append(realGC[:,3],realSC[:,1])
    real_ar = realGC[:,1]
    real_are = realGC[:,2]

    mock_R2d = np.append(mockGC[:,0],mockSC[:,0]) / (cm_in_au * dist)
    mock_pOld = np.append(mockGC[:,3],mockSC[:,1])
    mock_ar = mockGC[:,1]
    mock_are = mockGC[:,2]

    py.clf()
    py.hist(mockr,bins=100,color='g',normed=1)
    py.xlabel('3D r (pc)')
    py.ylabel('Normalied Histogram')
    py.savefig('/u/schappell/plots/mock_r_hist_'+label+'.png')

    py.clf()
    histm,bins,junk=py.hist(mock_R2d*mock_pOld,bins=binNum,histtype='step',label='Mock Data',color='g',normed=1)
    histr,bins,junk=py.hist(real_R2d*real_pOld,bins=bins,histtype='step',label='Real Data',color='b',normed=1)
    py.legend(loc=2)
    py.xlabel('Projected R (arcsec)')
    py.ylabel('Normalized Histogram')
    py.savefig('/u/schappell/plots/compareR2d_'+label+'.png')

    py.clf()
    histm,bins,junk=py.hist(mock_ar,bins=binNum,histtype='step',label='Mock Data',color='g',normed=1)
    histr,bins,junk=py.hist(real_ar,bins=bins,histtype='step',label='Real Data',color='b',normed=1)
    py.legend(loc=2)
    py.xlabel(r'$a_R$ (cm/s$^2$)')
    py.ylabel('Normalized Histogram')
    py.savefig('/u/schappell/plots/compare_aR_'+label+'.png')

    py.clf()
    histm,bins,junk=py.hist(mock_ar/mock_are,bins=binNum,histtype='step',label='Mock Data',color='g',normed=1)
    histr,bins,junk=py.hist(real_ar/real_are,bins=bins,histtype='step',label='Real Data',color='b',normed=1)
    py.legend(loc=2)
    py.xlabel(r'$a_R$ (sigma)')
    py.ylabel('Normalized Histogram')
    py.savefig('/u/schappell/plots/compare_aR_sigma_'+label+'.png')


