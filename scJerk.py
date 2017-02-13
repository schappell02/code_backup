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
import scipy.optimize as opter
from matplotlib.ticker import ScalarFormatter 
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
    origin_val = asciidata.open('/g/ghez/align/13_08_21/' + chainsDir + 'efit_summary.txt')
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
              align = 'align/align_d_rms_1000_abs_t', updateErr = True,
              poly='polyfit_nzj/fit', points='points_nz/', sigma=3.0,
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains_S0-2_newRV2/',
              magCut=22, nEpochs=24, verbose = True):
    #Operates the same as the histAccel function in scAccel.py, but using
    #jerk fit instead of accel fit to find significant accel and jerk

    #alnDir, root_tmp, align, poly, points, and chainsDir the date directory,
    #the root directory of the align, the align directory, and polyfit directory
    #the points directory, and the directory used with the latest fit of S0-2's
    #orbit to give the location anf velocity of Sgr A* in order to update
    #the errors
    #f_test runs the f test for vel/acc and acc/jerk
    #updateErr is whether you update the errors or not
    #p value is that for the f test
    #magCut and nEpochs is the cut in magnitude and epochs


    #First, let's look at the jerk fits and get significant accel and jerks
    s = starset.StarSet(root_tmp + alnDir + align)
    allnames = s.getArray('name')
    cnt = s.getArray('velCnt')
    mag = s.getArray('mag')
    
    fitRoot = root_tmp + alnDir + poly
    fitFile = fitRoot + '.accelFormal'
    t0File = fitRoot + '.t0'

    _fit = asciidata.open(fitFile)
    _t0 = asciidata.open(t0File)

    t0x = _t0[1].tonumpy()
    t0y = _t0[2].tonumpy()
    names = _fit[0].tonumpy()
#    pdb.set_trace()
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
    
    #Not all stars have jerk fits though
    
    r = np.hypot(x0,y0)
    
