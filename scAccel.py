#from gcwork 
import objects
#from gcwork 
import starset
import sc_accel_class as acc
#from gcwork 
import util
#from gcwork 
import orbits
#from gcwork 
import young
from pysqlite2 import dbapi2 as sqlite
import scipy
import pyfits
from scipy import stats
from scipy import special
from scipy import integrate
#from gcwork 
import starTables
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
from scipy.optimize import fsolve
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

def plotVelocityMap(alnDir='13_08_21/',
                    align='align/align_d_rms_1000_abs_t', poly='polyfit_c/fit'):

    outdir = home + 'plots/'
    s = starset.StarSet(root + alnDir + align)
    s.loadPolyfit(root + alnDir + poly)

    name = s.getArray('name')
    x = s.getArray('fitXv.p')
    y = s.getArray('fitYv.p')
    vx = s.getArray('fitXv.v')
    vy = s.getArray('fitYv.v')
    vxe = s.getArray('fitXv.verr')
    vye = s.getArray('fitYv.verr')
    vx *= 1e3
    vy *= 1e3
    vxe *= 1e3
    vye *= 1e3

    # temp
    cnt = s.getArray('velCnt')
    idx = np.where(cnt > 20)[0]
    print 'Ratio of X to Y proper motion errors: %4.3f' % (vxe[idx] / vye[idx]).mean()

    #1 outlier:
    #foo = np.where(vy > 500)[0]
    #print name[foo],x[foo],y[foo],vx[foo],vy[foo]

    py.clf()
    arrScale = 1.0
    qvr = py.quiver([x], [y], [vx], [vy], headwidth=3, color='black')
    py.quiverkey(qvr, 4.0, 5.0, 15, '15 mas/yr', coordinates='data')
    py.savefig(outdir + 'velMap.png')
    #py.savefig(outdir + 'velMap.eps')

    # To get young stars:
    #yng = youngNames.loadYoungStars()
    #print [yng[ii] for ii in range(len(yng))]

    py.clf()
    py.figure(figsize=(12,5))
    py.subplots_adjust(wspace=0.25, hspace=0.25, left=0.08, right=0.95, top=0.95)
    py.subplot(1,3,1)
    py.semilogy(np.hypot(x,y), np.abs(vx), 'r.')
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('|X Velocity| (mas/yr)')
    py.subplot(1,3,2)
    py.semilogy(np.hypot(x,y), np.abs(vy), 'b.')
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('|Y Velocity| (mas/yr)')
    py.subplot(1,3,3)
    py.semilogy(np.hypot(x,y), np.hypot(vx,vy), 'k.')
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('Proper Motion (mas/yr)')
    py.savefig(outdir + 'pm_vs_r2d.png')

    

def loadPop(alnDir='14_06_18/',
            align = 'align/align_d_rms_1000_abs_t', root_tmp='/g/ghez/align/',
            poly='polyfit_nz/fit',points='points_nz/', starlist='all',to_debug=False):
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


def nzErr(xerr, yerr, vxerr, vyerr, year_x, year_y, mag, root_tmp = '/g/ghez/align/', alnDir = '13_08_21/', 
          chainsDir = 'efit/chains/'):
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
    origin_val = asciidata.open(root_tmp + alnDir + chainsDir + 'efit_summary.txt')
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
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',outTolatex=False,outToAL=False,
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


    #pdb.set_trace()

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
#    py.clf()
#    pdb.set_trace()
#    r = r[idx]
#    ag1 = np.where((r < 1.7) & (cnt > 14))[0]
#    ag2 = np.where((r < 1.7) & (cnt > 14) & (mag < 15.5))[0]
#    co1 = plt.scatter(cnt,mag,facecolors='none',label='Entire Catalogue')
#    co2 = py.scatter(cnt[ag1],mag[ag1],label='Projected Radius and Epochs Cut')
#    py.scatter(mag[ag2],cnt[ag2],'ok')
#    py.xlabel('Number Epochs')
#    py.ylabel('K Mag')
#    py.axis([0,60,22,8])
#    py.legend((co1,co2),['Entire Catalogue','Proj. R + Epochs Cut'])
#    py.show()
#    pdb.set_trace()

#    co1 = py.scatter(cnt[ag1],mag[ag1],facecolors='none',label='Projected Radius and Epcohs Cut')
#    co2 = py.scatter(mag[ag2],cnt[ag2],'ok',label='+K Magnitude Cut')
#    py.xlabel('Number Epochs')
#    py.ylabel('K Mag')
#    py.show()

    #pdb.set_trace()


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


    for i in range(len(names)):
        if ((str(names[i]) == 'S0-16') | (str(names[i]) == 'S0-17')):
            pass_f[i] = 1
            pass_fj[i] = 1
            print 'S0-16 or S0-17 updated accel and jerk'
            #pdb.set_trace()
    idex = np.where(pass_fj ==1)[0]

    pdb.set_trace()

    # S0-16 and -17 have such a high inclination that it is being seen as a velocity, force the f test vel/acc to be acc

    idex = np.where((pass_fj ==1) & (pass_f == 1))[0]

    for i in idex:
        jdex = np.where(jnames == names[i])
        try:
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
        except:
            pass_fj[i] = 0 #star does not have jerk fit, set pass jerk F test to no 
            print names[i]+' did not have jerk fit'
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
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, root_tmp=root_tmp, alnDir=alnDir, chainsDir = chainsDir)
#        jx0e,jy0e,jvxe,jvye = nzErr(jx0e, jy0e, jvxe, jvye, jt0x, jt0y, jmag, alnDir=alnDir, chainsDir = chainsDir)

    if updateErr:
        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17.1])
        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5

        delta = mag*0.0


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

    #pdb.set_trace()

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
#    data.computeFTestJerk()

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

def plotVvsR(name,rmin,vx,vy,vxe,vye):

    # Construct an array of radii out to 4 arcsec in steps of 0.1''
    r = (arange(10 * 5) * 0.1) + 0.1
    r_au = r * dist
    r_pc = r_au / au_in_pc
    r_cm = r_au * cm_in_au
    r_km = r_cm / 1.e5

    # Compute total proper motion and error for this star
    v_kms = np.sqrt(vx**2+vy**2)
    ve_kms = np.sqrt(vx**2+vy**2)

    # NOTE: still need to account for the error in the 3D radius
    #       and error in the escape velocity curve

    # Find fast objects
    v_esc_cms = np.sqrt((2. * G * mass * msun) / r_cm) # for curve on plot
    v_esc_kms = v_esc_cms / 1.e5

     # Plot velocities vs. r, w/ error bars
    plot([rmin],[v_kms],'b.')
    errorbar([rmin],[v_kms],yerr=ve_kms,fmt='b.')
    plot(r_pc, v_esc_kms,'r-')
    if name == 'S1-13':
        vtot_kms = np.sqrt(vx**2+vy**2+745.**2)
        plot([rmin],[vtot_kms],'r.')
        text(rmin + 0.002,vtot_kms,name,fontsize=9)


def plotLimits(alnDir='13_08_21/',
               align='align/align_d_rms_1000_abs_t',
               poly='polyfit_c/fit', points='points_c/',
               youngOnly=True, sigma=3):
    """
    Find the 3 sigma radial acceleration limits for each star.
    """
    # Load GC constants
    #cc = objects.Constants()

    # Load up positional information from align. 
    s = starset.StarSet(root + alnDir + align)
    s.loadPolyfit(root + alnDir + poly, arcsec=1)
    s.loadPolyfit(root + alnDir + poly, arcsec=1, accel=1)

    names = s.getArray('name')

    # In arcsec
    x = s.getArray('x')
    y = s.getArray('y')
    # In arcsec/yr
    vx = s.getArray('fitXa.v')
    vy = s.getArray('fitYa.v')
    # In mas/yr^2
    ax = s.getArray('fitXa.a') * 1000.0
    ay = s.getArray('fitYa.a') * 1000.0
    axe = s.getArray('fitXa.aerr') * 1000.0
    aye = s.getArray('fitYa.aerr') * 1000.0
    r2d = s.getArray('r2d')
    cnt = s.getArray('velCnt')

    tag = '_all'
    tag2 = 'All'

    if (youngOnly == True):
        yngNames = young.youngStarNames()

        idx = []

        for ii in range(len(names)):
            #if (r2d[ii] > 0.8 and
            if (names[ii] in yngNames and cnt[ii] >= 20):
                idx.append(ii)

        names = [names[i] for i in idx]
        x = x[idx]
        y = y[idx]
        vx = vx[idx]
        vy = vy[idx]
        ax = ax[idx]
        ay = ay[idx]
        axe = axe[idx]
        aye = aye[idx]
        r2d = r2d[idx]
        print 'Found %d young stars' % len(names)
        tag = '_yng'
        tag2 = 'Young'

    # Lets do radial/tangential
    r = np.sqrt(x**2 + y**2) # arcsec
    if ('Radial' in poly):
        at = ax
        ar = ay
        ate = axe
        are = aye
    else:
        # In mas/yr^2
        ar = ((ax*x) + (ay*y)) / r 
        at = ((ax*y) - (ay*x)) / r
        are = np.sqrt((axe*x)**2 + (aye*y)**2) / r
        ate = np.sqrt((axe*y)**2 + (aye*x)**2) / r

    # Total acceleration
    atot = py.hypot(ax, ay)
    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot
    # Acceleration upper limit for each star
    accLim = ar - (sigma * are)

    # Calculate the acceleration limit set by the projected radius
    # Convert into cm
    r2d = r * dist * cm_in_au
    # acc1 in cm/s^2
    a2d = -G * mass * msun / r2d**2
    # acc1 in km/s/yr
    a2d *= sec_in_yr / 1.0e5
    # acc1 in mas/yr^2
    a2d *= 1000.0 / asy_to_kms
    accLimR2d = a2d

    _f = open(root + alnDir + 'tables/accel_limits' + tag + '.txt', 'w')

    _f.write('###   Radial Acceleration Limits for %s Stars   ###\n' % tag2)
    _f.write('%13s  %7s - (%d * %5s) =  %9s     %21s  Constrained?  Non-zero?\n' % \
             ('Name', 'a_rad', sigma, 'aerr', 'a_obs_lim', 'a_proj_lim'))
    
    fmt = '%13s  %7.3f - (%d * %5.3f) =    %7.3f  '
    fmt += 'vs %10.3f (mas/yr^2)  %8s  %8s'
    fmt2 = fmt + '\n'
    fmt3 = '%13s  %7.3f - (%d * %5.3f) =    %7.3f  vs %10.3f (mas/yr^2)  %5.3f\n'

    hiPri = []
    for i in range(len(ar)):
        constrained = ''
        ltzero = ''
        if (accLim[i] > accLimR2d[i]): # Below the curve (remember, negative values)
            constrained = '**'
            
            # Find those with significant accelerations below zero
            if ((ar[i] < 0.) & ((ar[i] + (sigma*are[i])) < 0.)):
                ltzero = '**'
                # Save these indices as high priority
                hiPri = np.concatenate([hiPri,[i]])

        # Print info for all young stars
        print fmt % (names[i], ar[i], sigma, are[i], accLim[i],
                     accLimR2d[i], constrained, ltzero)
        _f.write(fmt2 % (names[i], ar[i], sigma, are[i], accLim[i],
                         accLimR2d[i], constrained, ltzero))

    _f.close()

    # Make a high-priority list for stars w/ acc. detections
    _hipri = open(root + alnDir + 'tables/accel_limits' + tag + '_hiPriority.txt', 'w')
    _hipri.write('###   Radial Acceleration Limits for %s Stars   ###\n' % tag2)
    _hipri.write('%13s  %7s - (%d * %5s) =  %9s     %21s   %5s\n' % \
             ('Name', 'a_rad', sigma, 'aerr', 'a_obs_lim', 'a_proj_lim', 'r2d (")'))
    hiPri = [int(ii) for ii in hiPri]
    for h in hiPri:
        _hipri.write(fmt3 % (names[h], ar[h], sigma, are[h], accLim[h],
                         accLimR2d[h], r[h]))
    _hipri.close()

    # Stuff below confirms that these accelerations match those in accelLimit()
    # Convert accelerations to km/s/yr
    #ar = ar / 1000.0 * asy_to_kms 
    #are = are / 1000.0 * asy_to_kms 
    #accLim = ar - (sigma * are)
    #accLimR2d = accLimR2d / 1000.0 *asy_to_kms
    #clf()
    #plot(r2d/cm_in_pc,-accLim,'k.')
    #foo = argsort(r)
    #plot(r[foo]*0.04,-accLimR2d[foo],'--k')
    #axis('normal')
    #axis([0,0.18,0,40])
    #savefig(root + alnDir + 'plots/r2d_accel_limit_yng2.png')



def compareVelocity(alnDir='11_10_26/',
                    align = 'align/align_d_rms_1000_abs_t',
                    poly='polyfit_c/fit', points='points_c/'):
    """
    Compare velocities from the linear and acceleration fits.

    alnDir -- An align root directory such as '10_04_22/' 
    align -- Align root name 
    poly -- Polyfit root name 
    points -- Points directory
    """
    
    s = starset.StarSet(root + alnDir + align)
    s.loadPolyfit(root + alnDir + poly, accel=0, arcsec=1)
    s.loadPolyfit(root + alnDir + poly, accel=1, arcsec=1)

    names = s.getArray('name')
    x = s.getArray('x')
    y = s.getArray('y')
    r = hypot(x, y)

    vx_vel = s.getArray('fitXv.v')
    vy_vel = s.getArray('fitYv.v')
    vxe_vel = s.getArray('fitXv.verr')
    vye_vel = s.getArray('fitYv.verr')

    vx_acc = s.getArray('fitXa.v')
    vy_acc = s.getArray('fitYa.v')
    vxe_acc = s.getArray('fitXa.verr')
    vye_acc = s.getArray('fitYa.verr')
    
    # Calculate the residuals
    diffx = vx_vel - vx_acc
    diffxErr = np.sqrt(vxe_vel**2 + vxe_acc**2)
    
    diffy = vy_vel - vy_acc
    diffyErr = np.sqrt(vye_vel**2 + vye_acc**2)
    
    diff = py.hypot(diffx, diffy)
    diffErr = np.sqrt((diffx*diffxErr)**2 + (diffy*diffyErr)**2) / diff

    yngNames = young.youngStarNames()
    
    idx = (np.where((r > 0.8) & ((diff/diffErr) > 3.0)))[0]

    print '** Stars with large velocity discrepencies (N=%d): **' % len(idx)
    print '%15s  %5s  %5s  %5s  %3s' % ('Name', 'Sigma', 'X', 'Y', 'YNG?')
    for i in idx:
        if names[i] in yngNames:
            yngString = 'yng'
        else:
            yngString = ''
            
        print '%15s  %5.1f  %5.2f  %5.2f  %3s' % \
              (names[i], diff[i]/diffErr[i], x[i], y[i], yngString)

    # get the velocities for just the young stars
    vVx_yng = []
    vVy_yng = []
    aVx_yng = []
    aVy_yng = []
    vVxe_yng = []
    vVye_yng = []
    aVxe_yng = []
    aVye_yng = []
    for yy in range(len(yngNames)):
        #sidx = np.where(names == yngNames[yy])[0]
        if yngNames[yy] not in names:
            continue
        sidx = names.index(yngNames[yy])
        vVx_yng = np.concatenate([vVx_yng, [vx_vel[sidx]*1.e3]])
        vVy_yng = np.concatenate([vVy_yng, [vy_vel[sidx]*1.e3]])
        aVx_yng = np.concatenate([aVx_yng, [vx_acc[sidx]*1.e3]])
        aVy_yng = np.concatenate([aVy_yng, [vy_acc[sidx]*1.e3]])
        vVxe_yng = np.concatenate([vVxe_yng, [vxe_vel[sidx]*1.e3]])
        vVye_yng = np.concatenate([vVye_yng, [vye_vel[sidx]*1.e3]])
        aVxe_yng = np.concatenate([aVxe_yng, [vxe_acc[sidx]*1.e3]])
        aVye_yng = np.concatenate([aVye_yng, [vye_acc[sidx]*1.e3]])
        print 'Young star %6s: linear fit vx,vy = %6.2f, %6.2f mas/yr' % \
              (yngNames[yy], vx_vel[sidx]*1.e3, vy_vel[sidx]*1.e3)
        print '                   accel fit vx,vy = %6.2f, %6.2f mas/yr' % \
              (vx_acc[sidx]*1.e3, vy_acc[sidx]*1.e3)

    # Plot just young stars; X and Y seperately
    py.clf()
    py.subplot(2, 1, 1)
    py.errorbar(vVx_yng, aVx_yng, xerr=vVxe_yng, yerr=aVxe_yng, fmt='k.')
    rng = py.axis()
    py.plot(rng[0:2], rng[0:2], 'b--')
    py.ylabel('Accel Fit (mas/yr)')

    py.subplot(2, 1, 2)
    py.errorbar(vVy_yng, aVy_yng, xerr=vVye_yng, yerr=aVye_yng, fmt='k.')
    rng = py.axis()
    py.plot(rng[0:2], rng[0:2], 'b--')
    py.xlabel('Linear Fit (mas/yr)')
    py.ylabel('Accel Fit (mas/yr)')
    py.savefig(home + 'plots/vel_linear_vs_accel_fit_yngstars.png')
    
    

    # Plot all stars; X and Y seperately
    py.clf()
    py.subplot(2, 1, 1)
    py.errorbar(vx_vel*1.e3, vx_acc*1.e3, xerr=vxe_vel*1.e3,
                yerr=vxe_acc*1.e3, fmt='k.')
    rng = py.axis()
    py.plot(rng[0:2], rng[0:2], 'b--')
    py.ylabel('Accel Fit (mas/yr)')

    py.subplot(2, 1, 2)
    py.errorbar(vy_vel*1.e3, vy_acc*1.e3, xerr=vye_vel*1.e3,
                yerr=vye_acc*1.e3, fmt='k.')
    rng = py.axis()
    py.plot(rng[0:2], rng[0:2], 'b--')
    py.xlabel('Linear Fit (mas/yr)')
    py.ylabel('Accel Fit (mas/yr)')
    py.savefig(root + alnDir + 'plots/vel_linear_vs_accel_fit_allstars.png')
    
    
def highSigSrcs(radiusCut, sigmaCut, alnDir='13_08_21/',
                poly='polyfit_c/fit',
                save=False, verbose=True):
    
    """
    Make a list of all sources with significant accelerations.
    Assumes a plate scale of 9.95 mas/pixel.

    Inputs:
    radiusCut - the largest radius to include (arcsec)
    sigmaCut - the lowest a/a_err to include
    save - Save to tables/accelHighSigSrcs.txt (def=False)
    verbose - Print to screen (def=True)

    Outputs:
    srcNames - returns a list of stars with significant radial
    accelerations based on the passed in criteria.
    """
    # Load up the accelPolar file
    fitFile = root + alnDir + poly + '.accelPolar'
    scale = 0.00995  # arcsec/pixel

    tab = asciidata.open(fitFile)

    name = tab[0]._data
    radius = tab[2].tonumarray() * scale
    acc = tab[10].tonumarray() * scale
    accErr = tab[12].tonumarray() * scale
    sigma = acc / accErr

    # Make cuts in radius and sigma
    idx = ( where((radius < radiusCut) & (sigma > sigmaCut)) )[0]

    # Sort
    rdx = sigma[idx].argsort()
    idx = idx[rdx[::-1]]

    # Print out to the screen
    if (verbose == True):
        print '** Found %d significantly accelerating sources **' % len(idx)
        print ''
        print '%-15s  %8s  %10s  %8s' % \
              ('Name', 'Radius', 'Accel', 'Signif.')
        print '%-15s  %8s  %10s  %8s' % \
              ('', '(arcsec)', '(mas/yr^2)', '(sigma)')
        print '%-15s  %8s  %10s  %8s' % \
              ('---------------', '--------', '----------', '--------')

        for ii in idx:
            print '%-15s  %8.3f  %10.5f  %8.1f' % \
                  (name[ii], radius[ii], acc[ii], sigma[ii])


    # Save data to an output file
    if (save == True):
        outFile = 'tables/accelHighSigSrcs.txt'
        _out = open(outFile, 'w')

        _out.write('# Python: syAccel.highSigSrcs')
        _out.write('(%5.2f, %5.2f)\n' % \
                   (radiusCut, sigmaCut))
        _out.write('#\n')
        _out.write('%-15s  %8s  %10s  %8s\n' % \
                   ('# Name', 'Radius', 'Accel', 'Signif.'))
        _out.write('%-15s  %8s  %10s  %8s\n' % \
                   ('#', '(arcsec)', '(mas/yr^2)', '(sigma)'))
        _out.write('%-15s  %8s  %10s  %8s\n' % \
                   ('#--------------', '--------', '----------', '--------'))

        for ii in idx:
            _out.write('%-15s  %8.3f  %10.5f  %8.1f\n' % \
                       (name[ii], radius[ii], acc[ii], sigma[ii]))

        _out.close()

    # Return the list of significantly accelerating sources.
    return [name[ii] for ii in idx]
    
    
def velVsAcc():
    """
    Plot v/v_circular vs. a/a_bound.
    NOT COMPLETE?????
    """
    # Load up the accelPolar file
    fitFile = root + alnDir + poly + '.accelPolar'
    scale = 0.00995  # arcsec/pixel

    tab = asciidata.open(fitFile)

    name = tab[0]._data
    radius = tab[2].tonumarray() * scale

    velPhi = tab[5].tonumarray() * scale
    velRad = tab[6].tonumarray() * scale
    velPhiErr = tab[7].tonumarray() * scale
    velRadErr = tab[8].tonumarray() * scale
    acc = tab[10].tonumarray() * scale
    accErr = tab[12].tonumarray() * scale

    # Need to get the line-of-sight velocity from Paumard et al.
    

    vel = np.sqrt(velPhi**2 + velRad**2)
    velErr = np.sqrt((velPhi*velPhiErr)**2 + (velRad*velRadErr)**2) / vel

    # Determine the circular velocity 



def plotStar(starName,alnDir='14_06_18/',align='align/align_d_rms_1000_abs_t',
             poly='polyfit_nz/fit', points='points_nz/', radial=False, accel_fit=True):

#accel_fit = True the resulting plots of X and Y vs t and Y vs X show accel fit, if false show velocity fit
#either way outputs another figure showing residuals in X and Y from both fits

    outdir = home + 'plots/'

    s = starset.StarSet(rootDir + alnDir + align)
    
    s.loadPolyfit(rootDir + alnDir + poly, accel=0, arcsec=0)
    s.loadPolyfit(rootDir + alnDir + poly, accel=1, arcsec=0)

    names = s.getArray('name')
    
    ii = names.index(starName)
    star = s.stars[ii]

    pointsTab = asciidata.open(rootDir + alnDir + points + starName + '.points')

    time = pointsTab[0].tonumpy()
    x = pointsTab[1].tonumpy()
    y = pointsTab[2].tonumpy()
    xerr = pointsTab[3].tonumpy()
    yerr = pointsTab[4].tonumpy()

    fitxv = star.fitXv
    fityv = star.fitYv
    dtv = time - fitxv.t0
    fitLineXv = fitxv.p + (fitxv.v * dtv)
    fitSigXv = np.sqrt( fitxv.perr**2 + (dtv * fitxv.verr)**2 )

    fitLineYv = fityv.p + (fityv.v * dtv)
    fitSigYv = np.sqrt( fityv.perr**2 + (dtv * fityv.verr)**2 )

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
        

    diffXv = x - fitLineXv
    diffYv = y - fitLineYv
    #diff = np.hypot(diffX, diffY)
    diffv = (diffXv + diffYv) / 2.0
    rerrv = np.sqrt((diffXv*xerr)**2 + (diffYv*yerr)**2) / diffv
    sigXv = diffXv / xerr
    sigYv = diffYv / yerr
    sigv = diffv / rerrv

    #SC added, find accel fit and accel residuals
    fitxa = star.fitXa
    fitya = star.fitYa
    dta = time - fitxa.t0
    fitLineXa = fitxa.p + (fitxa.v * dta) + (fitxa.a * dta * dta / 2.0)
    print fitxa.t0
    print fitxa.p
    print fitxa.v
    print fitxa.a
    fitSigXa = np.sqrt(fitxa.perr**2 + (dta * fitxa.verr)**2 + (dta * dta * fitxa.aerr / 2.0)**2)
    fitLineYa = fitya.p + (fitya.v * dta) + (fitya.a * dta * dta / 2.0)
    print fitya.t0
    print fitya.p
    print fitya.v
    print fitya.a
    fitSigYa = np.sqrt(fitya.perr**2 + (dta * fitya.verr)**2 + (dta * dta * fitya.aerr / 2.0)**2)

    diffXa = x - fitLineXa
    diffYa = y - fitLineYa
    #diff = np.hypot(diffX, diffY)
    diffa = (diffXa + diffYa) / 2.0
    rerra = np.sqrt((diffXa*xerr)**2 + (diffYa*yerr)**2) / diffa
    sigXa = diffXa / xerr
    sigYa = diffYa / yerr
    siga = diffa / rerra

#    if (accel_fit == False):
    sigX=sigXv
    sigY=sigYv
    sig=sigv
    fitx=fitxv
    fity=fityv
    fitLineX=fitLineXv
    fitLineY=fitLineYv
    fitSigX=fitSigXv
    fitSigY=fitSigYv
    fit_tag='Vel Fit'

    if (accel_fit == True):
        sigX=sigXa
        sigY=sigYa
        sig=siga
        fitx=fitxa
        fity=fitya
        fitLineX=fitLineXa
        fitLineY=fitLineYa
        fitSigX=fitSigXa
        fitSigY=fitSigYa
        fit_tag='Accel Fit'

    # Determine if there are points that are more than 5 sigma off
    idxX = np.where(abs(sigX) > 5)
    idxY = np.where(abs(sigY) > 5)
    idx = np.where(abs(sig) > 5)

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

    figure(2, figsize=(7, 8))
    clf()

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
    plot(time, fitLineX, 'b-')
    plot(time, fitLineX + fitSigX, 'b--')
    plot(time, fitLineX - fitSigX, 'b--')
    errorbar(time, x, yerr=xerr, fmt='k.')
    rng = axis()
    axis(dateTicRng + [rng[2], rng[3]])
    xlabel('Date (yrs)')
    ylabel('X (arcsec)')
    title(starName+' '+fit_tag)
#    title('X chi2_red = %4.2f' % fitx.chi2red)
    #paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
    paxes.yaxis.set_major_formatter(fmtY)
    
#    paxes = subplot(3, 2, 2)
    paxes = subplot(2,1,2)
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

    savefig(outdir+'plotXY_t_' + starName + '.png')
    py.close(2)

    figure(3)
    clf()
    subplots_adjust(hspace=0.25,left=0.1, right=0.95,top=0.9, bottom=0.1, wspace=0.25)
    paxes=subplot(2,2,1)
    errorbar(time,diffXv*1e3,yerr=xerr*1e3,fmt='b.')
    plot([0,3000],[0,0],'k')
    axis(dateTicRng+[-10,10])
    xlabel('Date (yrs)')
    ylabel('X Residuals (mas)')
    title('Velocity Fit')
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

    paxes=subplot(2,2,2)
    errorbar(time,diffXa*1e3,yerr=xerr*1e3,fmt='b.')
    plot([0,3000],[0,0],'k')
    axis(dateTicRng+[-10,10])
    xlabel('Date (yrs)')
    ylabel('X Residuals (mas)')
    title('Acceleration Fit')
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

    paxes=subplot(2,2,3)
    errorbar(time,diffYv*1e3,yerr=yerr*1e3,fmt='b.')
    plot([0,3000],[0,0],'k')
    axis(dateTicRng+[-10,10])
    xlabel('Date (yrs)')
    ylabel('Y Residuals (mas)')
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

    paxes=subplot(2,2,4)
    errorbar(time,diffYa*1e3,yerr=yerr*1e3,fmt='b.')
    plot([0,3000],[0,0],'k')
    axis(dateTicRng+[-10,10])
    xlabel('Date (yrs)')
    ylabel('Y Residuals (mas)')
#    title('Acceleration Fit')
    paxes.xaxis.set_major_formatter(fmtX)
    paxes.get_xaxis().set_major_locator(dateTicLoc)
#    paxes.yaxis.set_major_formatter(fmtY)

    savefig(outdir+'plot_Residuals_'+starName+'.png')
    py.close(3)

    py.figure(4)
    py.clf()
    plot(fitLineX[idtmp], fitLineY[idtmp], 'b-')
    errorbar(x[idtmp], y[idtmp],xerr=xerr[idtmp], yerr=yerr[idtmp], fmt='k.')
    title(starName+' '+fit_tag)
#    rng = axis()
#    axis(xyplim)
    xlabel('X (arcsec)')
    ylabel('Y (arcsec)')
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
    savefig(outdir+'plotY_X_' + starName + '.png')
    py.close(4)
#    print outdir + 'plotStar_' + starName + '.png'
    #show()

    ##########
    #
    # Also plot radial/tangential
    #
    ##########
    if (radial == True):
        clf()

        dateTicLoc = MultipleLocator(3)
        #dateTicRng = [1995, 2011]
        dateTicRng = [2006, 2011]
        
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
        ylabel('R (arcsec)')
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
        ylabel('T (arcsec)')
        paxes.xaxis.set_major_formatter(fmtX)
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        paxes.yaxis.set_major_formatter(fmtY)
        
        paxes = subplot(3, 2, 3)
        plot(time, np.zeros(len(time)), 'b-')
        plot(time, fitSigR, 'b--')
        plot(time, -fitSigR, 'b--')
        errorbar(time, r - fitLineR, yerr=rerr, fmt='k.')
        #axis(dateTicRng + resTicRng)
        axis([2005,2011,-0.001,0.001])
        xlabel('Date (yrs)')
        ylabel('R Residuals (arcsec)')
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        
        paxes = subplot(3, 2, 4)
        plot(time, np.zeros(len(time)), 'b-')
        plot(time, fitSigT, 'b--')
        plot(time, -fitSigT, 'b--')
        errorbar(time, t - fitLineT, yerr=terr, fmt='k.')
        #axis(dateTicRng + resTicRng)
        axis([2005,2011,-0.001,0.001])
        xlabel('Date (yrs)')
        ylabel('T Residuals (arcsec)')
        paxes.get_xaxis().set_major_locator(dateTicLoc)
        
        bins = np.arange(-7, 7, 1)
        subplot(3, 2, 5)
        (n, b, p) = hist(sigR, bins)
        setp(p, 'facecolor', 'k')
        axis([-5, 5, 0, 20])
        xlabel('R Residuals (sigma)')
        ylabel('Number of Epochs')
        
        subplot(3, 2, 6)
        (n, b, p) = hist(sigT, bins)
        axis([-5, 5, 0, 20])
        setp(p, 'facecolor', 'k')
        xlabel('T Residuals (sigma)')
        ylabel('Number of Epochs')
        
        subplots_adjust(wspace=0.4, hspace=0.4, right=0.95, top=0.95)
        savefig(outdir+'plotStarRadial_' + starName + '.png')

    return fitLineXv,fitLineXa


def getResiduals(alnDir='13_08_21/',
                align='align/align_d_rms_1000_abs_t',
                poly='polyfit_c/fit', points='points_c/',
                additiveErr=True, trimOutliers=False, trimSigma=4,
                useAccFits=False, radCut=4.0,suffix=''):
    """Analyze the distribution of points relative to their best
    fit velocities. Optionally trim the largest outliers in each
    stars *.points file.  Optionally make radius cut with radCut flag.
    Set additiveErr=True to remove the additive error from positional
    errors in points files, and to plot this as a separate line.
    """

    outdir = home + 'plots/'

    s = starset.StarSet(root + alnDir + align)
    s.loadPolyfit(root + alnDir + poly, accel=0, arcsec=0)
    s.loadPolyfit(root + alnDir + poly, accel=1, arcsec=0)
    s.loadPoints(root + alnDir + points)
    names = s.getArray('name')
    x = s.getArray('x')
    y = s.getArray('y')

    addErr = 0
    if additiveErr == True:
        # Define the additive error:
        addErr = 0.1 # mas, AO data only

    # Make some empty arrays to hold all our results.
    sigmaX = np.arange(0, dtype=float)
    sigmaY = np.arange(0, dtype=float)
    sigma  = np.arange(0, dtype=float)
    diff_all = np.arange(0, dtype=float)
    diffX_all = np.arange(0, dtype=float)
    diffY_all = np.arange(0, dtype=float)
    xerr_all = np.arange(0, dtype=float)
    yerr_all = np.arange(0, dtype=float)
    mag_all = np.arange(0, dtype=float)
    xeom_all = np.arange(0, dtype=float)
    yeom_all = np.arange(0, dtype=float)
    xaln_all = np.arange(0, dtype=float)
    yaln_all = np.arange(0, dtype=float)
    starName_all = np.arange(0, dtype=float)

    nCnt = 0
    # Loop through all the stars and combine their residuals and their errors.
    for star in s.stars:
        starName = star.name
        cnt = star.pointsCnt        
        if cnt < 5:
            continue
        nCnt += 1

        ############
        # Points files - include err on mean, align err, and additive err
        #   -- only need these total errs when looking at residuals in units of 'sigma'
        ############

        t = np.array(star.years)
        errx_a_allEp = star.getArrayAllEpochs('xerr_a')*1e3 # mas
        erry_a_allEp = star.getArrayAllEpochs('yerr_a')*1e3
        errx_p_allEpTmp = star.getArrayAllEpochs('xerr_p')*1e3 # mas, includes additive
        erry_p_allEpTmp = star.getArrayAllEpochs('yerr_p')*1e3
        x_pts = star.getArrayAllEpochs('pnt_x')
        y_pts = star.getArrayAllEpochs('pnt_y')
        xe_pts = star.getArrayAllEpochs('pnt_xe') # include eom, align, additive err
        ye_pts = star.getArrayAllEpochs('pnt_ye')
        mag = star.getArrayAllEpochs('phot_mag')
        mag_e = star.getArrayAllEpochs('phot_mage')
        r2d = np.hypot(x_pts,y_pts)

        # Remove non-detections
        kp = np.where(x_pts > -999)[0]
        errx_a_allEp = errx_a_allEp[kp]
        erry_a_allEp = erry_a_allEp[kp]
        errx_p_allEpTmp = errx_p_allEpTmp[kp]
        erry_p_allEpTmp = erry_p_allEpTmp[kp]
        x_pts = x_pts[kp]
        y_pts = y_pts[kp]
        xe_pts = xe_pts[kp]
        ye_pts = ye_pts[kp]
        mag = mag[kp]
        mag_e = mag_e[kp]
        r2d = r2d[kp]
        t = t[kp]
        
        if ((radCut != None) and (r2d.mean() > radCut)):
            continue

        # Remove additive err
        if additiveErr == True:
            errx_p_allEp = np.sqrt(errx_p_allEpTmp**2 - addErr**2)  
            erry_p_allEp = np.sqrt(erry_p_allEpTmp**2 - addErr**2)  
        else:
            errx_p_allEp = errx_p_allEpTmp
            erry_p_allEp = erry_p_allEpTmp
            
        
        ############
        # Residuals - Best fit velocity model
        ############
        if (useAccFits == True):
            fitx = star.fitXa
            fity = star.fitYa
        else:
            fitx = star.fitXv
            fity = star.fitYv

        if np.isnan(fitx.t0) == True:
            pdb.set_trace()
        #    continue

        dt = t - fitx.t0
        fitLineX = fitx.p + (fitx.v * dt)
        fitSigX = np.sqrt( fitx.perr**2 + (dt * fitx.verr)**2 )

        fitLineY = fity.p + (fity.v * dt)
        fitSigY = np.sqrt( fity.perr**2 + (dt * fity.verr)**2 )

        if (useAccFits == True):
            fitLineX += (fitx.a * dt**2) / 2.0
            fitLineY += (fity.a * dt**2) / 2.0
            fitSigX = np.sqrt(fitSigX**2 + (dt**2 * fitx.aerr / 2.0)**2)
            fitSigY = np.sqrt(fitSigY**2 + (dt**2 * fity.aerr / 2.0)**2)

        # Residuals - 'sigma' includes err on mean, align err, and additive err
        diffX = np.abs(x_pts - fitLineX)
        diffY = np.abs(y_pts - fitLineY)
        #diff = np.hypot(diffX, diffY)
        diff = (diffX + diffY) / 2.0
        re_pts = np.sqrt((diffX*xe_pts)**2 + (diffY*ye_pts)**2) / diff
        sigX = diffX / xe_pts
        sigY = diffY / ye_pts
        sig = diff / re_pts
        #print nCnt, starName,mag.mean(), diff.mean()

        # Combine this stars information with all other stars.
        sigmaX = np.concatenate((sigmaX, sigX))
        sigmaY = np.concatenate((sigmaY, sigY))
        sigma = np.concatenate((sigma, sig))
        diffX_all = np.concatenate((diffX_all,diffX))
        diffY_all = np.concatenate((diffY_all,diffY))
        diff_all = np.concatenate((diff_all,diff))
        mag_all = np.concatenate((mag_all,mag))
        # error on mean from 3 submaps, does not include additive error:
        xeom_all = np.concatenate((xeom_all,errx_p_allEp))
        yeom_all = np.concatenate((yeom_all,erry_p_allEp))
        # align errors:
        xaln_all = np.concatenate((xaln_all,errx_a_allEp))
        yaln_all = np.concatenate((yaln_all,erry_a_allEp))

        # Concatenate the stars' names so we can group by star if needed
        for nn in range(len(diff)):
            starName_all = np.concatenate([starName_all,[starName]])
        

    # Get average of all the errors
    eom_all = (xeom_all + yeom_all) / 2.0
    aln_all = (xaln_all + yaln_all) / 2.0

    # Residuals should have a gaussian probability distribution
    # with a mean of 0 and a sigma of 1. Overplot this to be sure.
    ggx = np.arange(-7, 7, 0.25)
    ggy = normpdf(ggx, 0, 1)

    # Get average residual and average error per magnitude bin
    magStep = 1.0
    magBins = np.arange(10.0, 20.0, magStep)
    residMag = np.zeros(len(magBins), float)
    errMag = np.zeros(len(magBins), float)
    alnMag = np.zeros(len(magBins), float)
    print ''
    print '%4s  %8s  %9s  %9s  %11s' % \
          ('Mag', 'Residual', 'RMS Error', 'Aln Error','Total Error')

    for mm in range(len(magBins)):
        mMin = magBins[mm] - (magStep / 2.0)
        mMax = magBins[mm] + (magStep / 2.0)
        idx = (np.where((mag_all >= mMin) & (mag_all < mMax)))[0]

        # Use the average X and Y residuals
        if (len(idx) > 0):
            residMag[mm] = np.abs(np.median(diff_all[idx])*1e3)
            errMag[mm] = np.median(eom_all[idx])
            alnMag[mm] = np.median(aln_all[idx])
        if additiveErr==True:
            print '%4.1f  %8.3f  %9.3f  %9.3f  %11.3f ' % \
                  (magBins[mm], residMag[mm], errMag[mm],alnMag[mm],
                   np.sqrt(errMag[mm]**2+addErr**2+alnMag[mm]**2))
        else:
            print '%4.1f  %8.4f  %9.4f  %9.4f  %11.4f' % \
                  (magBins[mm], residMag[mm], errMag[mm],alnMag[mm],
                   np.sqrt(errMag[mm]**2+alnMag[mm]**2))


    # Sum all the error sources
    err_tot = np.sqrt(errMag**2 + addErr**2 + alnMag**2)

    # Get typical alignment error
    med_alnErr = np.median(aln_all)

    print ''
    print 'Median residual for K<14.5: %5.3f mas' % np.median(residMag[0:5])
    print 'Median RMS error for K<14.5: %5.3f mas' % np.median(errMag[0:5])
    print 'Median align error for K<14.5: %5.3f mas' % np.median(alnMag[0:5])
    if additiveErr==True:
        print 'Median quad sum of RMS err, align err, and additive err for K<14.5: %5.3f mas' % \
              np.median(np.sqrt(errMag**2+alnMag**2+addErr**2)[0:5])
        print ''
        print 'Median quad sum of RMS err and additive err for K<14.5: %5.3f mas' % \
              np.median(np.sqrt(errMag**2+addErr**2)[0:5])
    print ''


    # Also look at the errors in magnitude bins with equal number of stars,
    # instead of equally-spaced magnitude bins
    midx = mag_all.argsort()
    nMeasurements = 500
    magN = mag_all[midx]
    starN = starName_all[midx]
    eomN = eom_all[midx]
    alnN = aln_all[midx]
    diffN = diff_all[midx]

    Nbins = np.arange(0, len(midx)+nMeasurements, nMeasurements)
    residNmag = np.zeros(len(Nbins), float)
    errNmag = np.zeros(len(Nbins), float)
    alnNmag = np.zeros(len(Nbins), float)
    nMagBins = np.zeros(len(Nbins), float)
    numStars = np.zeros(len(Nbins), float)
    reqNmag = np.zeros(len(Nbins), float)
    print '%4s  %7s  %8s  %9s  %9s  %11s  %11s' % \
          ('#Mag', '# stars', 'Residual', 'RMS Error', 'Aln Error','Total Error','Resid-Total')
    for nn in range(1, len(Nbins)):
        nMin = Nbins[nn-1]
        nMax = Nbins[nn]
        residNmag[nn] = np.abs(np.median(diff_all[nMin:nMax])*1e3)
        errNmag[nn] = np.median(eom_all[nMin:nMax])
        alnNmag[nn] = np.median(aln_all[nMin:nMax])

        nMagBins[nn] = (magN[nMin:nMax]).mean()

        # How many stars on average go into these nMeasurements?
        numStars[nn] = len(np.unique(starN[nMin:nMax]))

        # What additive error is required to get quad sum of all errors
        # to total the residuals?
 #       if (residNmag[nn]**2 > (alnNmag[nn]**2 + errNmag[nn]**2)):
        reqNmag[nn] = np.sqrt(residNmag[nn]**2 - errNmag[nn]**2 - alnNmag[nn]**2)
 #       else:
 #           reqNmag[nn] = 0.0

        print '%4.1f  %7i  %8.4f  %9.4f  %9.4f  %11.4f  %11.4f' % \
              (nMagBins[nn], numStars[nn], residNmag[nn], errNmag[nn],alnNmag[nn],
               np.sqrt(errNmag[nn]**2+alnNmag[nn]**2), reqNmag[nn])
              
        
    residNmag = residNmag[np.nonzero(residNmag)[0]]
    errNmag = errNmag[np.nonzero(errNmag)[0]]
    alnNmag = alnNmag[np.nonzero(alnNmag)[0]]
    reqNmag = reqNmag[np.nonzero(reqNmag)[0]]
    nMagBins = nMagBins[np.nonzero(nMagBins)[0]]

    ##########
    # Plot
    ##########
    binsIn = np.arange(-7, 7, 0.5)
    py.clf()
    py.figure(figsize=(6,6))
    (nx, bx, ptx) = py.hist(sigmaX,bins=binsIn,color='k',histtype='step',linewidth=1)
    ggamp = ((sort(nx))[-2:]).sum() / (2.0 * ggy.max())
    plot(ggx, ggy*ggamp, 'k-', ms=5)
    py.xlabel('X Residuals (sigma)',fontsize=12)
    py.ylabel('N',fontsize=12)
    py.savefig(home+'plots/residualsX.png')

    py.clf()
    (ny, by, pty) = py.hist(sigmaY,bins=binsIn,color='k',histtype='step',linewidth=1)
    ggamp = ((sort(ny))[-2:]).sum() / (2.0 * ggy.max())
    plot(ggx, ggy*ggamp, 'k-', ms=5)
    py.xlabel('Y Residuals (sigma)',fontsize=12)
    py.ylabel('N',fontsize=12)
    py.savefig(home+'plots/residualsY.png')

    py.clf()
    sigmaA = []
    for ss in range(len(sigmaX)):
        sigmaA = np.concatenate([sigmaA,[sigmaX[ss]]])
        sigmaA = np.concatenate([sigmaA,[sigmaY[ss]]])
    (na, ba, pa) = py.hist(sigmaA, bins=binsIn,color='k',histtype='step',linewidth=1)
    ggamp = ((sort(na))[-2:]).sum() / (2.0 * ggy.max())
    py.plot(ggx, ggy*ggamp, 'k-')
    py.xlabel('Residuals (sigma)')
    py.ylabel('N',fontsize=12)
    py.savefig(outdir+'residualsXY.eps')
    py.savefig(outdir+'residualsXY.png')

    # Put them all in one plot
    py.clf()
    py.figure(figsize=(7,7))
    py.subplot(3, 1, 1)
    (nx, bx, px) = py.hist(sigmaX, bins=binsIn, color='k',histtype='step',linewidth=1)
    ggamp = ((sort(nx))[-2:]).sum() / (2.0 * ggy.max())
    py.plot(ggx, ggy*ggamp, 'k-')
    py.xlabel('X Residuals (sigma)')

    py.subplot(3, 1, 2)
    (ny, by, pty) = py.hist(sigmaY, bins=binsIn, color='k',histtype='step',linewidth=1)
    ggamp = ((sort(ny))[-2:]).sum() / (2.0 * ggy.max())
    py.plot(ggx, ggy*ggamp, 'k-')
    py.xlabel('Y Residuals (sigma)')

    py.subplot(3, 1, 3)
    (ny, by, pty) = py.hist(sigma, bins=np.arange(0, 7, 0.5),color='k',histtype='step',linewidth=1)
    py.xlabel('Total Residuals (sigma)')

    py.subplots_adjust(wspace=0.34, hspace=0.33, right=0.95, top=0.97)
    py.savefig(outdir+'residualsDistribution.png')

    # Also plot the various sources of error (averaged across epochs)
    # as a function of K mag
    # "Residuals" here is the offset in mas from the best fit
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=14)
    py.figure(figsize=(6,6))
    py.clf()
    py.plot(magBins, errMag, 'b-.', label=r'$\bf{\sigma_{rms}}$', lw=3)
    if additiveErr == True:
        py.plot([magBins.min(),magBins.max()], [addErr,addErr],'g--',lw=2,label=r'$\bf{\sigma_{add}}$')
    #py.plot([magBins.min(),magBins.max()], [med_alnErr,med_alnErr],'c--',lw=2,label=r'$\bf{\sigma_{aln}}$')
    py.plot(magBins, alnMag, 'c--', lw=2, label=r'$\bf{\sigma_{aln}}$')
    py.plot(magBins, err_tot, 'k-', label=r'$\bf{\sigma_{tot}}$')
    py.plot(magBins, residMag, 'r:', label='Residuals', lw=3)
    py.xlabel('K Magnitude', fontsize=16)
    py.ylabel('Positional Uncertainty (mas)', fontsize=16)
    py.legend(prop=prop,loc=2,markerscale=0.5)
    py.savefig(home + 'plots/aveErrVsMag%s.png' % suffix)
    py.savefig(home + 'plots/aveErrVsMag%s.eps' % suffix)

    # Plot residuals vs. K mag
    py.clf()
    py.semilogy(mag_all, diff_all*1e3,'k.',ms=0.5)
    py.xlabel('K Magnitude')
    py.ylabel('Residuals (mas)')
    py.savefig(outdir+'residVsKmag.png')


    # Also plot the various sources of error (averaged across epochs)
    # as a function of K mag, but binned by number of stars
    py.figure(figsize=(6,6))
    py.clf()
    py.plot(nMagBins, errNmag, 'b-.', label=r'$\bf{\sigma_{rms}}$', lw=3)
    #if additiveErr == True:
    # only plot the reqNmag where they are 'nan'
    foo = np.where(np.isnan(reqNmag) != True)[0]
    py.plot(nMagBins[foo], reqNmag[foo],'g--',lw=2,label=r'$\bf{\sigma_{add}}$')
    py.plot(nMagBins, alnNmag, 'c--', lw=2, label=r'$\bf{\sigma_{aln}}$')
    #py.plot(nMagBins, np.sqrt(errNmag**2+alnNmag**2), 'k-', label=r'$\bf{\sigma_{tot}}$')
    py.plot(nMagBins, residNmag, 'r:', label='Residuals', lw=3)
    py.xlabel('K Magnitude', fontsize=16)
    py.ylabel('Positional Uncertainty (mas)', fontsize=16)
    py.legend(prop=prop,loc=2,markerscale=0.5)
    py.savefig(home + 'plots/aveErrVsMag_Nbins%s.png' % suffix)
    py.savefig(home + 'plots/aveErrVsMag_Nbins%s.eps' % suffix)



def sigmaVsEpoch(alnDir='13_08_21/',
                 align='align/align_d_rms_1000_abs_t',
                 poly='polyfit_c/fit', useAccFits=False):
    """
    Plot the average offset (in sigma) from the best fit
    velocity as a function of epoch.
    """
    outdir = home + 'plots/'

    s = starset.StarSet(root + alnDir + align, relErr=1)
    s.loadPolyfit(root + alnDir + poly, accel=0, arcsec=1)
    s.loadPolyfit(root + alnDir + poly, accel=1, arcsec=1)

    numEpochs = len(s.stars[0].years)
    
    # Use only stars detected in all epochs
    epochCnt = s.getArray('velCnt')
    idx = (np.where(epochCnt == epochCnt.max()))[0]
    newStars = []
    for ii in idx:
        newStars.append(s.stars[ii])
    s.stars = newStars
    
    print 'Using %d out of %d stars detected in %d epochs' % \
          (len(newStars), len(epochCnt), epochCnt.max())

    # Make some empty arrays to hold all our results.
    sigmaX = np.zeros(numEpochs, float)
    sigmaY = np.zeros(numEpochs, float)
    sigma  = np.zeros(numEpochs, float)
    diffEpX = np.zeros(numEpochs, float)
    diffEpY = np.zeros(numEpochs, float)
    diffEp  = np.zeros(numEpochs, float)

    # Fetch the fit parameters for all the stars
    if (useAccFits == True):
        fitVarX = 'fitXa'
        fitVarY = 'fitYa'
    else:
        fitVarX = 'fitXv'
        fitVarY = 'fitYv'

    t0 = s.getArray(fitVarX + '.t0')
    x0 = s.getArray(fitVarX + '.p')
    vx = s.getArray(fitVarX + '.v')
    y0 = s.getArray(fitVarY + '.p')
    vy = s.getArray(fitVarY + '.v')

    x0e = s.getArray(fitVarX + '.perr')
    y0e = s.getArray(fitVarY + '.perr')
    vxe = s.getArray(fitVarX + '.verr')
    vye = s.getArray(fitVarY + '.verr')

    if (useAccFits == True):
        ax = s.getArray(fitVarX + '.a')
        ay = s.getArray(fitVarY + '.a')
        axe = s.getArray(fitVarX + '.aerr')
        aye = s.getArray(fitVarY + '.aerr')

    # Loop through all the epochs and determine average residuals
    for ee in range(numEpochs):
        # Observed data
        x = s.getArrayFromEpoch(ee, 'x')
        y = s.getArrayFromEpoch(ee, 'y')
        xerr_p = s.getArrayFromEpoch(ee, 'xerr_p')
        yerr_p = s.getArrayFromEpoch(ee, 'yerr_p')
        xerr_a = s.getArrayFromEpoch(ee, 'xerr_a')
        yerr_a = s.getArrayFromEpoch(ee, 'yerr_a')
        xerr = hypot(xerr_p, xerr_a)
        yerr = hypot(yerr_p, yerr_a)
        t = s.stars[0].years[ee]

        dt = t - t0
        fitLineX = x0 + (vx * dt)
        fitSigX = np.sqrt( x0e**2 + (dt * vxe)**2 )

        fitLineY = y0 + (vy * dt)
        fitSigY = np.sqrt( y0e**2 + (dt * vye)**2 )

        if (useAccFits == True):
            fitLineX += (ax * dt**2) / 2.0
            fitLineY += (ay * dt**2) / 2.0
            fitSigX = np.sqrt(fitSigX**2 + (dt**2 * axe / 2.0)**2)
            fitSigY = np.sqrt(fitSigY**2 + (dt**2 * aye / 2.0)**2)

        # Residuals
        diffX = x - fitLineX
        diffY = y - fitLineY
        #diff = hypot(diffX, diffY)
        diff = (diffX + diffY) / 2.0
        rerr = np.sqrt((diffX*xerr)**2 + (diffY*yerr)**2) / diff
        sigX = diffX / xerr
        sigY = diffY / yerr
        sig = diff / rerr

        #print 'Epoch %d' % ee

        #print '%8.5f +/- %8.5f   %8.5f +/- %8.5f  %8.5f (%5.1f sigma)' % \
        #      (x[0], xerr[0], fitLineX[0], fitSigX[0], diffX[0], sigX[0])
        #print '%8.5f +/- %8.5f   %8.5f +/- %8.5f  %8.5f (%5.1f sigma)' % \
        #      (y[0], yerr[0], fitLineY[0], fitSigY[0], diffY[0], sigY[0])

        sigmaX[ee] = sigX.mean()
        sigmaY[ee] = sigY.mean()
        sigma[ee] = median(sig)

        diffEpX[ee] = diffX.mean()
        diffEpY[ee] = diffY.mean()
        diffEp[ee] = median(diff)

                        
    ##########
    # Plot
    ##########
    clf()
    years = s.stars[0].years
    plot(years, sigmaX, 'rx')
    plot(years, sigmaY, 'bx')
    plot(years, sigma, 'ko')
    xlabel('Epoch (years)')
    ylabel('Median Residual Error (sigma)')
    legend(('X', 'Y', 'Total'))

    #savefig(outdir+'residualsVsEpoch.eps')
    savefig(outdir+'residualsVsEpoch.png')

    clf()
    plot(years, diffEpX*1000.0, 'rx')
    plot(years, diffEpY*1000.0, 'bx')
    plot(years, diffEp*1000.0, 'ko')
    xlabel('Epoch (years)')
    ylabel('Median Residual Error (mas)')
    legend(('X', 'Y', 'Total'))
    #savefig(outdir+'residualsVsEpochMAS.eps')
    savefig(outdir+'residualsVsEpochMAS.png')


    # Print out epochs with higher than 3 sigma median residuals
    hdx = (np.where(sigma > 3))[0]
    print 'Epochs with median residuals > 3 sigma:'
    for hh in hdx:
        print '%8.3f  residual = %4.1f' % (s.stars[0].years[hh], sigma[hh])


def chi2distrib(alnDir='13_08_21/',
                align='align/align_d_rms_1000_abs_t', poly='polyfit_c/fit',
                points='points_c/',starlist='all',epochCut=10):

    outdir = home + 'plots/'
    
    tag = ''
    if starlist=='all':
        tag = '_all'
    if starlist=='oldstars':
        tag = '_old'
    if starlist=='yngstars':
        tag = '_yng'

    s = loadPop(alnDir=alnDir,align=align,poly=poly,points=points,starlist=starlist)

    names = s.getArray('name')

    xchi2 = s.getArray('fitXv.chi2')
    ychi2 = s.getArray('fitYv.chi2')
    xchi2r = s.getArray('fitXv.chi2red')
    ychi2r = s.getArray('fitYv.chi2red')
    mag = s.getArray('mag')

    pntcnt = np.zeros(len(names))
    for ii in range(len(names)):
        pntFileName = '%s%s%s%s.points' % (root, alnDir, points, names[ii])
        pntFile = open(pntFileName)
        data = pntFile.readlines()
        pntcnt[ii] = len(data)
    # how many epochs are there?
    numEpochs = np.int(pntcnt[0]) # 16C should be in all epochs

    idx = np.where((pntcnt >= epochCut) & (mag < 16))[0]
    xchi2 = xchi2[idx]
    ychi2 = ychi2[idx]
    xchi2r = xchi2r[idx]
    ychi2r = ychi2r[idx]
    mag = mag[idx]
    pntcnt = pntcnt[idx]

    binsInX = np.arange(0.0, max(xchi2)+1, 2.0)
    binsInY = np.arange(0.0, max(ychi2)+1, 2.0)
    # Expected chi2 for n dof:
    dof = numEpochs - 2
    xpctdX = stats.chi2.pdf(binsInX,dof)

    # Only plot the stars detected in N=44 epochs
    all = np.where(pntcnt > 30)[0]
    print '%i stars with detections in >40 epochs' % len(all)
    py.clf()
    py.figure(figsize=(6,6))
    (nx,bx,p1x) = py.hist(xchi2[all],bins=binsInX,color='r',histtype='step',linewidth=2,label='X',normed=True)
    (ny,by,p1y) = py.hist(ychi2[all],bins=binsInY,color='b',histtype='step',linewidth=2,label='Y',normed=True)
    #print np.sum(nx*np.diff(bx))
    #print np.sum(ny*np.diff(by))
    py.plot(binsInX,xpctdX,'k--')
    py.xlabel('Total Chi Squared',fontsize=12)
    py.ylabel('N',fontsize=12)
    py.legend(numpoints=1)
    py.axis([0,100,0,0.03])
    py.savefig(home+'plots/hist_chi2%s.png' % tag)


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


def getOrbElements(alnDir='09_09_20/',
                   starType='all'):
    """
    For the stars with orbits from efit, read their .param output
    files and get their orbital elements
    """

    _list = asciidata.open(root + alnDir + 'efit/names.dat')
    names = _list[0]._data
    _type = _list[1].tonumpy()

    idx = []
    if starType == 'yng':
        for i in range(len(names)):
            if (_type[i] == 1):
                idx = np.concatenate([idx, [i]])
    elif starType == 'old':
        for i in range(len(names)):
            if (_type[i] == 0):
                idx = np.concatenate([idx, [i]])
    elif starType == 'all':
        idx = range(len(names))
    names = np.array([names[int(ii)] for ii in idx])
    _type = np.array([_type[int(ii)] for ii in idx])

    numStars = len(names)
    r0 = np.zeros(numStars)
    a = np.zeros(numStars)
    p = np.zeros(numStars)
    e = np.zeros(numStars)
    t0 = np.zeros(numStars)
    smallOmega = np.zeros(numStars)
    i = np.zeros(numStars)
    bigOmega = np.zeros(numStars)
    aUp = np.zeros(numStars)
    pUp = np.zeros(numStars)
    eUp = np.zeros(numStars)
    t0Up = np.zeros(numStars)
    smallOmegaUp = np.zeros(numStars)
    iUp = np.zeros(numStars)
    bigOmegaUp = np.zeros(numStars)
    aLo = np.zeros(numStars)
    pLo = np.zeros(numStars)
    eLo = np.zeros(numStars)
    t0Lo = np.zeros(numStars)
    smallOmegaLo = np.zeros(numStars)
    iLo = np.zeros(numStars)
    bigOmegaLo = np.zeros(numStars)

    search = 2

    outFile = root + alnDir + 'efit/new_orbits.dat'
    out = open(outFile, 'w')
    hdr = '%10s  %6s  %6s  %10s  %7s  %6s  %6s  %6s  %5s\n'
    fmt = '%10s  %6.2f  %6.1f  %10.4f  %7.4f  %6.1f  %6.1f  %6.1f  %5i\n'
    out.write(hdr % ('#Star','P','A','t0','e','i','Omega','omega','search'))
    out.write(hdr % ('#Name','(yrs)','(mas)','(yrs)','()','(deg)','(deg)','(deg)','(pix)'))

    outFile2 = root + alnDir + 'efit/cntrl_rva.dat'
    out2 = open(outFile2, 'w')
    out2.write('%15s    %5s     %12s      %14s     %12s  %8s\n' % \
          ('Name','r (")','v2d (mas/yr)','a_rad (mas/yr^2)','a_rad (sigma)','t0'))

    s = loadPop(alnDir=alnDir,starlist='all',align='align/align_p_rms_1000_abs_t',
                poly='polyfit_p/fit')
    allNames = s.getArray('name')

    # In mas, mas/yr, mas/yr^2
    x0 = s.getArray('x0')
    y0 = s.getArray('y0')
    x0e = s.getArray('x0e')
    y0e = s.getArray('y0e')
    vx = s.getArray('vx')
    vy = s.getArray('vy')
    v2d = np.sqrt(vx**2 + vy**2)
    vxe = s.getArray('vxe')
    vye = s.getArray('vye')
    ax = s.getArray('ax')
    ay = s.getArray('ay')
    a2d = np.sqrt(ax**2 + ay**2)
    axe = s.getArray('axe')
    aye = s.getArray('aye')
    cnt = s.getArray('cnt')
    t0 = s.getArray('t0x')
    mag = s.getArray('mag')

    # Get accelerations in the radial/tangential direction
    r = np.sqrt(x0**2 + y0**2)
    ar = ((ax*x0) + (ay*y0)) / r
    at = ((ax*y0) - (ay*x0)) / r
    are =  (axe*x0/r)**2 + (aye*y0/r)**2
    are += (y0*x0e*at/r**2)**2 + (x0*y0e*at/r**2)**2
    are =  np.sqrt(are)
    ate =  (axe*y0/r)**2 + (aye*x0/r)**2
    ate += (y0*x0e*ar/r**2)**2 + (x0*y0e*ar/r**2)**2
    ate =  np.sqrt(ate)

    r_stars = []
    orb_mags = []

    print '%15s   %5s    %5s     %12s      %14s     %12s  %8s' % \
          ('Name','KMag','r (")','v2d (mas/yr)','a_rad (mas/yr^2)','a_rad (sigma)','t0')
    for ss in range(numStars):
        _fileo = open(root + alnDir + 'efit/orbit.' + names[ss] + '.output.param1')
        _file = _fileo.readline().split()
        _elem = open(root + alnDir + 'efit/orbit.' + names[ss] + '.output.error')
        elem = _elem.readline().split()

        r0[ss] = float(_file[0])
        a[ss] = float(elem[0])
        aUp[ss] = float(elem[1])
        aLo[ss] = float(elem[2])
        p[ss] = float(elem[3])
        pUp[ss] = float(elem[4])
        pLo[ss] = float(elem[5])
        e[ss] = float(elem[6])
        eUp[ss] = float(elem[7])
        eLo[ss] = float(elem[8])
        t0[ss] = float(elem[9])
        t0Up[ss] = float(elem[10])
        t0Lo[ss] = float(elem[11])
        smallOmega[ss] = float(elem[12])
        smallOmegaUp[ss] = float(elem[13])
        smallOmegaLo[ss] = float(elem[14])
        i[ss] = float(elem[15])
        iUp[ss] = float(elem[16])
        iLo[ss] = float(elem[17])
        bigOmega[ss] = float(elem[18])
        bigOmegaUp[ss] = float(elem[19])
        bigOmegaLo[ss] = float(elem[20])

        # Write it out to an orbits.dat type file
        out.write(fmt % (names[ss], p[ss], a[ss], t0[ss], e[ss], i[ss],
                         bigOmega[ss], smallOmega[ss], search))

        # Get other information on these stars (r, v, a, etc)
        nn = allNames.index(names[ss])
        print '%15s   %5.2f   %5.3f"   %8.3f mas/yr   %9.3f mas/yr^2   %8.2f sigma    %8.3f' % \
              (names[ss], mag[nn], r[nn], v2d[nn]*1e3, ar[nn]*1e3, ar[nn]/are[nn], t0[nn])
        out2.write('%15s   %5.3f"   %8.3f mas/yr   %9.3f mas/yr^2   %8.2f sigma    %8.3f\n' % \
              (names[ss], r[nn], v2d[nn]*1e3, ar[nn]*1e3, ar[nn]/are[nn], t0[nn]))

        # Get the r2d's for these stars:
        r_stars = np.concatenate([r_stars,[r[nn]]])

        orb_mags = np.concatenate([orb_mags,[mag[nn]]])
        
    out.close()
    out2.close()

    # Convert semi-major axis from mas to pc
    a *= dist / (1.e3*206265.)
    aUp *= dist / (1.e3*206265.)
    aLo *= dist / (1.e3*206265.)

    # Make errors symmetric by taking average of upper and lower errors
    a_err = (aUp + aLo) / 2.0
    p_err = (pUp + pLo) / 2.0
    e_err = (eUp + eLo) / 2.0
    t0_err = (t0Up + t0Lo) / 2.0
    smallOmega_err = (smallOmegaUp + smallOmegaLo) / 2.0
    i_err = (iUp + iLo) / 2.0
    bigOmega_err = (bigOmegaUp + bigOmegaLo) / 2.0

    # Get peri- and apo- distances 
    amin = a * (1. - e)
    amax = a * (1. + e)

    inner = np.where(r_stars < 0.8)[0]
    outer = np.where(r_stars >= 0.8)[0]
    yng = np.where(_type == 1)[0]
    old = np.where(_type == 0)[0]
    unk = np.where(_type < 0)[0]

    iDisk = 115.
    iDiske = 3.
    WDisk = 100.
    WDiske = 3.
    iSig = [np.hypot(iDiske, ii) for ii in i_err]
    WSig = [np.hypot(WDiske, ii) for ii in bigOmega_err]
    onD = np.where((abs(i - iDisk) < iSig) | (abs(bigOmega - WDisk) < WSig))[0]

    fig = py.figure(1)
    fig.clear()
    fig.hold(True)
    ax = fig.add_subplot(111)
    #loga = np.log10(a)
    #logaEU = np.log10(a+a_err)-loga
    #logaEL = -np.log10(a-a_err)+loga
    #logaE = np.array((logaEL, logaEU))
    #logaE_old = np.array((logaEL[old], logaEU[old]))
    #logaE_yng = np.array((logaEL[yng], logaEU[yng]))
    #logaE_unk = np.array((logaEL[unk], logaEU[unk]))
    p1 = ax.errorbar(a[old],e[old],xerr=a_err[old],yerr=e_err[old],fmt='r.',ms=8)
    p2 = ax.errorbar(a[yng],e[yng],xerr=a_err[yng],yerr=e_err[yng],fmt='b.',ms=8)
    p3 = ax.errorbar(a[unk],e[unk],xerr=a_err[unk],yerr=e_err[unk],fmt='k.',ms=8)
    #p1 = py.errorbar(loga[old],e[old],xerr=logaE_old,yerr=e_err[old],fmt='r.',ms=8)
    #p2 = py.errorbar(loga[yng],e[yng],xerr=logaE_yng,yerr=e_err[yng],fmt='b.',ms=8)
    #p3 = py.errorbar(loga[unk],e[unk],xerr=logaE_unk,yerr=e_err[unk],fmt='k.',ms=8)
    vi = ax.xaxis.get_view_interval().astype(np.int)
    ax.set_xscale('log')
    py.xlabel('Semi-major Axis (pc)')
    py.ylabel('Eccentricity')
    lgd1 = 'Old'
    lgd2 = 'Young'
    lgd3 = 'Unknown'
    lgd = ax.legend((p1[0], p2[0], p3[0]), (lgd1,lgd2,lgd3), numpoints=1, fancybox=True, loc=2,
                    labelspacing=0.01,handletextpad=0.01,borderpad=0.1,handlelength=0.75)
    for t in lgd.get_texts():
        t.set_fontsize('x-small')
    py.savefig(home + 'plots/cntrl_eccVsA_test.png')
    py.close(1)

    fig = py.figure(2)
    fig.clear()
    fig.hold(True)
    usetexTrue()
    ax = fig.add_subplot(111)
    pIn = ax.errorbar(i[inner],bigOmega[inner],xerr=i_err[inner],yerr=bigOmega_err[inner],fmt='r.')
    pOut = ax.errorbar(i[outer],bigOmega[outer],xerr=i_err[outer],yerr=bigOmega_err[outer],fmt='b.')
    # Also plot a point for the i,Omega of the disk (per JLu)
    ax.errorbar([115],[100],yerr=3,xerr=3,fmt='ko',ms=5)
    lgda = r'r $<$ 0.8"'
    lgdb = r'r $\geq$ 0.8"'
    lgd2 = ax.legend((pIn[0], pOut[0]), (lgda,lgdb), numpoints=1)
    lgdLines = lgd2.get_lines()
    py.xlabel('Inclination (deg)')
    py.ylabel(r'$\Omega$ (deg)')
    ax.axis([0,180,0,360])
    py.savefig(home + 'plots/cntrl_OmegaInc.png')
    py.close(2)

    py.figure(3)
    py.clf()
    py.plot(amin,amax,'k.')
    py.xlabel('Periastron Distance (pc)')
    py.ylabel('Apoastron Distance (pc)')
    py.savefig(home + 'plots/cntrl_AminVsAmax.png')
    py.close(3)

    py.figure(4)
    py.clf()
    binsIn = np.arange(0.0001, 1.0, 0.05)
    py.figure(figsize=(6,6))
    (bins,data)=histNofill.hist(binsIn,e)
    py.plot(bins,data)
    py.xlabel('Eccentricity',fontsize=12)
    py.ylabel('N',fontsize=12)
    py.axis([0,1,0,max(data)+1])
    py.savefig(home + 'plots/hist_ecc.png')
    py.close(4)

    # Also plot up the cumulative distribution function for eccentricities
    e_sort = sort(e)
    cdf = range(1, len(e) + 1)
    cdf = [cdf[cc] / (len(e)*1.0) for cc in range(len(cdf))]
    e0 = [1.0, 2.0]
    ecc_err = array([1.0*len(e_sort)])
    eccfit = fitPowerLawMP(e0,[e_sort, cdf, ecc_err],1)
    eccParams = eccfit.params
    eline = np.arange(0.01, 1.0, 0.01)
    cdfline = eccParams[0]*eline**eccParams[1]
    py.figure(5)
    py.clf()
    py.plot(e_sort, cdf,'k.')
    py.plot(eline,cdfline,'k-')
    py.xlabel('Eccentricity',fontsize=12)
    py.ylabel('CDF',fontsize=12)
    py.text(0.05,0.7, r'cdf ${\sim} e^{%4.2f}$' % eccParams[1])
    py.axis([0,1,0,1])
    py.savefig(home + 'plots/cdf_ecc.png')
    usetexFalse()
    py.close(5)


    py.figure(6)
    py.clf()
    p1 = py.errorbar(a[old],orb_mags[old],xerr=a_err[old],fmt='r.',ms=15)
    p2 = py.errorbar(a[yng],orb_mags[yng],xerr=a_err[yng],fmt='b.',ms=15)
    p3 = py.errorbar(a[unk],orb_mags[unk],xerr=a_err[unk],fmt='k.',ms=15)
    py.xlabel('Semi-major Axis (pc)')
    py.ylabel('K Magnitude')
    lgd1 = 'Old'
    lgd2 = 'Young'
    lgd3 = 'Unknown'
    lgd = py.legend((p1, p2, p3), (lgd1,lgd2,lgd3), numpoints=1)
    lgdLines = lgd.get_lines()
    lgdLines[0].set_marker('.')
    lgdLines[1].set_marker('.')
    lgdLines[2].set_marker('.')
    lgdLines[0].set_mfc('r')
    lgdLines[1].set_mfc('b')
    lgdLines[2].set_mfc('k')
    lgdLines[0].set_ms(12)
    lgdLines[1].set_ms(12)
    lgdLines[2].set_ms(12)
    py.savefig(root + alnDir + 'plots/mag_A.png')

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
              chainsDir = 'efit/chains/', starlist='all', 
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
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, root_tmp=root_tmp, alnDir=alnDir, chainsDir = chainsDir)

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


def makeNSFtable(alnDir='14_06_18/', root_tmp='/g/ghez/align/',
                 align = 'align/align_d_rms_1000_abs_t', updateErr = True,radCut = 100.0,
                 poly='polyfit_nz/fit', points='points_nz/', sigma=3.0, polyj='polyfit_nzj/fit',
                 f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',
                 starlist='all',magCut=22, nEpochs=14):
    
    names,x0,y0,mag,r,ar,are,sigmaR,z0,z0e,az,aze,r3d,r3de,mext,cnt,vz_cms,vze_cms,pass_fj,jnames,jx,jy,jxe,jye,jxsig,jysig,xchi2r,ychi2r=histAccel(alnDir=alnDir,root_tmp=root_tmp,align=align,updateErr=updateErr,radCut=radCut,poly=poly,points=points,sigma=sigma,polyj=polyj,f_test=f_test,pvalue=pvalue,chainsDir=chainsDir,starlist=starlist,magCut=magCut,nEpochs=nEpochs,verbose=False,outTolatex=True)
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
              f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',
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
             f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',
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
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, root_tmp=root_tmp, alnDir=alnDir, chainsDir = chainsDir)

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
                 f_test=True, pvalue=4.0, chainsDir = 'efit/chains/',outTolatex=False,outToAL=False,
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
                 chainsDir = 'efit/chains/', slope=[-10.0,10.0,0.1],inter=[0.0,10.0,0.1],
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
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, mag, root_tmp=root_tmp, alnDir=alnDir, chainsDir = chainsDir)
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

    #maxR = np.max(r3d_pc[fmdex])+0.000001
    #minR = np.min(r3d_pc[fmdex])-0.000001
    #binSize = (maxR-minR)/binNum
    #radLim = np.array([minR+binSize*i for i in range(binNum+1)])
#    radLim = np.array([binSize*i for i in range(binNum+1)])
    binNum = 3
    maxR = 1.0
    minR = 0.01
    radLim = np.array([0.01,0.025,0.05,0.075])
    number = np.zeros(binNum)
    volume = np.zeros(binNum)
    radius3D = np.zeros(binNum)
    binSize = np.zeros(binNum)

    for i in range(binNum):
        pdb.set_trace()
#        tmpDex = np.where((np.log10(r3d_pc) >= radLim[i]) & (np.log10(r3d_pc) < radLim[i+1]) & (sig_r3d >= sig3d) & (r3de_pc < errCut))[0]
        tmpDex = np.where((oldProb > 0.0) & (r3d_pc >= radLim[i]) & (r3d_pc < radLim[i+1]) & (sig_r3d >= sig3d) & (r3de_pc < errCut))[0]
        tmpSum = oldProb[tmpDex]
        number[i] = tmpSum.sum()
#        volume[i] = 4.0*pi*((10.0**radLim[i+1])**3 - (10.0**radLim[i])**3)/3.0
        volume[i] = 4.0*pi*(radLim[i+1]**3 - radLim[i]**3)/3.0
        radius3D[i] = np.median((r3d_pc[tmpDex]))
        #pdb.set_trace()
        binSize[i] = radLim[i+1] - radLim[i]

    densityBin = number/volume #number density in pc^-3
    numerr = np.sqrt(number)
    err_pl = np.log10((numerr+number)/number)
    error_density = numerr/volume
    for i in range(binNum):
        if (error_density[i] >= densityBin[i]):
            error_density[i] = densityBin[i] - 1.0

#    pdb.set_trace()
    sval,ival=scLeast2(radius3D,densityBin,error=error_density,slope=slope,inter=inter,guess=guess)

    print number

    plOld = np.where((oldProb > 0.0) & (sig_r3d >= sig3d) & (r3de_pc < errCut) & (r3d_pc <= r3dMax))[0]
    ploz = np.zeros(len(plOld))+10.0**0.05
    
    rho_pl = np.log10(densityBin)
    r_pl = np.log10(radius3D)
    rLeft = np.array([radLim[i] for i in range(binNum)])
    #errpatch = np.array([[densityBin[0]-1.0,error_density[1],error_density[2]],error_density])
 #   fitLine = r_pl*sval + ival
    py.clf()
    py.bar(rLeft,densityBin,width=binSize,bottom=1.0,color='w')
    py.errorbar(radius3D,densityBin,yerr=error_density,fmt='o')
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
    like_final = np.exp(lnL - np.max(lnL))
    X,Y = np.meshgrid(gamma,(rcut/cm_in_pc))
    import scOld
    levels = scOld.getContourLevels(like_final)
    left,width = 0.1,0.65
    bottom,height = 0.1,0.65
    bottom_h = left_h = left+width+0.02
    rect_cont = [left,bottom,width,height]
    rect_histx = [left,bottom_h,width,0.2]
    rect_histy = [left_h,bottom,0.2,height]
    py.clf()
    axCont = py.axes(rect_cont)
    axHistx = py.axes(rect_histx)
    axHisty = py.axes(rect_histy)
    
    axCont.contour(Y,X,like_final,levels=levels,colors=('limegreen','gray','k'))
    py.xlabel('Radius Cut (pc)')
    py.ylabel('Gamma')
    py.savefig('/u/schappell/plots/single_power_like_contour.png')
    pdb.set_trace()




def likeGamma(alnDir='14_06_18/',root_tmp='/g/ghez/align/',updateErr=True,align='align/align_d_rms_1000_abs_t',
              poly='polyfit_nz/fit',points='points_nz/',R2Dcut = 1.7,nEpochs = 14.,magCut=15.5,globalt0 = 2013.318,dataCut=100.,
              starlist='all',chainsDir='efit/chains/',grange=[-14.0,2.0,100.],rrange=[1.71,10.0,100.],numz=100.):

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
    py.savefig('/u/schappell/plots/rough2D.png')
    #pdb.set_trace()

    if (updateErr == True):
        origin_val = asciidata.open('/g/ghez/align/14_06_18/' + chainsDir + 'efit_summary.txt')
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



def forShoko(stars):

    for i in range(len(stars)):
        tmpName = str(stars[i])
        pointsTab = asciidata.open('/g/ghez/align/14_06_18/points_nz/' + tmpName + '.points')
        xp = pointsTab[1].tonumpy()
        yp = pointsTab[2].tonumpy()

        radius = np.sqrt(xp**2+yp**2)
        print tmpName,np.min(radius)
