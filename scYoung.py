##---------------------------------------------------------------------------
## This code loads IPython but modifies a few things if it detects it's running
## embedded in another IPython session (helps avoid confusion)
#1
#try:
#    __IPYTHON__
#except NameError:
#    argv = ['']
#    banner = exit_msg = ''
#else:
#    # Command-line options for IPython (a list like sys.argv)
#    argv = ['-pi1','In <\\#>:','-pi2','   .\\D.:','-po','Out<\\#>:']
#    banner = '*** Nested interpreter ***'
#    exit_msg = '*** Back in main IPython ***'
#
## First import the embeddable shell class
#from IPython.Shell import IPShellEmbed
## Now create the IPython shell instance. Put ipshell() anywhere in your code
## where you want it to open.
#ipshell = IPShellEmbed(argv,banner=banner,exit_msg=exit_msg)
##---------------------------------------------------------------------------

from gcwork import objects
#from gcwork import starset
import sc_starset as starset
from gcwork import util
from gcwork import orbits
from pysqlite2 import dbapi2 as sqlite
from scipy import stats
from gcwork import starTables
from gcwork import young
import scipy
import pyfits
import nmpfit_sy
import asciidata, os, sys, pickle
from pylab import *
import numpy as np
import pylab as py
import math
import histNofill
import pdb

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

pi = math.pi

# Latest Mass and Ro 
mass = 4.6e6
dist = 8232.9
G = 6.6726e-8
msun = 1.99e33
GM = G * mass * msun
sec_in_yr = 3.1557e7
cm_in_au = 1.496e13
cm_in_pc = 3.086e18
km_in_pc = 3.086e13
au_in_pc = 206265.0
asy_to_kms = dist * cm_in_au / (1e5 * sec_in_yr)

rootDir = '/g/ghez/align/'
home = '/u/schappell/'
#rootDir = '/u/syelda/research/gc/aligndir/'
#rootDir = '/u/syelda/research/gc/aligndir/test_schoedel/'

def loadPop(alnDir='11_07_05/', align = 'align/align_d_rms_1000_abs_t',
            poly='polyfit_c/fit',points='points_c/', starlist='yngstars'):
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

    s = starset.StarSet(rootDir + alnDir + align)
    s.loadPolyfit(rootDir + alnDir + poly, accel=0, arcsec=1)
    s.loadPolyfit(rootDir + alnDir + poly, accel=1, arcsec=1)
    s.loadPoints(rootDir + alnDir + points)

    # Get the most up-to-date young stars from sql data base.
    yng = young.youngStarNames()

    # Pull out old stars
    yngstars = []
    names = [star.name for star in s.stars]

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

    # Build oldstars 
    for star in oldstars:
        # In arcsec, mas/yr, mas/yr^2
        # Stars may not have accel fit do try instead
        try:
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



def nzErr(xerr, yerr, vxerr, vyerr, year_x, year_y, alnDir = '13_08_21/', chainsDir = 'efit/chains_S0-2_newRV2/'):
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

    #Update errors
    xerr = np.sqrt(xerr**2 + ori_x0e**2 + ((year_x - t_0)*ori_vxe)**2)
    yerr = np.sqrt(yerr**2 + ori_y0e**2 + ((year_y - t_0)*ori_vye)**2)
    vxerr = np.sqrt(vxerr**2 + ori_vxe**2)
    vyerr = np.sqrt(vyerr**2 + ori_vye**2)

    return xerr, yerr, vxerr, vyerr




def plotAccelOnImage(alnDir='13_08_21/', align='align/align_d_rms_1000_abs_t',
              poly='polyfit_c/fit', points='points_c/',nEpochs=30,magCut=15,
              starlist='oldstars',suffix='_old_rmConfused',sigma=5.0, markRef=False,
              refFile='source_list/label_restrict_noAccelSources_fullAlign.dat',
              specMap=False, updateErr = True, chainsDir='efit/chains_S0-2_newRV2/'):
    """
    Plot positions of young stars over a GC image, with
    symbols indicating whether star has detectable accelerations.
    """

    s = loadPop(alnDir=alnDir,align=align,poly=poly,points=points,starlist=starlist)
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
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

    t0x = s.getArray('t0x')
    t0y = s.getArray('t0y')

    #Update errors in position
    if updateErr:
        x0e,y0e,vxe,vye = nzErr(x0e, y0e, vxe, vye, t0x, t0y, alnDir=alnDir, chainsDir = chainsDir)

    # Make an epochs cut
    #idx = np.where(cnt > nEpochs)[0]
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
#        atBins=np.array([9, 12.73, 13.78, 14.56, 15.18, 15.39, 15.595, 15.88, 17])
#        deltaArr=np.array([1.5302, 2.0025, 2.9809, 3.8496, 4.6642, 4.6273, 5.0453, 5.2388])*1e-5
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
#        for i in range(len(mag)):
#            for j in range(len(deltaAt)):
#                if ((mag[i] > atBins[j]) & (mag[i] <= atBins[j+1])):
#                    delta[i]=deltaAt[j]
#                    for k in range(len(deltaRa)):
#                        if ((r[i] > radBins[k]) & (r[i] <= radBins[k+1])):
#                            if deltaRa[k] > delta[i]:
#                                delta[i]=deltaRa[k]
        
        magBins = np.array([9,14.32,15.44,15.86,18])
        radBins = np.array([0.0,1.9434,2.9501,3.8977,6.5])
#        deltaArr = np.array([[14.601,36.187,10.737,19.217],[24.742,21.776,35.729,13.194],[22.587,41.706,44.864,12.631],
#                             [49.032,46.48,50.812,28.274]])*1e-6
        deltaArr = np.array([[11.49,2.85,11.09,28.78],[22.64,15.93,37.53,5.35],[22.59,25.97,44.86,42.86],
                             [49.43,51.14,6.54,28.28]])*1e-6
        for i in range(len(mag)):
            for j in range(len(magBins)-1):
                if ((mag[i] > magBins[j]) & (mag[i] <= magBins[j+1])):
                    idxm = j
            for p in range(len(radBins)-1):
                if ((r[i] > radBins[p]) & (r[i] <= radBins[p+1])):
                    idxr = p
            delta[i] = deltaArr[idxm,idxr]
#            print mag[i],r[i],delta[i]

        ate = np.sqrt(ate**2 + delta**2)
        are = np.sqrt(are**2 + delta**2)


    sigmaR = ar/are
    sigmaT = at/ate

    phys = np.where((ar/are) < -sigma)[0]
    posRad = np.where((ar/are) > (sigma))[0]
    tan = np.where(np.abs(at/ate) > (sigma))[0]

    # Total acceleration
    atot = py.hypot(ax, ay)
    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

    if markRef == True:
        # Read in label file and mark the reference stars on the image
        inFile = asciidata.open(rootDir + alnDir + refFile)
        inID = inFile[0].tonumpy()
        inX = inFile[2].tonumpy()
        inY = inFile[3].tonumpy()
        use = inFile[11].tonumpy()

        ref = np.where(use != '0')[0]
        inID = [inID[ii] for ii in ref]
        inX = inX[ref]
        inY = inY[ref]

    if specMap == True:
        imgFile = '/u/ghezgroup/data/gc/04jul/combo/mag04jul.fits'
        sgra = [490., 613.]
    else:
        imgFile='/u/ghezgroup/data/gc/09maylgs/combo/mag09maylgs_sgra_dim_kp.fits'
        sgra=[624.5,726.3]
#        imgFile = '/u/ghezgroup/data/gc/14maylgs2/combo/mag14maylgs2_kp.fits'
#        sgra = [627., 695.]

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

    hdr = '%15s  %5s  %6s  %17s  %14s  %17s  %14s  %8s  %5s %5s '
    fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %8.3f mas/yr^2  %8.2f sigma  %2d epochs  %5.2f  %5.2f'

    print
    print '%s SC Significant physical Accelerations' % (len(phys))
    print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)','a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs','X','Y')

    if specMap == True:
        ax.imshow(np.log10(img+1), aspect='equal', interpolation='bicubic',
                  extent=[max(xL), min(xL), min(yL), max(yL)],vmin=0.8,vmax=2.0,
                  origin='lowerleft', cmap=py.cm.gray_r)
    else:
        ax.imshow(np.log10(img+1), aspect='equal', interpolation='bicubic',
                  extent=[max(xL), min(xL), min(yL), max(yL)],vmin=2.2,vmax=5,
                  origin='lowerleft', cmap=py.cm.gray_r)
    for ii in range(len(x0)):
        if ii in phys:
            p1 = ax.plot(x0[ii],y0[ii],color='blue',marker='o',mfc='blue',ms=5)
            if ((ii not in posRad) and (ii not in tan)): 
                print fmt % (names[ii],mag[ii],r[ii],ar[ii]*1e3,sigmaR[ii],at[ii]*1e3,sigmaT[ii],cnt[ii],x0[ii],y0[ii])
        if ii in posRad:
            p2 = ax.plot(x0[ii],y0[ii],color='red',marker='o',mfc='None',mec='red',mew=1,ms=7)
        if ii in tan:
            p3 = ax.plot(x0[ii],y0[ii],color='green',marker='o',mfc='None',mec='green',mew=1,ms=11)
        if ((ii not in phys) and (ii not in posRad) and (ii not in tan)):
            p4 = ax.plot(x0[ii],y0[ii],color='black',marker='x',mfc='black',ms=5)
        if ((names[ii] == 'S1-51') | (names[ii] == 'S0-33') | (names[ii] == 'S1-92') | (names[ii] == 'S1-97')  | (names[ii] == 'S0-40')):
            ax.plot(x0[ii],y0[ii],color='m',marker='H')

    # Temporary
    #foo = names.index('S2-23')
    #ax.plot(x0[foo],y0[foo],'rs',ms=8)

    # ipshell()
    if markRef == True:
        p5 = ax.plot(inX,inY,'mx',ms=7,mew=1.5)

    # Overplot a circle at radius 3"
    #an = linspace(0,2*pi,100)
    #ax.plot(2.8*cos(an),2.8*sin(an),'k')
    #ax.plot([3,-2,-2,3,3],[-2,-2,2,2,-2],'k--')

    # Set up legend
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=10)

    if markRef == False:
        legend_items = ['Physical','Pos. Radial','Tangential','None']
        ax.legend((p1,p2,p3,p4),legend_items, numpoints=1, loc=3, prop=prop)
    else:
        legend_items = ['Physical','Pos. Radial','Tangential','None','Ref Star']
        ax.legend((p1,p2,p3,p4,p5),legend_items, numpoints=1, loc=3, prop=prop)
    py.plot([0],[0],'k+',ms=7,mew=2)
    lgdLines = ax.get_legend().get_lines()
    py.setp(lgdLines, visible=False)
    ax.set_xlabel('RA (arcsec)')
    ax.set_ylabel('Dec (arcsec)')
    ax.axis([np.floor(max(x0))+0.5,np.round(min(x0))-0.5,
             np.floor(min(y0))-0.5,np.round(max(y0))+0.5])
    fig.savefig(home + 'plots/accel_on_image%s.png' % suffix)

    py.clf()
    py.figure(figsize=(6,6))
    fig.subplots_adjust(left=0.1,right=0.95,top=0.95)
    py.plot(r, ar/are, 'r.', label='Radial')
    py.plot(r, at/ate, 'b.', label='Tangential')
    py.plot([0,7],[0,0],'k--')
    py.plot([0,7],[sigma,sigma],'k--')
    py.plot([0,7],[-sigma,-sigma],'k--')
    #py.plot([2.8,2.8],[-100,100],'k--')
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('Acceleration Significance (sigma)')
    #py.axis([0.0, 7.0, -30, 30])
    py.legend(numpoints=1,prop=prop,fancybox=True,loc=4)
    py.savefig(home + 'plots/accel_vs_radius%s.png' % suffix)


    #plot points on mag vs radius and show which ones have significant physical and non-physical accel
    py.clf()
    for ii in range(len(x0)):
        if ii in phys:
            p1 = py.plot(r[ii],mag[ii],color='blue',marker='o',mfc='blue',ms=5)
#            if ((ii not in posRad) and (ii not in tan)): 
#                print fmt % (names[ii],mag[ii],r[ii],ar[ii]*1e3,sigmaR[ii],at[ii]*1e3,sigmaT[ii],cnt[ii],x0[ii],y0[ii])
        if ii in posRad:
            p2 = py.plot(r[ii],mag[ii],color='red',marker='o',mfc='None',mec='red',mew=1,ms=7)
        if ii in tan:
            p3 = py.plot(r[ii],mag[ii],color='green',marker='o',mfc='None',mec='green',mew=1,ms=11)
        if ((ii not in phys) and (ii not in posRad) and (ii not in tan)):
            p4 = py.plot(r[ii],mag[ii],color='black',marker='o',mfc='black',ms=5)
    
    legend_items = ['Physical','Pos. Radial','Tangential','None']
 #   py.legend((p1,p2,p3,p4),legend_items, numpoints=1, loc=3, prop=prop)
 #   lgdLines = ax.get_legend().get_lines()
 #   py.setp(lgdLines, visible=False)
    py.xlabel('Radius (")')
    py.ylabel('Mag')
    py.axis([np.floor(min(r)),np.round(max(r)),
             np.round(max(mag)),np.floor(min(mag))])
    py.savefig(home + 'plots/mag_radius_phy_non%s.png' % suffix)

def histAccel(alnDir='11_10_26/',align='align/align_d_rms_1000_abs_t',
              poly='polyfit_c/fit', points='points_c/',
              starlist='yngstars', plotSigAcc=False, sigma=5.0,
              magCut=15, nEpochs=30):
    """
    Make a histogram of the accelerations in the radial/tangential,
    X/Y, and inline/perp. direction of motion. Also, the large outliers
    (> 3 sigma) are printed out.

    Inputs:
    root   = The root directory of an astrometry analysis
                (e.g. '08_02_16/' or './' if you are in the directory).
    align     = The align root file name (including the directory relative
                to root). Make sure that polyfit was run on this align
		output.
    poly      = The polyfit root file name (including the directory relative
                to root). This should be run on the same align as above.
    points    = The points directory.
    starlist  = Only plot specific subset of stars. Must be 'oldstars', 
                'yngstars', 'all'.
    plotSigAcc = Set to True to plot the fits for stars with significant 
                 accelerations.
    sigma      = Significance value for plotting (used with plotSigAcc)
    nEpochs    = Number of epochs a star must have been detected in for this analysis 

    Output:
    plots/polyfit_hist_accel.eps (and png)
    -- Contains the histograms of the accelerations.

    plots/polyfit_hist_accel_nepochs.eps (and png)
    -- Contains a plot of number of stars vs. acceleration significance.
    """
    outdir = home + 'plots/'
    
    s = loadPop(alnDir=alnDir,align=align,
                poly=poly,points=points,starlist=starlist)

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
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')

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
    foo = np.where((mag < 13.5) & (cnt > nEpochs))[0]
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

    foo = np.where((mag < 15.5) & (cnt > nEpochs))[0]
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
    sigmaAll = np.concatenate([sigmaR, sigmaT])
    (ksdR, kspR) = stats.stats.kstest(sigmaR, 'norm', N=len(sigmaR))
    (ksdT, kspT) = stats.stats.kstest(sigmaT, 'norm', N=len(sigmaT))
    (ksdA, kspA) = stats.stats.kstest(sigmaAll, 'norm', N=len(sigmaAll))
    #print ''
    #print 'KS Test for normality (prob that observed is gaussian):'
    #print '\tRadial Acc. KS Prob. = %5.3f' % (kspR)
    #print '\tTangen Acc. KS Prob. = %5.3f' % (kspT)
    #print '\tCombo Acc. KS Prob. = %5.3f' % (kspA)
    #print ''

    tag = ''
    if starlist=='all':
        tag = '_all'
    if starlist=='oldstars':
        tag = '_old'
    if starlist=='yngstars':
        tag = '_yng'

    py.clf()

    fig = py.figure(figsize=(8,8))
    fig.subplots_adjust(left=0.1,right=0.95,top=0.95,hspace=0.25)
    fig.clf()

    print ''
    print 'Number of stars plotted: %i' % len(ax)

    def makePlots(val, acc, num, label):
        ax1 = fig.add_subplot(3,2,num)
        ax1.hist(val, bins=range(-25, 25, 1), color='b')
        #ax1.hist(val, bins=range(-8, 8, 1), color='b')
        #ax1.axis([-8, 8, 0, (len(val)/3) + 2])
        ax1.axis([-25, 25, 0, (len(val)/3) + 2])
        ax1.set_xlabel(label)
        ax1.set_ylabel('N')

        print ''
        fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %2d epochs'
        hdr = '%15s  %5s  %6s  %17s  %14s  %8s '
        if label == 'Radial':
            idx = (np.where((val < -sigma)))[0]
            print '%s Significant %s (physical) Accelerations' % (len(idx), label)
            print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)', 'Nepochs')
            if len(idx) > 0:
                for i in idx:
                    print fmt  % (names[i], mag[i], r[i], acc[i]*1e3, val[i], cnt[i])

            # Get significant unphysical accelerations (positive radial)
            idx = (np.where((val > sigma)))[0]
            print
            print '%s Significant Positive %s (unphysical) Accelerations' % \
                  (len(idx), label)
            print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)',
                         'Nepochs')
            if len(idx) > 0:
                for i in idx:
                    print fmt % (names[i], mag[i], r[i], acc[i]*1e3, val[i], cnt[i])

        if (label == 'Tangential'):
            # Get significant unphysical accelerations (tangential)
            idx = (np.where((np.abs(val) > sigma)))[0]
            print
            print '%s Significant %s (unphysical) Accelerations' % (len(idx), label)
            print hdr % ('Name','K','r (")', 'a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs')
            if len(idx) > 0:
                for i in idx:
                    print fmt % (names[i], mag[i], r[i], acc[i]*1e3, val[i], cnt[i])


    # XY (in units of sigma)
    makePlots(ax / axe, ax, 1, 'X')
    makePlots(ay / aye, ay, 2, 'Y')
    
    # Radial/Tangential (in units of sigma)
    makePlots(ar / are, ar, 3, 'Radial')
    makePlots(at / ate, at, 4, 'Tangential')

    # In line of Motion and out (in units of sigma)
    makePlots(an / ane, an, 5, 'Perpendicular to v')
    makePlots(am / ame, am, 6, 'Parallel to v')

    fig.savefig(outdir + 'polyfit_hist_accel%s.png' % tag)


    # Plot up just the radial and tangential accels
    usetexTrue()
    fig = py.figure(figsize=(8,4))
    fig.subplots_adjust(hspace=0.2, left=0.1, right=0.95,
                       top=0.95, bottom=0.15)
    fig.clf()
    ax1 = fig.add_subplot(121)
    ax1.hist(ar/are, bins=np.arange(-25, 25, 1.0), histtype='step', color='k',lw=2)
    ax1.xaxis.set_minor_locator(MultipleLocator(5))
    ax1.axis([-30,30,0,23])
    #ax1.axis([-10,10,0, (len(ar)/4)+2])
    ax1.set_xlabel(r'{\bf Radial Acceleration ($\sigma$)}')
    ax1.set_ylabel(r'{\bf N}')
    ax2 = fig.add_subplot(122)
    ax2.hist(at/ate, bins=np.arange(-25, 25, 1.0), histtype='step', color='k',lw=2)
    ax2.xaxis.set_minor_locator(MultipleLocator(5))
    #ax2.axis([-10,10,0, (len(ar)/4)+2])
    ax2.axis([-30,30,0,23])
    ax2.set_xlabel(r'{\bf Tangential Acceleration ($\sigma$)}')
    fig.savefig(outdir + 'polyfit_hist_accel%s_radtan.png' % tag)
#    fig.savefig(outdir + 'polyfit_hist_accel%s_radtan.eps' % tag)
    usetexFalse()
    
    # Plot up accel histograms binned by radius
    fig = py.figure(figsize=(6,12))
    fig.subplots_adjust(hspace=0.3, wspace=0.2, left=0.1, right=0.95,
                       top=0.95, bottom=0.1)
    fig.clf()
    rlabel = ['r < 1"', 'r = 1-2"', 'r = 2-3"', 'r = 3-4"', 'r = 4-5"', ]
    for ii in range(5):
        rr = np.where((r >= (ii)) & (r < (ii+1)))[0]
        if len(rr) > 0:
            ax1 = fig.add_subplot(5,1,ii+1)
            ax1.hist(ar[rr]/are[rr], bins=range(-25, 25, 1), histtype='step', color='b')
            ax1.axis([-25,25,0,50])
            ax1.text(-15,5,rlabel[ii])
            ax1.set_ylabel('N')
            if ii == 0:
                ax1.set_title('Radial Accelerations (sigma) binned by Radius')
    fig.savefig(outdir + 'polyfit_hist_radial_accel%s_radius.png' % tag)

    # Temporary
    #foo = names.index('S0-6')
    #print 'S0-6:'
    #print 'a_rad = %6.3f +- %5.3f mas/yr/yr' % (ar[foo]*1.e3, are[foo]*1.e3)
    #print 'a_tan = %6.3f +- %5.3f mas/yr/yr' % (at[foo]*1.e3, ate[foo]*1.e3)
    # end temporary

    # Analyze the non-central arcsec sources for Nepochs threshold
    idx = (np.where(r > 0.8))[0]

    fig = py.figure(figsize=(6,6))
    fig.subplots_adjust(hspace=0.2, left=0.1, right=0.95,
                       top=0.95, bottom=0.1)
    fig.clf()
    ax1 = fig.add_subplot(211)
    ax1.plot(ar[idx] / are[idx], cnt[idx], 'k.')
    ax1.axis([-20, 20, 0, 45])
    ax1.set_xlabel('Radial Acc. Sig. (sigma)')
    ax1.set_ylabel('Number of Epochs Detected')

    ax2 = fig.add_subplot(212)
    ax2.plot(at[idx] / ate[idx], cnt[idx], 'k.')
    ax2.axis([-20, 20, 0, 45])
    ax2.set_xlabel('Tangent. Acc. Sig. (sigma)')
    ax2.set_ylabel('Number of Epochs Detected')

    fig.savefig(outdir + 'polyfit_hist_accel_nepochs%s.png' % tag)

    # Plot significant accelerations?
    if plotSigAcc == True:
        # Physical accelerations (negative radial)
        idx1 = np.where((ar/are) < -1.0*sigma)[0]

        if len(idx1) > 0:
            for i in range(len(idx1)):
                plotStar(names[i].strip(),alnDir=alnDir,align=align,
                         poly=poly,points=points,radial=True)

        # Unphysical accelerations (positive radial)
        idx2 = np.where((ar/are) > sigma)[0]

        if len(idx2) > 0:
            for i in range(len(idx2)):
                # If already plotted this, skip it
                if i in idx1:
                    continue
                else:
                    plotStar(names[i].strip(),alnDir=alnDir,align=align,
                             poly=poly,points=points,radial=True)

        # Unphysical accelerations (tangential)
        idx3 = np.where(np.abs(at/ate) > sigma)[0]

        if len(idx3) > 0:
            for i in range(len(idx3)):
                # If already plotted this, skip it
                if ((i in idx1) | (i in idx2)):
                    continue
                else:
                    plotStar(names[i].strip(),alnDir=alnDir,align=align,
                             poly=poly,points=points,radial=True)

    ##############################################################
    #SC added
    ##############################################################
    hdr = '%15s  %5s  %6s  %17s  %14s  %17s  %14s  %8s  %3s '
    fmt = '%15s  %5.2f   %5.3f"  %8.3f mas/yr^2  %8.2f sigma  %8.3f mas/yr^2  %8.2f sigma  %2d epochs  %3s'
    idsig = (np.where((sigmaR < -sigma) & (np.abs(sigmaT) < sigma)))[0] #sig accel
    idnp = (np.where((sigmaR > sigma) | (np.abs(sigmaT) > sigma)))[0] #non physical accel

    py.figure(6)
    py.clf()
    py.plot(ar*1e3,at*1e3,'k.')
    leg1=py.errorbar(ar[idnp]*1e3,at[idnp]*1e3,xerr=are[idnp]*3e3,yerr=ate[idnp]*3e3,fmt='.',label='Non Phys. Accel.')
    leg2=py.errorbar(ar[idsig]*1e3,at[idsig]*1e3,xerr=are[idsig]*3e3,yerr=ate[idsig]*3e3,fmt='.',label='Sig. Accel.')
    py.plot([-2,5],[0,0],'k')
    py.plot([0,0],[-1,1],'k')
    py.xlabel('Radial Acceleration (mas/yr/yr)')
    py.ylabel('Tangent Acceleration (mas/yr/yr)')
    py.legend()
    py.savefig(outdir + 'tan_rad.png')

    print
    print '%s SC Significant physical Accelerations' % (len(idsig))
    print hdr % ('Name','K','r (")', 'a_rad (mas/yr^2)','a_rad (sigma)','a_tan (mas/yr^2)','a_tan (sigma)', 'Nepochs', 'Fit')
    if len(idsig) > 0:
        for i in idsig:
            if pass_f[i] == 1:
                fit = 'Acc'
            elif pass_f[i] == 0:
                fit = 'Vel'
            print fmt % (names[i], mag[i], r[i], ar[i]*1e3, sigmaR[i], at[i]*1e3, sigmaT[i], cnt[i], fit)

    #####################################################################


        
def plotStar(starName,alnDir='13_08_21/',align='align/align_d_rms_1000_abs_t',
             poly='polyfit_c/fit', points='points_c/', radial=False, accel_fit=True):

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

def accelLimit(alnDir='11_07_05/',
              align = 'align/align_d_rms_1000_abs_t',
              poly='polyfit_c/fit', points='points_c/',
              starlist='yngstars', nEpochs=10):
    """
    Plot plane-of-the-sky acceleration limits as a function of projected distance.

    Inputs:
    root   = The root directory of an astrometry analysis
                (e.g. '08_02_16/' or './' if you are in the directory).
    align     = The align root file name (including the directory relative
                to root). Make sure that polyfit was run on this align
		output.
    poly      = The polyfit root file name (including the directory relative
                to root). This should be run on the same align as above.
    points    = The points directory.
    starlist  = Only plot specific subset of stars. Must be 'oldstars', 
                'yngstars', 'all'. 
    nEpochs   = Minimum number of epochs a star must be detected in (default=10).

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
    s = loadPop(alnDir=alnDir,align=align,
                poly=poly,points=points,starlist=starlist)

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
    t0 = s.getArray('t0x')

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

    # Determine # of data points per star
    pntcnt = np.zeros(len(names))
    for ii in range(len(names)):
        pntFileName = '%s%s%s%s.points' % (rootDir, alnDir, points, names[ii])
        pntFile = open(pntFileName)
        data = pntFile.readlines()
        pntcnt[ii] = len(data)

    # Only plot acceleration upper limits for stars detected in Nepochs:
    idx = np.where(pntcnt >= nEpochs)[0]

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
    r2d = r2d[idx]
    r2d_pc = r2d_pc[idx]
    rcgs = rcgs[idx]
    t0 = t0[idx]

    # Total acceleration
    atot = py.hypot(ax, ay)
    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

    # Construct an array of radii out to 4 arcsec in steps of 0.1''
    r = (arange(10 * 5) * 0.1) + 0.1
    r_au = r * dist
    r_pc = r_au / au_in_pc
    r_cm = r_au * cm_in_au

    # Determine the theoretical curve for a vs. r
    a_cm_s2 = -GM/ r_cm**2
    a_km_s_yr = a_cm_s2 * sec_in_yr / 1.0e5


    # What is the polyfit acceleration upper limit
    # along the radial component? Convert 2D cartesian to circular:
    (ar, at, are, ate) = util.xy2circErr(x0, y0, ax, ay,
                                         x0e, y0e, axe, aye)
    
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
    _out = open(rootDir + alnDir + 'tables/kinematics'+tag+'.txt','w')
    hdr = '#%10s  %7s  %12s  %12s  %12s  %7s  %7s  %15s   %9s\n'
    fmt = '%10s  %7.2f  %12.4e  %12.4e  %12.4e  %7.2f  %7.2f  %15.2f  %11.3f\n'
    _out.write(hdr % ('Name','t_ref','X (km)', 'Y (km)', 'Zmin (km)','Vx (km/s)',\
                      'Vy (km/s)','a_rad (km/s/yr)','Zmin (pc)'))

    print
    print 'Stars with (radial acc + 3 sigma) below curve:'
    print '%8s  %7s  %7s  %7s   %16s' % \
          ('Star','X (pc)', 'Y (pc)', 'Zmin (pc)', 'a_rad (km/s/yr)')
    print '--------------------------------------------------------'

    # Plot zmin vs. 2d radius
    # If stars fall below the curve in the above plot, then their 3D distance must
    # be greater than their projected distance.  Can therefore get z_min.
    min_cnt = 0
    for ii in range(len(rcgs)):
        # Upper limits to accelerations give lower limits to 3D radii
        # Figure out which stars are below the curve, use ar + 3 sigma to get sources (alimcgs)
        rmin2 = (-GM * rcgs[ii] / alimcgs[ii])**(2.0/3.0) # Gives (minimum r)-squared
        if ((rmin2 - rcgs[ii]**2) > 0): # Below the curve (at 3 sigma)
            # Get the indices of these stars for later
            sig_id = concatenate([sig_id,[ii]])

            zmincgs = np.sqrt(rmin2 - rcgs[ii]**2)
            r_p = rcgs / cm_in_pc
            zmin = zmincgs / (cm_in_pc) # pc
            zmin_km = zmincgs / 1.e5

            py.plot([r_p[ii]], [zmin], 'k.')
            py.xlabel('Projected Radius (pc)')
            py.ylabel('Minimum Z (pc)')

            print '%8s  %7.3f  %7.3f  %7.3f     %6.3f' % \
                  (names[ii], x0_pc[ii], y0_pc[ii], zmin, ar[ii])

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
        if (-ar[ii] - (3.*are[ii]) > 0.0): # and they're significantly different from zero
            det = py.errorbar(r2d_pc[ii], -ar[ii], yerr=(3.*are[ii]), fmt='b.',ms=8)
            xt = r2d_pc[ii]
            yt = -ar[ii]
            nm = names[int(ii)]
            #text(xt,yt,nm,fontsize=9)
            dtct = np.concatenate([dtct,[ii]])
    # Plot upper limits (from top of this program)
    for ii in range(len(r2d_pc)):
        if ii in dtct:
            continue
        else:
            uppr = py.plot([r2d_pc[ii]], [-accLimRadial[ii]], 'k_')
            py.arrow(r2d_pc[ii], -accLimRadial[ii], 0, -5, hold=True, color='black',\
                     width=0.0005, head_length=1, linewidth=0.002) #, \
            #          units='x',width=0.0005, color='black')
    print 'Number of stars below curve and significantly different from zero: %i' % len(dtct)
    py.plot(r_pc, -a_km_s_yr, '--k')
    lgd1 = 'Upper Limit'
    lgd2 = 'Detection'
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=10)
    lgd = py.legend((uppr[0], det[0]), (lgd1,lgd2), fancybox=True, numpoints=1,prop=prop)
    py.text(0.05, 80, r'{\bf $|$a$|_{max} = \frac{GM}{R^2}$}')
    py.xlabel(r'{\bf Projected Radius (pc)}')
    py.ylabel(r'{\bf $|$a$_R|$ (km/s/yr)}')
    py.title('Acceleration Detections \& Upper Limits')
    py.axis([0,0.2,0,100])
    py.savefig(outdir+'r2d_accel_detectUpLim%s.png' % tag)
    py.savefig(outdir+'eps/r2d_accel_detectUpLim%s.eps' % tag)
    py.close(2)

    print ''
    unphys = np.where(accLimRadial > 0.0)[0]
    print 'Number of positive radial accelerations: %i' % len(unphys)

    #fig=figure(3)
    #clf()
    #ax = p3.Axes3D(fig)
    #py.figure(4,figsize=(6,6))
    #py.clf()
    #for ii in range(len(rcgs)):
        #rmin2 = (-GM * rcgs[ii] / alimcgs1[ii])**(2.0/3.0) # Gives (minimum r)-squared
        #if ((rmin2 - rcgs[ii]**2) > 0):
            #if (-ar[ii] - (3.*are[ii]) > 0.0): # significant acceleration (above zero)
                #zmincgs = np.sqrt(rmin2 - rcgs[ii]**2)
                #zmin = zmincgs / (cm_in_pc)

                #figure(3)
                #ax.scatter3D([x0_pc[ii]],[zmin],[y0_pc[ii]],facecolor='blue')
                #if tag == '_old':
                #    ax.plot3D([x0_pc[ii],x0_pc[ii]],[zmin,zmin],[y0_pc[ii],-0.06],\
                #              color='g',ls='--',lw='1')
                #else:
                #    ax.plot3D([x0_pc[ii],x0_pc[ii]],[zmin,zmin],[y0_pc[ii],-0.08],\
                #              color='g',ls='--',lw='1')
                #rmin = np.sqrt(x0_pc[ii]**2+y0_pc[ii]**2+zmin**2)

                #figure(4)
                #plotVvsR(names[ii],rmin,vx_kms[ii],vy_kms[ii],vxe_kms[ii],vye_kms[ii])
                

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

    #py.figure(4)
    #py.xlabel('3D Radius (pc)')
    #py.ylabel('Proper Motion (km/s)')
    #py.title('Proper Motion vs. Radius')
    #py.savefig(outdir+'VvsR%s.png' % tag)

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


def table_young_accel(alnDir='11_04_24_oldmm/', align='align/align_d_rms_1000_abs_t',
              poly='polyfit_c/fit', points='points_c/'):
    """
    Creates a latex table of known young stars and their
    accelerations, significance, etc...
    """

    s = loadPop(alnDir=alnDir,align=align,
                poly=poly,points=points,starlist='yngstars')
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
    ax = s.getArray('ax') * 1.e3
    ay = s.getArray('ay') * 1.e3
    axe = s.getArray('axe') * 1.e3
    aye = s.getArray('aye') * 1.e3
    cnt = s.getArray('cnt')
    mag = s.getArray('mag')
    t0 = s.getArray('t0x') # t0x and t0y are the same

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

    phys = np.where((ar/are) < -3.0)[0]
    posRad = np.where((ar/are) > 3.0)[0]
    tan = np.where(np.abs(at/ate) > 3.0)[0]

    # Which stars have a significant acceleration (physical or not)
    _sig = np.concatenate([phys,posRad])
    _sig = np.concatenate([_sig,tan])
    _sig = np.unique(_sig)

    # Total acceleration
    atot = py.hypot(ax, ay)
    atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

    # Create the latex file
    out = open(rootDir + alnDir + 'tables/young_accel.tex','w')
    out.write('\\documentclass{aastex} \n')
    out.write('\\begin{singlespace} \n')
    out.write('\\begin{deluxetable}{lccccccccr} \n')
    #out.write('\\rotate \n')
    out.write('\\tabletypesize{\\scriptsize} \n')
    out.write('\\setlength{\\tabcolsep}{1.5mm} \n')
    out.write('\\tablewidth{0pt} \n')
    out.write('\\begin{document} \n')
    out.write('\\tablecaption{Young Star Accelerations}\n')
    out.write('\\tablehead{ \n')
    out.write('  \\colhead{Star Name} & \n')
    out.write('  \\colhead{$Kp$} & \n')
    out.write('  \\colhead{$T_o$} & \n')
    out.write('  \\colhead{$X_o$} & \n')
    out.write('  \\colhead{$Y_o$} & \n')
    out.write('  \\colhead{$a_r$} & \n')
    out.write('  \\colhead{$a_r$} & \n')
    out.write('  \\colhead{$a_t$} & \n')
    out.write('  \\colhead{$a_t$} & \n')
    out.write('  \\colhead{Signif-} \\\\ \n')
    out.write('%\n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(mag)} & \n')
    out.write('  \\colhead{(yr)} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{(mas/yr$^2$)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
    out.write('  \\colhead{(mas/yr$^2$)} & \n')
    out.write('  \\colhead{(sigma)} & \n')
    out.write('  \\colhead{icant?} \n')
    out.write('} \n')
    out.write('\\startdata \n')

    fmt = '%15s  %1s  %5.2f  %1s  %7.2f  %1s  %6.2f  %1s  %6.2f  %1s  %7.3f  %5s  '
    fmt += '%6.3f  %1s  %5.2f  %1s  %7.3f  %5s  %6.3f  %1s  %5.2f  %1s  %1s  %4s\n'

    for ii in range(len(x0)):
        if ii in _sig:
            sig = '*'
        else:
            sig = ''
        out.write(fmt % (names[ii], '&', mag[ii], '&', t0[ii], '&', x0[ii], '&',
                         y0[ii], '&', ar[ii], '$\pm$', are[ii], '&', np.abs((ar/are))[ii],
                         '&', at[ii], '$\pm$', ate[ii], '&', np.abs((at/ate))[ii], '&', sig, '\\\\'))

    out.write('\\\\\n')
    out.write('\\enddata \n')
    out.write('\\end{deluxetable} \n')
    out.write('\\end{singlespace} \n')
    out.write('\\end{document} \n')

    out.close()


def plotCentral(alnDir='11_07_05/', align='align/align_d_rms_1000_abs_t',
                poly='polyfit_c/fit', points='points_c/',starlist='all'):
    """
    Calls plotStar() function above on central arcsecond sources.
    """
    
    # Read in the stars and find the central arcsecond stars
    # Load data
    s = loadPop(alnDir=alnDir,align=align,
                poly=poly,points=points,starlist=starlist)

    names = s.getArray('name')

    for name in names:
        if ('S0' in name):
            print 'Plotting %s' % name
            plotStar(name,alnDir=alnDir,align=align,poly=poly,points=points)
        else:
            continue
        




def velocity_vs_accel(alnDir='11_10_26/',align='align/align_d_rms_1000_abs_t',
                      poly='polyfit_c/fit', points='points_c/',
                      pvalue=4.0, youngOnly=True, verbose=True,
                      dumpPickle=False, loadPickle=False,
                      extraPlots=False):
    """
    Computes the F statistic to determine whether velocity
    or acceleration fits are better for each star.
    Calls Tuan's accelclass code to do the F test.
    Input:
    pvalue -- Significance threshold for F test

    Returns names of stars with F statistic indicating acceleration fit
    is better than velocity fit, and the associated probabilities.
    """

    signif = scipy.special.erfc(pvalue/np.sqrt(2.0))
    print 'Significance threshold: %5.3f' % signif

    root = rootDir + alnDir

    if loadPickle == True:
        pFile = open(root + 'tables/ftest_data.pickle', 'r')
        data = pickle.load(pFile)
    else:
        import accel_class as acc
        data = acc.accelClass(rootDir=root,align=align,poly=poly,points=points,
                              verbose=verbose)
        data.computeFTest()
        data.computeFTestJerk()

    xFprob = data.xFProb # F test p-value for accel vs. vel
    yFprob = data.yFProb

    xFprobJ = data.xFProbJrk # F test p-value for jerk vs. accel
    yFprobJ = data.yFProbJrk

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

    # Plot up a histogram of the F probabilities
    step = 0.025
    binsIn = np.arange(0, 1.0+step, step)
    py.figure(1)
    py.figure(figsize=(6,6))
    py.clf()
    (nx,bx,ppx) = py.hist(xFprob,bins=binsIn,color='r',histtype='step',label='X',normed=True)
    (ny,by,ppy) = py.hist(yFprob,bins=binsIn,color='b',histtype='step',label='Y',normed=True)
    py.xlabel('F Probability',fontsize=12)
    py.ylabel('N',fontsize=12)
    py.legend(numpoints=1)
    py.savefig(outdir+'hist_Fprob_all.png')

    # Find young stars and determine whether velocity or accel
    # fit is better based on F probability
    _yng = young.youngStarNames()
    all_yng = []

    for yy in range(len(_yng)):
        yngstar = str(_yng[yy])
        yng = np.where(names == yngstar)[0]
        if len(yng) > 0:
            all_yng = np.concatenate([all_yng, yng])

    all_yng = [int(yy) for yy in all_yng]
    if youngOnly == True:
        names = [names[nn] for nn in all_yng]
        nEpochs = nEpochs[all_yng]
        mag = mag[all_yng]
        x = x[all_yng]
        y = y[all_yng]
        xe = xe[all_yng]
        ye = ye[all_yng]
        vx = vx[all_yng]
        vy = vy[all_yng]
        vxe = vxe[all_yng]
        vye = vye[all_yng]
        ax = ax[all_yng]
        ay = ay[all_yng]
        axe = axe[all_yng]
        aye = aye[all_yng]
        r2d = r2d[all_yng]
        xFprob = xFprob[all_yng]
        yFprob = yFprob[all_yng]
        xFprobJ = xFprobJ[all_yng]
        yFprobJ = yFprobJ[all_yng]

    # Trim out the central arcsecond stars
    idx = np.where(r2d > 0.0)[0]
    #idx = np.where(r2d > 0.8)[0]
    names = [names[nn] for nn in idx]
    nEpochs = nEpochs[idx]
    mag = mag[idx]
    x = x[idx]
    y = y[idx]
    xe = xe[idx]
    ye = ye[idx]
    vx = vx[idx]
    vy = vy[idx]
    vxe = vxe[idx]
    vye = vye[idx]
    ax = ax[idx]
    ay = ay[idx]
    axe = axe[idx]
    aye = aye[idx]
    r2d = r2d[idx]
    xFprob = xFprob[idx]
    yFprob = yFprob[idx]
    xFprobJ = xFprobJ[idx]
    yFprobJ = yFprobJ[idx]


    # Convert to radial and tangential accelerations (for plotting)
    ar = ((ax*x) + (ay*y)) / r2d
    at = ((ax*y) - (ay*x)) / r2d
    are =  (axe*x/r2d)**2 + (aye*y/r2d)**2
    are += (y*xe*at/r2d**2)**2 + (x*ye*at/r2d**2)**2
    are =  np.sqrt(are)
    ate =  (axe*y/r2d)**2 + (aye*x/r2d)**2
    ate += (y*xe*ar/r2d**2)**2 + (x*ye*ar/r2d**2)**2
    ate =  np.sqrt(ate)

    sigmaR = ar/are
    sigmaT = at/ate

    # How many epochs was each star detected in?
    cnt = np.zeros((len(names)),dtype=int)
    for nn in range(len(names)):
        pts = asciidata.open('%s%s%s.points' % (root, points, names[nn]))
        ep = pts[0].tonumpy()
        cnt[nn] = len(ep)

    if verbose == True:
        print '%i young stars' % len(r2d)
        print ''
    # Which young stars should we be using acceleration fits for?
    # See which ones have p-values from F test below the significance value chosen
    acc = np.where((xFprob < signif) | (yFprob < signif))[0]
    accX = np.where(xFprob < signif)[0]
    accY = np.where(yFprob < signif)[0] 
    # What about jerk fits?
    jrk = np.where((xFprobJ < signif) | (yFprobJ < signif))[0]
    jrkX = np.where(xFprobJ < signif)[0]
    jrkY = np.where(yFprobJ < signif)[0] 

    #if ((len(accX) > 0) | (len(accY) > 0)):
    #    print 'Young stars with p < %7.2e (%3.1f sigma) that should have acceleration fits (N=%i):' % \
    #          (signif, pvalue, len(accX)+len(accY))
    #    hdr = '%12s   %2s   %7s   %7s   %5s   %5s   %5s %6s %6s  %5s %6s %6s'
    #    fmt = '%12s   %2i   %7.2e   %7.2e   %5.2f   %5.2f   %6.3f %5.3f (%6.2f)   %6.3f %5.3f (%6.2f)'
    #    print ''
    #    print 'Passed F test in X direction:'
    #    print hdr % ('Name', 'N', 'X_Fprob', 'Y_Fprob', 'ax', 'ay', 'ar','are', '(sig)', 'at', 'ate', '(sig)')
    #    for aa in accX:
    #        print fmt % (names[aa], cnt[aa], xFprob[aa], yFprob[aa],
    #                     ax[aa]*1.e3, ay[aa]*1.e3, ar[aa]*1.e3, are[aa]*1.e3,ar[aa]/are[aa],
    #                     at[aa]*1.e3, ate[aa]*1.e3, at[aa]/ate[aa])
#
#        print ''
#        print 'Passed F test in Y direction:'
#        print hdr % ('Name', 'N', 'X_Fprob', 'Y_Fprob', 'ax', 'ay', 'ar', 'are', '(sig)', 'at', 'ate', '(sig)')
#        for aa in accY:
#            print fmt % (names[aa], cnt[aa], xFprob[aa], yFprob[aa],
#                         ax[aa]*1.e3, ay[aa]*1.e3, ar[aa]*1.e3, are[aa]*1.e3, ar[aa]/are[aa],
#                         at[aa]*1.e3, ate[aa]*1.e3, at[aa]/ate[aa])
#    print
#    print '*******************'
#    if ((len(accX) > 0) | (len(accY) > 0)):
#        print 'Young stars with p < %7.2e (%3.1f sigma) that should have jerk fits (N=%i):' % \
#              (signif, pvalue, len(jrkX)+len(jrkY))
#        #hdr = '%12s   %2s   %7s   %7s   %5s   %5s   %5s %6s  %5s %6s'
#        #fmt = '%12s   %2i   %7.2e   %7.2e   %5.2f   %5.2f   %6.3f (%5.3f)   %6.3f (%5.3f)'
#        print ''
#        print 'Passed F test in X direction:'
#        print hdr % ('Name', 'N', 'X_Fprob', 'Y_Fprob', 'ax', 'ay', 'ar','are', '(sig)', 'at', 'ate', '(sig)')
#        for aa in jrkX:
#            print fmt % (names[aa], cnt[aa], xFprob[aa], yFprob[aa],
#                         ax[aa]*1.e3, ay[aa]*1.e3, ar[aa]*1.e3, are[aa]*1.e3, ar[aa]/are[aa],
#                         at[aa]*1.e3, ate[aa]*1.e3, at[aa]/ate[aa])
#
#        print ''
#        print 'Passed F test in Y direction:'
#        print hdr % ('Name', 'N', 'X_Fprob', 'Y_Fprob', 'ax', 'ay', 'ar','are', '(sig)', 'at', 'ate', '(sig)')
#        for aa in jrkY:
#            print fmt % (names[aa], cnt[aa], xFprob[aa], yFprob[aa],
#                         ax[aa]*1.e3, ay[aa]*1.e3, ar[aa]*1.e3, are[aa]*1.e3, ar[aa]/are[aa],
#                         at[aa]*1.e3, ate[aa]*1.e3, at[aa]/ate[aa])
        

    py.figure(2)
    py.clf()
    (nx,bx,ppx) = py.hist(xFprob,bins=binsIn,color='r',histtype='step',label='X',normed=True)
    (ny,by,ppy) = py.hist(yFprob,bins=binsIn,color='b',histtype='step',label='Y',normed=True)
    py.xlabel('F Probability',fontsize=12)
    py.ylabel('N',fontsize=12)
    py.legend(numpoints=1)
    py.savefig(outdir+'hist_Fprob_young.png')

    # Plot accelerations (x vs. y)
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=10)
    py.figure(3)
    py.clf()
    py.figure(figsize=(10,5))
    py.subplots_adjust(left=0.1,right=0.95,top=0.9,bottom=0.1,hspace=0.3,wspace=0.3)
    py.subplot(1,2,1)
    py.plot(ax*1.e3,ay*1.e3,'k.')
    # Mark the stars that passed the F test
    if len(accX) > 0:
        for aa in accX:
            p1 = py.plot(ax[aa]*1.e3,ay[aa]*1.e3,'ro',mfc='None',mec='r')
    if len(accY) > 0:
        for aa in accY:
            p2 = py.plot(ax[aa]*1.e3,ay[aa]*1.e3,'ro',mfc='None',mec='b')
    legend_items = ['Pass F (X)','Pass F (Y)']
    py.legend((p1,p2),legend_items, numpoints=1, loc=3, prop=prop)
    py.xlabel('X Acceleration (mas/yr^2)')
    py.ylabel('Y Acceleration (mas/yr^2)')
    py.subplot(1,2,2)
    py.plot(ax*1.e3,ay*1.e3,'k.')
    # Mark the stars that passed the F test
    if len(accX) > 0:
        for aa in accX:
            py.plot(ax[aa]*1.e3,ay[aa]*1.e3,'ro',mfc='None',mec='r')
    if len(accY) > 0:
        for aa in accY:
            py.plot(ax[aa]*1.e3,ay[aa]*1.e3,'ro',mfc='None',mec='b')
    # Zoom in a little
    py.axis([-0.25,0.25,-0.25,0.25])
    py.title('Zoomed In')
    py.xlabel('X Acceleration (mas/yr^2)')
    py.ylabel('Y Acceleration (mas/yr^2)')
    py.savefig(outdir+'xyAccel_ftest.png')

    accnames = [names[nn] for nn in acc]
    accXFp = [xFprob[nn] for nn in acc]
    accYFp = [yFprob[nn] for nn in acc]

    if dumpPickle == True:
        pFile = open(root + 'tables/ftest_data.pickle', 'w')
        pickle.dump(data, pFile)

    if extraPlots == True:
        py.figure(4)
        py.clf()
        py.figure(figsize=(10,5))
        py.subplots_adjust(left=0.1,right=0.95,top=0.9,bottom=0.1,hspace=0.3,wspace=0.3)
        py.subplot(1,2,1)
        py.plot(ar*1.e3,at*1.e3,'k.')
        # Mark the stars that passed the F test in either X or Y
        if len(acc) > 0:
            for aa in acc:
                p1 = py.plot(ar[aa]*1.e3,at[aa]*1.e3,'ro',mfc='None',mec='r')
        #py.legend((p1,p2),legend_items, numpoints=1, loc=3, prop=prop)
        py.xlabel('Radial Acceleration (mas/yr^2)')
        py.ylabel('Tangential Acceleration (mas/yr^2)')
        py.plot(0.4,-1.0,'ro',mfc='None', mec='r')
        py.text(0.425,-1.02,'Pass F in',fontsize=10)
        py.text(0.425,-1.1,'either X or Y',fontsize=10)

        py.subplot(1,2,2)
        py.plot(ar*1.e3,at*1.e3,'k.')
        # Mark the stars that passed the F test
        if len(acc) > 0:
            for aa in acc:
                py.plot(ar[aa]*1.e3,at[aa]*1.e3,'ro',mfc='None',mec='r')
        # Zoom in a little
        py.axis([-0.25,0.25,-0.25,0.25])
        py.title('Zoomed In')
        py.xlabel('Radial Acceleration (mas/yr^2)')
        py.ylabel('Tangential Acceleration (mas/yr^2)')
        py.savefig(root+'plots/rad_tan_Accel_ftest.png')
        
        # Plot up the significance of the accelerations vs. F probability
        py.figure(5)
        py.clf()
        py.figure(figsize=(10,5))
        py.subplots_adjust(left=0.1,right=0.95,top=0.9,bottom=0.1,hspace=0.3,wspace=0.3)
        py.subplot(1,2,1)
        py.plot(xFprob*1.e5, ax/axe, 'r.')
        for nn in range(len(names)):
            if names[nn] == 'S1-14':
                py.text(xFprob[nn]*1.e5+0.05, ax[nn]/axe[nn]-1, names[nn], fontsize=8)
            else:
                py.text(xFprob[nn]*1.e5+0.05, ax[nn]/axe[nn], names[nn], fontsize=8)
        py.plot([-0.5, (xFprob[nn]*1.e5).max()], [4.0, 4.0], 'k--')
        py.plot([-0.5, (xFprob[nn]*1.e5).max()], [-4.0, -4.0], 'k--')
        # Just show the ones that passed the F test
        py.axis([-0.5,signif*1.e5,(ax/axe).min()-0.5,(ax/axe).max()+0.5])
        py.xlabel('F probability (X) * 1e5')
        py.ylabel('X Acceleration Significance (mas/yr^2)')
        py.title('Stars passing F test for accelerations')

        py.subplot(1,2,2)
        py.plot(yFprob*1.e5, ay/aye, 'b.')
        for nn in range(len(names)):
            py.text(yFprob[nn]*1.e5+0.05, ay[nn]/aye[nn], names[nn], fontsize=8)
        py.plot([-0.5, (yFprob[nn]*1.e5).max()], [4.0, 4.0], 'k--')
        py.plot([-0.5, (yFprob[nn]*1.e5).max()], [-4.0, -4.0], 'k--')
        # Just show the ones that passed the F test
        py.axis([-0.5,signif*1.e5,(ay/aye).min()-0.5,(ay/aye).max()+0.5])
        py.xlabel('F probability (Y) * 1e5')
        py.ylabel('Y Acceleration Significance (mas/yr^2)')
        py.title('Stars passing F test for accelerations')
        py.savefig(outdir+'xyAccel_vs_Fprob.png')

    return accnames, accXFp, accYFp


########################################################################
####SC added######################
def starInfo(starName, alnDir='13_08_21/',
              align = 'align/align_d_rms_1000_abs_t',
              poly='polyfit_c/fit', points='points_c/',
              starlist='all', nEpochs=30):
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

    print 'Projected radius of '+ starName + ': ' + r2d_pc[ii] + ' pc'
    print 'Radial acceleration of ' + starName + ': ' + ar[ii] + ' km/s/yr'
    print '3 sig error: ' + 3.*are[ii] + ' km/s/yr'
