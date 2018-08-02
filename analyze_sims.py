import asciidata, pyfits
import os, math, random, shutil
from decimal import Decimal
from gcwork import starset
from pysqlite2 import dbapi2 as sqlite
import pylab as py
import numpy as np
import matplotlib.nxutils as nx
from scipy import optimize
import nmpfit_sy
import pdb

#root = '/u/leo/TMTsims/iris_code_20120130/code/' 
root = '/u/syelda/research/gc/TMT/TMTsims/iris_code_20120130/code/'

def usetexTrue():
    py.rc('text', usetex=True)
    py.rc('font', **{'family':'sans-serif', 'size':16})
    #py.rc('axes', titlesize=20, labelsize=20)
    py.rc('xtick', labelsize=16)
    py.rc('ytick', labelsize=16)

def usetexFalse():
    py.rc('text', usetex=False)
    py.rc('font', family='sans-serif', size=14)
    #py.rc('axes', titlesize=16, labelsize=16)
    py.rc('xtick', labelsize=14)
    py.rc('ytick', labelsize=14)


def go(num=10, alnDir='align_18092012/', runSTF=False, runCAL=False,
       runALN=False, mainCorr=0.8, subCorr=0.6):
    """
    Run astrometric analysis on TMT simulations of GC images.

    Input:
    num (int):	      Number of individual frames.
    alnDir (str):     Align directory.
    runSTF (bool):    Run starfinder (def=False).
    runCAL (bool):    Run calibrate on starfinder star list (def=False).
    runALN (bool):    Run align and align_rms to get astrometric
    		      errors (def=False).
    mainCorr (float): Starfinder correlation threshold for main maps (def=0.8).
    subCorr (float):  Starfinder correlation threshold for submaps (def=0.6).
    """
    rootdir = root + alnDir
    stfdir = rootdir + 'starfinder/'
    alndir = stfdir + 'align/'

    dirStart = os.getcwd()

    mainmap = '%sdeepmap.fits' % rootdir
    submap1 = '%ssubmap_1.fits' % rootdir
    submap2 = '%ssubmap_2.fits' % rootdir
    submap3 = '%ssubmap_3.fits' % rootdir

    if runSTF == True:
        print 'COMBO starfinder (main + submaps)'
        os.chdir(stfdir)
        
        # Write an IDL batch file
        fileIDLbatch = 'idlbatch_combo_tmt'
        fileIDLlog = fileIDLbatch + '.log'
        _batch = open(fileIDLbatch, 'w')
        _batch.write("find_stf_tmt, ")
        _batch.write("'" + mainmap + "', ")
        _batch.write("'%3.1f', " % mainCorr)
        _batch.write("/makePsf, /makeRes, /makeStars")
        _batch.write("\n")
        _batch.write("find_stf_tmt, ")
        _batch.write("'" + submap1 + "', ")
        _batch.write("'%3.1f', " % subCorr)
        _batch.write("/makePsf, /makeRes, /makeStars")
        _batch.write("\n")
        _batch.write("find_stf_tmt, ")
        _batch.write("'" + submap2 + "', ")
        _batch.write("'%3.1f', " % subCorr)
        _batch.write("/makePsf, /makeRes, /makeStars")
        _batch.write("\n")
        _batch.write("find_stf_tmt, ")
        _batch.write("'" + submap3 + "', ")
        _batch.write("'%3.1f', " % subCorr)
        _batch.write("/makePsf, /makeRes, /makeStars")
        _batch.write("\n")
        _batch.write("exit\n")
        _batch.close()
        cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog            
        os.system(cmd)

        os.chdir(dirStart)

    if runCAL == True:
        ########
        # Run calibrate on the main and sub maps
        ########
        print '*** calibrating '
        os.chdir(stfdir)

        calSrcs = '-S 16C,16NW,16CC'
        cmd = 'calibrate_new %s -f 1 -M %i -c 12 -R -T %.1f  ' % \
                  (calSrcs, 6, 0.0)

        # Main Map
        main = 'deepmap_0.8_stf.lis' #% stfdir
        os.system(cmd + main)
        
        # Sub Map
        for s in range(1, 4):
                sub = 'submap_%d_0.8_stf.lis' %  (s)

                os.system(cmd + sub)

        os.chdir(dirStart)

    if runALN == True:
        ########
        # Run align and align_rms on star lists.
        ########
    
        print '*** align_rms'
        os.chdir(alndir)
    
        # Make align.list
        _alignlist = open('align.list', 'w')
    
        # Main -- add to list
        main = '../deepmap_0.8_stf_cal.lis'
        _alignlist.write('%s 30\n' % (main))
    
        # Subs -- add to list
        for s in range(1, 4):
            sub = '../submap_%d_%s_stf_cal.lis' % (s, 0.8)
            _alignlist.write('%s 30\n' % (sub))
    
        _alignlist.close()
    
        # Call java align
        cmd = 'java -Xmx1024m align -a 0 -R 5 -p -v '
        cmd += 'align.list'
        os.system(cmd)
    
        # Call align_rms
        os.system('align_rms -m align %d %d' % (3, 3))
    
        # Copy over output to new starlist BUT modify errors first.
        stars = asciidata.open('align_rms.out')
    
        # Do some cleanup of formatting
        stars[0].reformat('%-13s') # name
        stars[1].reformat('%6.3f') # mag
        stars[2].reformat('%8.3f') # epoch
        stars[3].reformat('%11.5f') # x
        stars[4].reformat('%11.5f') # y
        stars[5].reformat('%8.5f') # xerr
        stars[6].reformat('%8.5f') # yerr
        stars[7].reformat('%7.2f') # snr
        stars[8].reformat('%5.2f') # corr
        stars[9].reformat('%5d')   # numFrames
        stars[10].reformat('%8.3f') # ??
    
    
        stars.writeto('magTMT_rms.lis')
    
        os.chdir(dirStart)

def align_rms(alnDir, eom=True, minSub=3, totSub=3, magCutOff=17.0):
    """
    Computes the rms error or error on the mean from an
    aligned set of main+sub maps. Similar to align_rms.c,
    but just does a quick and dirty calculation and plots
    up the results against K mag.
    """

    # Read in the align
    wdir = root + alnDir
    align = wdir + '/starfinder/align/align'
    plotdir = wdir + 'starfinder/align/'

    # Loads position data from align file
    s = starset.StarSet(align)

    starCnt = len(s.stars)
    epochCnt = len(s.stars[0].e)
    
    names = s.getArray('name')
    mag = s.getArray('mag')

    x = np.zeros((totSub+1, len(names)), float)
    y = np.zeros((totSub+1, len(names)), float)
    for ii in range(totSub+1):
        x[ii,:] = s.getArrayFromEpoch(ii,'xpix')
        y[ii,:] = s.getArrayFromEpoch(ii,'ypix')
        

    xerr = []
    yerr = []
    mags = []
    # Find the stars in minSub submaps
    for ss in range(len(names)):
        x_a = []
        y_a = []
        for mm in range(totSub+1):
            if x[mm,ss] < 0:
                continue
            else:
                # Place this star's position in an array
                x_a = np.concatenate([x_a, [x[mm,ss]]])
                y_a = np.concatenate([y_a, [y[mm,ss]]])

        if (len(x_a) == minSub + 1): # main map + submaps 
            x_a_mean = x_a[1:minSub+1].mean()
            y_a_mean = y_a[1:minSub+1].mean()
            
            # Save this star's magnitude as well
            mags = np.concatenate([mags,[mag[ss]]])

            xdev = 0
            ydev = 0
            for ii in range(1,minSub+1):
                xdev += (x_a[ii] - x_a_mean)**2
                ydev += (y_a[ii] - y_a_mean)**2
            xstd = np.sqrt(xdev/(minSub-1))
            ystd = np.sqrt(ydev/(minSub-1))

            if eom == True:
                xstd /= np.sqrt(minSub)
                ystd /= np.sqrt(minSub)

            xerr = np.concatenate([xerr,[xstd]])
            yerr = np.concatenate([yerr,[ystd]])

    magStep = 1.0
    magBins = np.arange(10.0, 23.0+magStep, magStep)
    
    errMag = np.zeros(len(magBins), float)

    err = (xerr + yerr) / 2.0

    ##########
    # Compute errors in magnitude bins
    ########## 
    #print '%4s  %s' % ('Mag', 'Err (mas)')
    for jj in range(len(magBins)):
        mMin = magBins[jj] - (magStep / 2.0)
        mMax = magBins[jj] + (magStep / 2.0)
        idx = (np.where((mags >= mMin) & (mags < mMax)))[0]

        if (len(idx) > 0):
            errMag[jj] = np.median(err[idx])
        

    usetexTrue()
    py.clf()
    py.subplots_adjust(left=0.15, top=0.92, right=0.95, bottom=0.1)
    py.figure(figsize=(6,6))
    py.semilogy(mags,xerr,'r.')
    py.semilogy(mags,yerr,'b.')
    py.semilogy(magBins, errMag, 'g.-', lw=2)
    py.xlabel('K')
    py.ylabel('Positional Uncertainty (pix)')
    py.axis([6,24,0.0001,2.0])
    py.savefig(plotdir+'posErr_vs_Kmag_pix.png')
    py.close()

    py.clf()
    py.subplots_adjust(left=0.15, top=0.92, right=0.95, bottom=0.1)
    py.figure(figsize=(6,6))
    py.semilogy(mags,xerr*4.,'r.')
    py.semilogy(mags,yerr*4.,'b.')
    py.semilogy(magBins, errMag*4., 'g.-', lw=2)
    py.xlabel('K')
    py.ylabel('Positional Uncertainty (mas)')
    py.axis([6,24,0.0001,8.0])
    py.savefig(plotdir+'posErr_vs_Kmag_mas.png')
    py.close()
    usetexFalse()

    ##########
    # 
    # Save relevant numbers to an output file.
    #
    ##########
    # Print out some summary information
    brt = np.where(magBins < magCutOff)[0]

    print 'Number of detections: %4d' % len(mags)
    print 'Average, Median Pos Error (uas) for K < %2i:  %7.5f, %7.5f uas' % \
          (magCutOff, errMag[brt].mean()*4000., np.median(errMag[brt])*4000.)



def plotPosError(starlist, raw=True, suffix='', radius=20, magCutOff=17.0,
                 title=True):
    """
    Make three standard figures that show the data quality 
    from a *_rms.lis file; errors determined using 3 submaps.

    1. astrometric error as a function of magnitude.
    2. photometric error as a function of magnitude.
    3. histogram of number of stars vs. magnitude.

    Use raw=True to plot the individual stars in plots 1 and 2.
    """
    # Load up the starlist
    lis = asciidata.open(starlist)

    # Plate scale
    scale = 0.004 # arcsec
    
    name = lis[0]._data
    mag = lis[1].tonumpy()
    x = lis[3].tonumpy()
    y = lis[4].tonumpy()
    xerr = lis[5].tonumpy()
    yerr = lis[6].tonumpy()
    snr = lis[7].tonumpy()
    corr = lis[8].tonumpy()

    merr = 1.086 / snr

    # Convert into arsec offset from field center
    # We determine the field center by assuming that stars
    # are detected all the way out the edge.
    xhalf = x.max() / 2.0
    yhalf = y.max() / 2.0
    x = (x - xhalf) * scale
    y = (y - yhalf) * scale
    xerr *= scale * 1.e3
    yerr *= scale * 1.e3

    r = np.hypot(x, y)
    #err = np.hypot(xerr, yerr)
    err = (xerr + yerr) / 2.0

    magStep = 1.0
    radStep = 1.0
    magBins = np.arange(10.0, 21.0, magStep)
    radBins = np.arange(0.5, 9.5, radStep)
    
    errMag = np.zeros(len(magBins), float)
    errRad = np.zeros(len(radBins), float)
    merrMag = np.zeros(len(magBins), float)
    merrRad = np.zeros(len(radBins), float)

    ##########
    # Compute errors in magnitude bins
    ########## 
    #print '%4s  %s' % ('Mag', 'Err (mas)')
    for mm in range(len(magBins)):
        mMin = magBins[mm] - (magStep / 2.0)
        mMax = magBins[mm] + (magStep / 2.0)
        idx = (np.where((mag >= mMin) & (mag < mMax) & (r < radius)))[0]

        if (len(idx) > 0):
            errMag[mm] = np.median(err[idx])
            merrMag[mm] = np.median(merr[idx])
        
        #print '%4.1f  %5.2f' % (magBins[mm], errMag[mm])
        
                       
    ##########
    # Compute errors in radius bins
    ########## 
    for rr in range(len(radBins)):
        rMin = radBins[rr] - (radStep / 2.0)
        rMax = radBins[rr] + (radStep / 2.0)
        idx = (np.where((r >= rMin) & (r < rMax) & (mag < magCutOff)))[0]

        if (len(idx) > 0):
            errRad[rr] = np.median(err[idx])
            merrRad[rr] = np.median(err[idx])

    idx = (np.where((mag < magCutOff) & (r < radius)))[0]
    errMedian = np.median(err[idx])

    ##########
    #
    # Plot astrometry errors
    #
    ##########
 
    # Remove figures if they exist -- have to do this
    # b/c sometimes the file won't be overwritten and
    # the program crashes saying 'Permission denied'
    if os.path.exists('plotPosError%s.png' % suffix):
        os.remove('plotPosError%s.png' % suffix)
    if os.path.exists('plotMagError%s.png' % suffix):
        os.remove('plotMagError%s.png' % suffix)
    if os.path.exists('plotNumStars%s.png' % suffix):
        os.remove('plotNumStars%s.png' % suffix)

    py.figure(figsize=(6,6))
    py.clf()
    py.subplots_adjust(left=0.15, top=0.92, right=0.95)
    if (raw == True):
        idx = (np.where(r < radius))[0]
        py.semilogy(mag[idx], err[idx], 'k.')
        
    py.semilogy(magBins, errMag, 'g.-', lw=2)
    #py.axis([8, 20, 1e-3, 1.0])
    py.axis([8, 22, 1e-3, 30.0])
    py.xlabel('K Magnitude')
    py.ylabel('Positional Uncertainty (mas)')
    if title == True:
        py.title(starlist)
    
    py.savefig('plotPosError%s.png' % suffix)
    py.savefig('plotPosError%s.eps' % suffix)

    ##########
    #
    # Plot photometry errors
    #
    ##########
    py.clf()
    if (raw == True):
        idx = (np.where(r < radius))[0]
        py.plot(mag[idx], merr[idx], 'k.')
        
    py.plot(magBins, merrMag, 'g.-')
    py.axis([8, 22, 0, 0.15])
    py.xlabel('K Magnitude')
    py.ylabel('Photo. Uncertainty (mag)')
    py.title(starlist)
    
    py.savefig('plotMagError%s.png' % suffix)

    ##########
    # 
    # Plot histogram of number of stars detected
    #
    ##########
    py.clf()
    idx = (np.where(r < radius))[0]
    (n, bb, pp) = py.hist(mag[idx], bins=np.arange(9, 22, 0.5), histtype='step',
                          color='k')
    py.xlabel('K Magnitude')
    py.ylabel('Number of Stars')

    py.savefig('plotNumStars%s.png' % suffix)

    # Find the peak of the distribution
    maxHist = n.argmax()
    maxBin = bb[maxHist]

    ##########
    # 
    # Save relevant numbers to an output file.
    #
    ##########
    # Print out some summary information
    brt = np.where(mag < 17.0)[0]
    print
    print 'Average Pos Error for K < 17.0: %6.4f mas' % err[brt].mean()
    print
    print 'Number of detections: %4d' % len(mag)
    print 'Median Pos Error (mas) for K < %2i, r < %4.1f (N=%i):  %7.5f mas' % \
          (magCutOff, radius, len(brt), errMedian)
    print 'Median Mag Error (mag) for K < %2i, r < %4.1f (N=%i):  %5.3f mag' % \
          (magCutOff, radius, len(brt), np.median(merr[idx]))
    #print 'Turnover mag = %4.1f' % (maxBin)


    out = open('plotPosError%s.txt' % suffix, 'w')
    out.write('Number of detections: %4d\n' % len(mag))
    out.write('Median Pos Error (mas) for K < %2i, r < %4.1f (N=%i):  %7.5f\n' % \
          (magCutOff, radius, len(brt), errMedian))
    out.write('Median Mag Error (mag) for K < %2i, r < %4.1f (N=%i):  %5.3f\n' % \
          (magCutOff, radius, len(brt), np.median(merr[idx])))
    #out.write('Turnover mag = %4.1f\n' % (maxBin))
    out.close()
    

def compare_in_vs_out(alnDir='align_18092012/',single=False,makeInFile=False,
                      removeEdges=False,suffix='',gridTest=False,getHiPrec=False):
    """
    Compare input positions with output positions extracted by
    starfinder.  Requires that an align was run, where the two
    starlists (in vs. out) were aligned.

    If makeInFile == True, then only the input list is made in a
    format consistent with starfinder star lists.
    If makeInFile == False, it is assumed that the input and
    output star lists have been aligned already.

    If gridTest == True, positions are assumed to be of a grid, and not the GC.
    If getHiPrec == True, selects only stars with input positions known to
    	10uas level. This requires a certain input list in the alnDir.
    """
    rootdir = root + alnDir
    if gridTest == True:
        outdir = rootdir 
    else:
        #outdir = rootdir + 'starfinder/in_vs_out/'
        #outdir = rootdir + 'starfinder_knownPSF/'
        outdir = rootdir + 'starfinder_knownPSF/in_vs_out/'

    if makeInFile == True:
        # Get the input positions
        dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'

        scale = 0.004 # arcsec per pixel
        size = [4096, 4096]
        radius = size[0]/2*scale*np.sqrt(2.0)  # actual simulated radius

        # Create a connection to the database file
        connection = sqlite.connect(dbfile)

        # Create a cursor object
        cur = connection.cursor()
        cur.execute('SELECT name,kp,x,y,r2d,ddate from stars where r2d < ?', [radius])

        name = []
        kp = []
        x = []
        y = []
        r2d = []
        ddate = []
        fmt = '%12s  %6.3f  %8.3f  %10.5f  %10.5f  %5.2f  %4.2f  %1i  %8.3f\n'
        for row in cur:
            name = np.concatenate([name, [row[0]]])
            kp = np.concatenate([kp, [row[1]]])
            x = np.concatenate([x, [row[2]]]) # arcsec
            y = np.concatenate([y, [row[3]]]) # arcsec
            r2d = np.concatenate([r2d, [row[4]]]) # arcsec
            ddate = np.concatenate([ddate, [row[5]]])
       
        x = -x/scale+size[0]/2.0
        y = y/scale+size[1]/2.0

        # Write pixel positions to a file
        _infile = rootdir+'input_gc_positions_database.lis'
        #_infile = rootdir+'input_gc_positions.lis'
        infile = open(_infile,'w')
        dum = 1.0
        # Include Sgr A* at the top
        for ii in range(len(name)):
            if ii == 0:
                infile.write(fmt % ('SgrA', 16.7, 2011.0, 2048.0, 2048.0,
                                    dum, dum, dum, dum))
                infile.write(fmt % (name[ii],kp[ii],ddate[ii],x[ii],y[ii],
                                    dum, dum, dum, dum))
            else:
                infile.write(fmt % (name[ii],kp[ii],ddate[ii],x[ii],y[ii],
                                    dum, dum, dum, dum))
        infile.close()
    

    # If we had to make the input file, then we probably can't run the next steps
    # until we align the starfinder star list with the input list
    if makeInFile == False:
        # plot the input vs. output aligned positions (must have run align already)
        if single == True:
            adir = 'single_maps'
            #asuffix = '1'
            asuffix = '2_a3'
            #adir = 'test4'
            #adir = ''
            #asuffix = '_probeArm_distorted'
            #asuffix = '_pa_to_orig'
            #asuffix = ''
            #asuffix = '_image'
        else:
            adir = 'deepmap'
            asuffix = ''
        #aln = asciidata.open(rootdir + 'starfinder/in_vs_out/%s/align%s.pos' % (adir,asuffix))
        aln = asciidata.open(outdir + '%s/align%s.pos' % (adir,asuffix))
        star = aln[0].tonumpy()
        inX = aln[1].tonumpy()
        inY = aln[2].tonumpy()
        outX = aln[3].tonumpy()
        outY = aln[4].tonumpy()
    
        #_cnt = asciidata.open(rootdir + 'starfinder/in_vs_out/%s/align%s.vel' % (adir,asuffix))
        _cnt = asciidata.open(outdir + '%s/align%s.vel' % (adir,asuffix))
        cnt = _cnt[1].tonumpy()
        mag = _cnt[2].tonumpy()
        magOrig = mag # keep track of all mags for plotting stars not detected by STF
        good = np.where(cnt > 1)[0]
        star = star[good]
        inX = inX[good]
        inY = inY[good]
        outX = outX[good]
        outY = outY[good]
        mag = mag[good]
        dx = inX - outX
        dy = inY - outY

        print 'Total number of stars: %i' % len(dx)
        print asuffix
        print

        # Which planted stars were not detected/matched?
        bad = np.where(cnt == 1)[0]
        py.clf()
        py.figure(figsize=(6,6))
        py.hist(magOrig[bad],bins=200,histtype='step')
        py.xlabel('K')
        py.ylabel('N')
        py.title('Planted stars not detected by STF')
        py.savefig('%s%s/planted_not_detected_mags%s%s.png' % 
            (outdir, adir, suffix, asuffix))

        if removeEdges == True:
            # Remove edges
            good = np.where((inX > 100.) & (inX < 3996.) &
                            (inY > 100.) & (inY < 3996.))[0]
            inX = inX[good]
            inY = inY[good]
            outX = outX[good]
            outY = outY[good]
            mag = mag[good]
            star = star[good]
            dx = inX - outX
            dy = inY - outY
    
            print 'Number of stars after edge removal: %i' % len(good) 

  
        if getHiPrec == True:
            # Also get the original positions in arcsec
            orig = asciidata.open(rootdir + 'input_gc_positions_database_arcsec.lis')
            oNames = orig[0].tonumpy()
            oNames = np.array([oo.strip() for oo in oNames])
            oNames = [oNames[i].replace('irs', '') for i in range(len(oNames))]
            oX = np.sqrt(orig[3].tonumpy()**2)
            oY = np.sqrt(orig[4].tonumpy()**2)
            # Find the ones with positions given to 3 decimals and those to 5 decimals
            oXd = np.array([Decimal(str(oo)) % 1 for oo in oX])
            oYd = np.array([Decimal(str(oo)) % 1 for oo in oY])
    
            oN1 = []
            oN2 = []
            oN3 = [] # 3 decimals
            oN5 = [] # 5 decimals
            oHi = [] # all with greater than 3 decimals
            for ii in range(len(oXd)):
                if ((len(str(oXd[ii])) == 3) | (len(str(oYd[ii])) == 3)):
                    # Find the matching star in the align:
                    foo1 = np.where(star == oNames[ii])[0]
                    oN1 = np.concatenate([oN1, foo1])
                if ((len(str(oXd[ii])) == 4) | (len(str(oYd[ii])) == 4)):
                    # Find the matching star in the align:
                    foo2 = np.where(star == oNames[ii])[0]
                    oN2 = np.concatenate([oN2, foo2])
                if ((len(str(oXd[ii])) == 5) | (len(str(oYd[ii])) == 5)):
                    #oN3 = np.concatenate([oN3, [oNames[ii]]])
                    # Find the matching star in the align:
                    foo3 = np.where(star == oNames[ii])[0]
                    oN3 = np.concatenate([oN3, foo3])
                if ((len(str(oXd[ii])) == 7) | (len(str(oYd[ii])) == 7)):
                    #oN5 = np.concatenate([oN5, [oNames[ii]]])
                    # Find the matching star in the align:
                    foo5 = np.where(star == oNames[ii])[0]
                    oN5 = np.concatenate([oN5, foo5])
                if ((len(str(oXd[ii])) > 5) & (len(str(oYd[ii])) > 5)): # hi precision
                    # Find the matching star in the align:
                    foo = np.where(star == oNames[ii])[0]
                    oHi = np.concatenate([oHi, foo])
    
            oN1 = [np.int(jj) for jj in oN1]
            oN2 = [np.int(jj) for jj in oN2]
            oN3 = [np.int(jj) for jj in oN3]
            oN5 = [np.int(jj) for jj in oN5]
            oHi = [np.int(jj) for jj in oHi]
        
        else:
            oHi = np.arange(0,len(inX),1)
            oHi = [int(ii) for ii in oHi]

        py.close('all')

        sub = np.where((mag > 1) & (mag < 16) & (np.abs(dy) < 0.15) & (np.abs(dx) < 0.15))[0]
        #pdb.set_trace()
        print len(sub)
        top3x = outX[0:3]
        top3y = outY[0:3]
        usetexTrue()
        py.clf()
        py.figure(figsize=(6,6))
        py.subplots_adjust(left=0.15, top=0.92, right=0.95)
        #qvr = py.quiver([outX],[outY],[dx],[dy],units='x',angles='xy',scale=0.0001) # tail of arrow
        qvr = py.quiver([outX[sub]],[outY[sub]],[dx[sub]],[dy[sub]],units='x',angles='xy',scale=0.0001) # tail of arrow
        							                    # at outX,outY
        #py.plot(top3x, top3y, 'r.')
        py.quiverkey(qvr,4020,4240,0.01,'0.01 pix',coordinates='data')
        py.xlabel('X (pixels)')
        py.ylabel('Y (pixels)')
        #py.title('Residual Distortion')
        py.axis([-50,4146,-50,4146])
        py.savefig(outdir + r'%s/delta_in_vs_out%s%s.png' % 
            (adir, suffix, asuffix))
        #py.savefig(outdir + r'%s/delta_in_vs_out%s%s.eps' % 
        #    (adir, suffix, asuffix))
        py.close()

        import matplotlib.font_manager
        prop = matplotlib.font_manager.FontProperties(size=14)

        #binsIn = np.arange(-0.001, 0.001, 0.00001)
        #binsIn = np.arange(-0.03, 0.03, 0.0005)
        #binsIn = np.arange(-1, 1, 0.01)
        binsIn = 500
        py.clf()
        py.figure(figsize=(6,6))
        n1,b1,p1 = py.hist(dx, bins=binsIn, histtype='step',color='r',label='dx')
        n2,b2,p2 = py.hist(dy, bins=binsIn, histtype='step',color='b',label='dy')
        py.legend(numpoints=1,fancybox=True,prop=prop)
        py.xlabel('In - Out (pixels)')
        py.ylabel('N')
        py.savefig(outdir + r'%s/hist_delta_in_vs_out_pix%s%s.png' % 
            (adir, suffix, asuffix))
        #py.savefig(outdir + r'%s/hist_delta_in_vs_out_pix%s%s.eps' % 
        #    (adir, suffix, asuffix))
        py.close()

        xphase = inX - np.floor(inX)
        yphase = inY - np.floor(inY)
#        py.clf()
#        py.figure(figsize=(6,6))
#        py.subplots_adjust(left=0.17, top=0.95, right=0.95)
#        py.plot(xphase,dx,'r.',label='x')
#        py.plot(yphase,dy,'b.',label='y')
#        py.xlabel('Pixel Phase')
#        py.ylabel('In - Out (pix)')
#        py.legend(numpoints=1,fancybox=True,prop=prop)
#        py.axis([-0.01,1.01,-0.1,0.1])
#        py.savefig(outdir + r'%s/delta_in_vs_out_vs_pixPhase%s%s.png' % 
#            (adir, suffix, asuffix))
#        #py.savefig(outdir + r'%s/delta_in_vs_out_vs_pixPhase%s%s.eps' % 
#        #    (adir, suffix, asuffix))
#        py.close()
#
        py.clf()
        py.figure(figsize=(6,6))
        py.semilogy(mag,np.abs(dx),'r.',label='dx')
        py.semilogy(mag,np.abs(dy),'b.',label='dy')
        py.legend(loc=2,numpoints=1,fancybox=True,prop=prop)
        #py.axis([6,20,1e-4,10])
        py.xlabel('K')
        py.ylabel(r'$|$In - Out$|$ (pix)')
        #py.title(r'%s' % suffix)
        py.savefig(outdir + r'%s/in_vs_out_Kmag_pix%s%s.png' % 
            (adir, suffix, asuffix))
        py.close()
        # End Temporary


        print 'All stars: %i' % len(dx)
        print 'Average X offset: %8.5f +/- %8.5f pix (1 sigma = %5.2f uas)' % \
              (dx.mean(), dx.std(ddof=1), dx.std(ddof=1)*4000.)
        print 'Average Y offset: %8.5f +/- %8.5f pix (1 sigma = %5.2f uas)' % \
              (dy.mean(), dy.std(ddof=1), dy.std(ddof=1)*4000.)
        print
        print 'Median |X offset|: %8.5f pix' % np.median(np.abs(dx))
        print 'Median |Y offset|: %8.5f pix' % np.median(np.abs(dy))
        print
        print 'Std. dev. of |X offsets|: %8.5f pix' % np.abs(dx).std(ddof=1)
        print 'Std. dev. of |Y offsets|: %8.5f pix' % np.abs(dy).std(ddof=1)
        kp = np.where((np.abs(dx) < 0.15) & (np.abs(dy) < 0.15))[0]
        kp = np.array([int(ii) for ii in kp])
        #dxp = np.array([np.abs(dx[ii]) for ii in kp])
        #dyp = np.array([np.abs(dy[ii]) for ii in kp])
        dxp = np.array([dx[ii] for ii in kp])
        dyp = np.array([dy[ii] for ii in kp])
        print
        print 'Remove outliers, left with %i stars:' % len(kp)
        print 'Average X offset: %8.5f +/- %8.5f pix (1 sigma = %5.2f uas)' % \
              (dxp.mean(), dxp.std(ddof=1), dxp.std(ddof=1)*4000.)
        print 'Average Y offset: %8.5f +/- %8.5f pix (1 sigma = %5.2f uas)' % \
              (dyp.mean(), dyp.std(ddof=1), dyp.std(ddof=1)*4000.)
        print
        print 'Std. dev. of |X offsets|: %8.5f pix (1 sigma = %5.2f uas)' % \
              (np.abs(dxp).std(ddof=1), np.abs(dxp).std(ddof=1)*4000.)
        print 'Std. dev. of |Y offsets|: %8.5f pix (1 sigma = %5.2f uas)' % \
              (np.abs(dyp).std(ddof=1), np.abs(dyp).std(ddof=1)*4000.)
        print '---------------------------------------'

        if getHiPrec == True:
            print
            print '---------------------------------------'
    
            print 'Stars known to >=5 decimal places: %i' % len(oHi)
            print 'Average X offset: %8.5f +/- %8.5f pix (1 sigma = %5.2f uas)' % \
                  (dx[oHi].mean(), dx[oHi].std(ddof=1), dx[oHi].std(ddof=1)*4000.)
            print 'Average Y offset: %8.5f +/- %8.5f pix (1 sigma = %5.2f uas)' % \
                  (dy[oHi].mean(), dy[oHi].std(ddof=1), dy[oHi].std(ddof=1)*4000.)
            print
            print 'Median |X offset|: %8.5f pix' % np.median(np.abs(dx[oHi]))
            print 'Median |Y offset|: %8.5f pix' % np.median(np.abs(dy[oHi]))
            print
            print 'Std. dev. of |X offsets|: %8.5f pix' % np.abs(dx[oHi]).std(ddof=1)
            print 'Std. dev. of |Y offsets|: %8.5f pix' % np.abs(dy[oHi]).std(ddof=1)
    
            kp = np.where((np.abs(dx[oHi]) < 0.025) & (np.abs(dy[oHi]) < 0.025))[0]
            kp = np.array([int(ii) for ii in kp])
            #dxp = np.array([np.abs(dx[oHi[ii]]) for ii in kp])
            #dyp = np.array([np.abs(dy[oHi[ii]]) for ii in kp])
            dxp = np.array([dx[oHi[ii]] for ii in kp])
            dyp = np.array([dy[oHi[ii]] for ii in kp])
            print
            print 'Remove outliers (of the set of remaining stars), left with %i stars:' % len(kp)
            print 'Std. dev. of X offsets: %8.5f pix' % dxp.std(ddof=1)
            print 'Std. dev. of Y offsets: %8.5f pix' % dyp.std(ddof=1)

        # Put everything in milli-arcseconds for the histogram
        dx *= 4.0
        dy *= 4.0
    
        binsIn = 500
        #binsIn = np.arange(-0.12, 0.12, 0.001)
        #binsIn = np.arange(-4, 4, 0.1)
        py.clf()
        py.figure(figsize=(6,6))
        py.hist(dx, bins=binsIn, histtype='step',color='r',label='dx')
        py.hist(dy, bins=binsIn, histtype='step',color='b',label='dy')
        py.legend(numpoints=1,fancybox=True,prop=prop)
        py.xlabel('In - Out (mas)')
        py.ylabel('N')
        #py.title(r'%s' % suffix)
        py.savefig(outdir + r'%s/hist_delta_in_vs_out%s%s.png' % 
            (adir, suffix, asuffix))
        py.close()
    
        gdX = np.where(np.abs(dx) < 0.5)[0]
        gdY = np.where(np.abs(dy) < 0.5)[0]
        py.clf()
        py.figure(figsize=(6,6))
        #py.semilogy(mag[gdX],np.abs(dx[gdX]),'r.',label='dx')
        #py.semilogy(mag[gdY],np.abs(dy[gdY]),'b.',label='dy')
        py.semilogy(mag,np.abs(dx),'r.',label='dx')
        py.semilogy(mag,np.abs(dy),'b.',label='dy')
        py.legend(loc=2,numpoints=1,fancybox=True,prop=prop)
        #py.axis([6,20,1e-4,10])
        py.xlabel('K')
        py.ylabel(r'$|$In - Out$|$ (mas)')
        #py.title(r'%s' % suffix)
        py.savefig(outdir + r'%s/in_vs_out_Kmag%s%s.png' % 
            (adir, suffix, asuffix))
        py.close()

        #py.clf()
        #py.figure(figsize=(6,6))
        #
        #py.savefig()
        #py.close()
    
        usetexFalse()

        brt = np.where(mag < 18.0)[0]
        print
        print 'Bright stars (Kp < 18):'
        print 'Average X offset: %5.3f +- %5.3f mas' % \
              ((dx[brt]).mean(), (dx[brt]).std(ddof=1))
        print 'Average Y offset: %5.3f +- %5.3f mas' % \
              ((dy[brt]).mean(), (dy[brt]).std(ddof=1))
        print
        print 'Median X offset: %5.3f mas' % np.median(dx[brt])
        print 'Median Y offset: %5.3f mas' % np.median(dy[brt])

def confusion(alnDir='align_30042013/',single=True, interactive=False,
              suffix='',oplotKeck=False,plotRadialConfusion=False,scale=0.01):
    """
    Set interactive to True to select specific region of image
    and/or specific magnitude stars to plot deltas for.
    """

    rootdir = root + alnDir
    #outdir = rootdir + 'starfinder/in_vs_out/'
    outdir = rootdir + 'starfinder_knownPSF/in_vs_out/'
    #outdir = rootdir

    #font = {'family':'sans-serif','size': 18, 'weight': 'bold'}
    #py.rc('font', **font)

    if single == True:
        adir = 'single_maps'
        #adir = ''
        asuffix = '1'
        #asuffix = ''
        #asuffix = '1_corr'
    else:
        adir = 'deepmap'
        asuffix = '_m'
        #asuffix = ''

    aln = asciidata.open('%s%s/align%s.origpos' % (outdir, adir,asuffix))
    #aln = asciidata.open('%s%s/align%s.pos' % (outdir, adir,asuffix))
    starOrig = aln[0].tonumpy()
    inXOrig = aln[1].tonumpy()
    inYOrig = aln[2].tonumpy()
    outXOrig = aln[3].tonumpy()
    outYOrig = aln[4].tonumpy()
    
    _cnt = asciidata.open('%s%s/align%s.vel' % (outdir, adir,asuffix))
    cntOrig = _cnt[1].tonumpy()
    magOrig = _cnt[2].tonumpy()

    _mag = asciidata.open('%s%s/align%s.mag' % (outdir, adir,asuffix))
    starNmag = _mag[0].tonumpy()
    inMagOrig = _mag[4].tonumpy()
    outMagOrig = _mag[5].tonumpy()

    good = np.where(cntOrig > 1)[0]
    star = starOrig[good]
    inX = inXOrig[good]
    inY = inYOrig[good]
    outX = outXOrig[good]
    outY = outYOrig[good]
    mag = magOrig[good]
    inMag = inMagOrig[good]
    outMag = outMagOrig[good]
    cnt = cntOrig[good]
    dx = inX - outX
    dy = inY - outY

    print 'Total number of stars: %i' % len(dx)
    print asuffix
    print

    mismatch = np.where(np.abs(inMag - outMag) >= 1000.0)[0]
    match = np.where(np.abs(inMag - outMag) < 1000.0)[0]
    print 'Number of possible mismatches: %i' % len(mismatch)
    print
    star = star[match]
    inX = inX[match]
    inY = inY[match]
    outX = outX[match]
    outY = outY[match]
    mag = mag[match]
    inMag = inMag[match]
    outMag = outMag[match]
    cnt = cnt[match]
    dx = inX - outX
    dy = inY - outY

    if plotRadialConfusion == True:
        # Put the positions in arcsec from SgrA*
        sgra = [512.0, 612.0] # pixel position in simulation
        xasec = -1.*(inX - sgra[0]) * scale
        yasec = (inY - sgra[1]) * scale
        r2d_asec = np.hypot(xasec,yasec)

    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=16)

    if interactive == True:
        xctr = raw_input("Define box center, X position (pix): ")
        yctr = raw_input("Define box center, Y position (pix): ")
        boxsize = raw_input("Define size of box (pix): ")
        deltaLo = raw_input("Size of offsets to examine, lower limit (pix)? ")
        deltaHi = raw_input("Size of offsets to examine, upper limit (pix)? ")

        inReg = open('%s/inputPos_%s_%s_%s.reg' % \
                     (outdir, str(xctr), str(yctr), str(boxsize)), 'w')
        outReg = open('%s/outputPos_%s_%s_%s.reg' % \
                      (outdir, str(xctr), str(yctr), str(boxsize)), 'w')
        fmt2 = 'image;circle(%10.5f,%10.5f,5)\n'

        xctr = float(xctr)
        yctr = float(yctr)
        rng = float(boxsize)/2.
        deltaLo = float(deltaLo)
        deltaHi = float(deltaHi)

        # For the stars that fall in this box, find those with offsets of size delta
        xys = np.column_stack((inX, inY))
        xvrt = np.array([xctr-rng,xctr-rng,xctr+rng,xctr+rng])
        yvrt = np.array([yctr+rng,yctr-rng,yctr-rng,yctr+rng])
        verts = np.column_stack((xvrt, yvrt))
        mask = nx.points_inside_poly(xys, verts)

        # Write out a regions file
        inReg.write('# Region file format: DS9 version 3.0\n')
        inReg.write('global color=black font="helvetica 10 normal" edit=1 move=1 delete=1 include=1 fixed=0 width=2\n')
        outReg.write('# Region file format: DS9 version 3.0\n')
        outReg.write('global color=white font="helvetica 10 normal" edit=1 move=1 delete=1 include=1 fixed=0 width=2\n')
        for rr in range(len(inX)):
            if mask[rr] == True:
                #if ((dx[rr] > delta) | (dy[rr] > delta)):
                if ((np.abs(dx[rr]) > deltaLo) & (np.abs(dy[rr]) > deltaLo) &
                    (np.abs(dx[rr]) < deltaHi) & (np.abs(dy[rr]) < deltaHi)):
                    inReg.write(fmt2 % (inX[rr]+1.,inY[rr]+1.)) # add 1 pix for ds9
                    outReg.write(fmt2 % (outX[rr]+1.,outY[rr]+1.)) # add 1 pix for ds9
                else:
                    continue

        inReg.close()
        outReg.close()

    usetexTrue()

    # Which planted stars were not detected/matched?
    bad = np.where(cntOrig == 1)[0]
    py.figure(figsize=(6,6))
    py.subplots_adjust(left=0.15, top=0.9, right=0.95, bottom=0.12)
    py.clf()
    py.hist(magOrig[bad],bins=np.arange(6,25,1),histtype='step',color='k')
    py.xlabel('K')
    py.ylabel('N')
    py.title('Planted stars not detected by STF')
    py.savefig('%s%s/confusion_planted_not_detected_mags%s%s.png' % 
        (outdir, adir, suffix, asuffix))
    py.close()

    magBins = np.arange(10.,mag.max()+1)
    print mag.max()
    #rows = np.round(len(magBins)/4.)
    #cols = np.round(len(magBins)/4.)
    rows = np.round(len(magBins)/3.) # keck
    cols = np.round(len(magBins)/3.) # keck
    rmsX = np.zeros(len(magBins),dtype=float)
    rmsY = np.zeros(len(magBins),dtype=float)
    rmsR = np.zeros(len(magBins),dtype=float)
    comp = np.zeros(len(magBins),dtype=float)
    #py.figure(1)
    fig1 = py.figure(1,figsize=(12,12))
    py.clf()
    py.subplots_adjust(left=0.02, top=0.97, right=0.98, bottom=0.02,
                       wspace=0.05, hspace=0.2)
    #py.figure(2)
    fig2 = py.figure(2,figsize=(12,12))
    py.clf()
    py.subplots_adjust(left=0.05, top=0.92, right=0.98, bottom=0.02,
                       wspace=0.2, hspace=0.25)
    for mm in range(len(magBins)):
        print 'Plotting subplot %s' % str(mm+1)
        if mm == 0:
            midx = np.where(mag < magBins[0]+1)[0]
            tl = r'K $<$ %2i (N = %i)' % (magBins[mm]+1, len(midx))
            midxOrig = np.where(magOrig < magBins[0]+1)[0]
        elif mm == len(magBins)-1:
            midx = np.where(mag > magBins[-1])
            tl = r'K $>$ %2i (N = %i)' % (magBins[-1], len(midx))
            midxOrig = np.where(magOrig > magBins[-1])
        else:
            midx = np.where((mag > magBins[mm]) & (mag < magBins[mm+1]))[0]
            tl = 'K = %2i-%2i (N = %i)' % (magBins[mm], magBins[mm+1], len(midx))
            midxOrig = np.where((magOrig > magBins[mm]) & (magOrig < magBins[mm+1]))[0]
        print tl
        #print len(midx), magBins[mm]

        #rmsX[mm] = dx[midx].std(ddof=1)
        #rmsY[mm] = dy[midx].std(ddof=1)
        #rmsR[mm] = (rmsX[mm] + rmsY[mm]) / 2.0 # take the average
        #print '1 sigma |in - out| = %6.2f uas' % (rmsR[mm]*4000.) # uas

        # Do some sigma-clipping
        dxM = dx[midx]
        dyM = dy[midx]
        sigX = dxM.std(ddof=1)
        sigY = dyM.std(ddof=1)
        idx = np.where(np.abs(dxM) < (dxM.mean()+3.0*sigX))[0]
        idy = np.where(np.abs(dyM) < (dyM.mean()+3.0*sigY))[0]
        dxM = dxM[idx]
        dyM = dyM[idy]
        rmsX[mm] = dxM.std(ddof=1)
        rmsY[mm] = dyM.std(ddof=1)
        rmsR[mm] = (rmsX[mm] + rmsY[mm]) / 2.0 # take the average
        print '1 sigma |in - out| = %6.2f uas' % (rmsR[mm]*4000.) # uas
        
        if len(idx) < 2:
            print 
            print 'Too few stars in this magnitude bin (%s)' % (tl)
            continue

        # completeness
        comp[mm] = len(midx)*1.0 / len(midxOrig)*1.0

        py.figure(1)
        fig1.add_subplot(rows,cols,mm+1)
        qvr = py.quiver([outX[midx]],[outY[midx]],[dx[midx]],[dy[midx]],units='x',
                        angles='xy',scale=0.03) # tail of arrow at outX,outY
        thePlot = py.gca()
        py.title(tl,fontsize=10)
        py.axis([0,4096,0,4096])
        py.setp( thePlot.get_xticklabels(), fontsize=12, visible=False)
        py.setp( thePlot.get_yticklabels(), fontsize=12, visible=False)

        py.figure(2)
        fig2.add_subplot(rows,cols,mm+1)
        dxmin = dxM.min()
        dxmax = dxM.max()
        dymin = dyM.min()
        dymax = dyM.max()
        xinc = (dxmax - dxmin) / 100.0
        yinc = (dymax - dymin) / 100.0
        py.rc('xtick', labelsize=10)
        py.rc('ytick', labelsize=10)
        #thePlot = py.gca()
        if mm == 0:
            py.xticks([dxmin+3.*xinc, 0, dxmax+5.*xinc])
        elif mm < 5:
            py.xticks([dxmin+3.*xinc, 0, dymax-3.*yinc])
        else:
            py.xticks([dxmin+7.*xinc, 0, dxmax-7.*xinc])
        #thePlot.xaxis.set_major_formatter(py.FormatStrFormatter('%6.3e'))
        #n1,b1,p1 = py.hist(dx[midx],bins=np.arange(dxmin,dxmax,xinc),
        n1,b1,p1 = py.hist(dxM,bins=np.arange(dxmin,dxmax,xinc),
                histtype='step',color='r')
        #n2,b2,p2 = py.hist(dy[midx],bins=np.arange(dymin,dymax,yinc),
        n2,b2,p2 = py.hist(dyM,bins=np.arange(dymin,dymax,yinc),
                histtype='step',color='b')
        py.title(tl,fontsize=10)
        py.axis([min(dxM.min(),dyM.min()), max(dxM.max(),dyM.max()),
                 0,max(n1.max(),n2.max())+2])

    if plotRadialConfusion == True:
        # Get confusion error in radial bins
        radStep = 1.0
        radBins = np.arange(0.0, 8.0, radStep)
        rmsXR = np.zeros(len(radBins),dtype=float)
        rmsYR = np.zeros(len(radBins),dtype=float)
        rmsRR = np.zeros(len(radBins),dtype=float)
        for rr in range(len(radBins)):
            rMin = radBins[rr] - (radStep / 2.0)
            rMax = radBins[rr] + (radStep / 2.0)
            # the following condition is specific for Boehle et al. 2014 (S0-38):
            #idxR = np.where((r2d_asec >= rMin) & (r2d_asec < rMax) & (inMag > 16.5) & (inMag < 17.5))[0]
            # the following condition is specific for Boehle et al. 2014 (S0-2):
            idxR = np.where((r2d_asec >= rMin) & (r2d_asec < rMax) & (inMag > 13.5) & (inMag < 14.5))[0]
    
            dxR = dx[idxR]
            dyR = dy[idxR]
            rmsXR[rr] = dxR.std(ddof=1)
            rmsYR[rr] = dyR.std(ddof=1)
            rmsRR[rr] = (rmsXR[rr] + rmsYR[rr]) / 2.0 # take the average
            print 'R2D = (%4.1f arcsec): 1 sigma |in - out| = %6.2f mas (N=%i)' % \
                  ((radBins[rr]+0.5), (rmsRR[rr]*4.), len(idxR)) # mas
            

    py.figure(1)
    py.savefig('%s%s/confusion_delta_in_vs_out%s%s.png' % 
        (outdir, adir, suffix, asuffix))
    py.close(1)

    py.figure(2)
    py.suptitle('In - Out (pixel) binned by K mag')
    py.savefig('%s%s/confusion_hist_delta_in_vs_out%s%s.png' % 
        (outdir, adir, suffix, asuffix))
    py.close(2)

    # Do we want to overplot the Keck confusion estimates?
    # This must have been printed out to a text file already and in uas
    if oplotKeck == True:
        keck = asciidata.open('%s/KeckSims/align_20052013/starfinder/in_vs_out/single_maps/rms_in_vs_out_mag.txt' %
                              root)
        magBinKeck = keck[0].tonumpy()
        rmsRkeck = keck[1].tonumpy()
        keck = asciidata.open('%s/KeckSims/align_20052013/starfinder/in_vs_out/single_maps/completeness.txt' %
                              root)
        magCompKeck = keck[0].tonumpy()
        compKeck = keck[1].tonumpy()
        keck = asciidata.open('%s/KeckSims/align_20052013/starfinder/in_vs_out/single_maps/allmags.txt' %
                              root)
        magKeck = keck[0].tonumpy()

    usetexTrue()
    py.figure(3)
    py.clf()
    py.figure(figsize=(6,6))
    py.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.12)
    #py.semilogy(magBins+0.5, rmsR*4000., 'r-', lw=2, label='TMT')
    #py.semilogy(magBins+0.5, rmsR*4000., 'r.', ms=8)
    #py.plot(magBins[:-2]+0.5, rmsR[:-2]*4., 'r-', lw=2, label='NGAO/Keck')
    #py.plot(magBins[:-2]+0.5, rmsR[:-2]*4., 'r.', ms=8)
    py.plot(magBins[:-2]+0.5, rmsR[:-2]*4., 'k-', lw=2, label='Keck') # speckle holography simulation
    py.plot(magBins[:-2]+0.5, rmsR[:-2]*4., 'k.', ms=8)
    if oplotKeck == True:
        py.plot(magBinKeck+0.5, rmsRkeck/1.e3, 'b--', lw=1.5, label='Current LGSAO/Keck')
        py.plot(magBinKeck+0.5, rmsRkeck/1.e3, 'b.', ms=8)
        #py.semilogy(magBinKeck+0.5, rmsRkeck/1.e3, 'b--', lw=1.5, label='Current AO')
        #py.semilogy(magBinKeck+0.5, rmsRkeck/1.e3, 'b.', ms=8)
        py.legend(numpoints=1,fancybox=True,loc=2,prop=prop)
    thePlot = py.gca()
    py.setp( thePlot.get_xticklabels(), fontsize=18, weight='heavy')
    py.setp( thePlot.get_yticklabels(), fontsize=18, weight='heavy')
    py.axis([10,21,1.e-2,3])
    #py.axis([10,19,0,1.0]) # spec holography
    py.xlabel(r'{\bf K mag}',fontsize=18)
    #py.ylabel(r'$\sigma_{\rm confusion}$ ($\mu$as)')
    #py.ylabel(r'$\sigma_{\rm confusion}$ (mas)')
    py.ylabel(r'{\bf Astrometric Uncertainty from Source Confusion (mas)}',
              fontsize=14,multialignment='center')
    py.savefig('%s%s/rms_in_vs_out_mag%s%s.png' % \
               (outdir, adir, suffix, asuffix))
    #py.savefig('%s%s/rms_in_vs_out_mag%s%s.eps' % \
    #           (outdir, adir, suffix, asuffix))
    py.close(3)
    print ''
    for mm in range(len(magBins)):
        print '%4.1f   %10.5f  %5.3f' % (magBins[mm]+0.5, rmsR[mm]*4., comp[mm])
    print
    if oplotKeck == True:
        print 'Current AO'
        for mm in range(len(magBinKeck)):
            print '%4.1f   %10.5f  %5.3f' % (magBinKeck[mm]+0.5, rmsRkeck[mm]/1.e3, compKeck[mm])

    if plotRadialConfusion == True:
        py.figure(8)
        py.clf()
        py.figure(figsize=(6,6))
        py.plot(radBins+0.5, rmsRR*4., 'k-', lw=2, label='Keck') # speckle holography simulation
        py.plot(radBins+0.5, rmsRR*4., 'k.', ms=8)
        thePlot = py.gca()
        py.setp( thePlot.get_xticklabels(), fontsize=18, weight='heavy')
        py.setp( thePlot.get_yticklabels(), fontsize=18, weight='heavy')
        #py.axis([0,10,0.2,1.5]) # for log scale
        py.axis([0,10,0,1.2])
        py.xlabel(r'{\bf Radius (arcsec)}',fontsize=18)
        py.ylabel(r'{\bf Astrometric Uncertainty from Source Confusion (mas)}',
                  fontsize=14,multialignment='center')
        py.savefig('%s%s/rms_in_vs_out_radial%s%s.png' % \
                   (outdir, adir, suffix, asuffix))
        py.savefig('%s%s/rms_in_vs_out_radial%s%s.eps' % \
                   (outdir, adir, suffix, asuffix))
        py.close(8)

    py.figure(4)
    py.clf()
    py.figure(figsize=(6,6))
    py.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.12)
    #n1,b1,p1 = py.hist(mag,bins=np.arange(6,26,0.25),histtype='step',color='b',label='Matched')
    #n2,b2,p2 = py.hist(magOrig[bad],bins=np.arange(6,26,0.25),histtype='step',color='r',label='Not Matched')
    #n2,b2,p2 = py.hist(magOrig,bins=np.arange(6,24,0.5),histtype='step',color='r',
    #                   label='Planted')
    n1,b1,p1 = py.hist(mag,bins=np.arange(6,24.5,0.5),histtype='step',color='k',
                       label='TMT')
    if oplotKeck == True:
        n2,b2,p2 = py.hist(magKeck,bins=np.arange(6,24.5,0.5),histtype='step',color='b',
                           label='Keck')
    py.xlabel('K Magnitude')
    py.ylabel('Number of Detected Stars')
    #py.legend(numpoints=1,fancybox=True,loc=2,prop=prop)
    py.axis([8,22,0,350])
    #py.axis([5,24,0,max(n1.max(),n2.max())+500])
    #py.savefig(r'%s%s/plots/mag_detected_notdetected%s%s.png' % \
    #           (outdir, adir, suffix, asuffix))
    #py.savefig(r'%s%s/plots/mag_detected_notdetected%s%s.eps' % \
    #           (outdir, adir, suffix, asuffix))
    py.savefig(r'%s%s/mag_detected%s%s.png' % \
               (outdir, adir, suffix, asuffix))
    #py.savefig(r'%s%s/plots/mag_detected_planted%s%s.png' % \
    #           (outdir, adir, suffix, asuffix))
    #py.savefig(r'%s%s/plots/mag_detected_planted%s%s.eps' % \
    #           (outdir, adir, suffix, asuffix))
    py.close(4)

    # Plot completeness curve
    py.figure(6)
    py.clf()
    py.figure(figsize=(6,6))
    py.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.12)
    py.plot(magBins+0.5, comp, 'k.')
    py.plot(magBins+0.5, comp, 'k-', label='NGAO/Keck')
    if oplotKeck == True:
        py.plot(magCompKeck+0.5, compKeck, 'b--', lw=1.5, label='Current LGSAO/Keck')
        py.plot(magCompKeck+0.5, compKeck, 'b.', ms=8)
        py.legend(numpoints=1,fancybox=True,loc=3,prop=prop)
    py.xlabel('Brightness (K-Magnitude)',fontsize=20)
    py.ylabel('Completeness',fontsize=20)
    py.axis([9,25,0,1.05])
    py.savefig(r'%s%s/completeness%s%s.png' % \
               (outdir, adir, suffix, asuffix))
    py.close(6)

    if oplotKeck == True:
        py.figure(7)
        py.clf()
        py.figure(figsize=(6,6))
        py.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.12)
        n1,b1,p1 = py.hist(mag,bins=np.arange(6,26,0.25),histtype='step',color='r',label='NGAO/Keck',lw=2)
        #n2,b2,p2 = py.hist(magOrig,bins=np.arange(6,26,0.25),histtype='step',color='r',linestyle='dashed',lw=1.5)
        n3,b3,p3 = py.hist(magKeck,bins=np.arange(6,26,0.25),histtype='step',color='b',label='Current LGSAO/Keck',lw=2)
        thePlot = py.gca()
        py.setp( thePlot.get_xticklabels(), fontsize=18, weight='heavy')
        py.setp( thePlot.get_yticklabels(), fontsize=18, weight='heavy')
        #thePlot.set_xticklabels(thePlot.get_xticks(), font)
        #thePlot.set_yticklabels(thePlot.get_yticks(), font)
        py.xlabel(r'{\bf K mag}',fontsize=18)
        py.ylabel(r'{\bf Number of Stars Detected}',fontsize=18)
        py.legend(numpoints=1,fancybox=True,loc=2,prop=prop)
        py.axis([5,25,0,max(n1.max(),n3.max())+100])
        py.savefig(r'%s%s/histMag_compareKeck%s%s.png' % \
                   (outdir, adir, suffix, asuffix))
        py.close(7)
    

    # Read in the planted star information to make a KLF of what we planted
    #plnt = asciidata.open('%s/known_and_simulated_2012.0_pixels.lis' % rootdir)
    #magPlnt = plnt[1].tonumpy()
    #py.figure(5)
    #py.clf()
    #py.figure(figsize=(6,6))
    #py.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.12)
    #nn,bb,pp = py.hist(magPlnt,bins=np.arange(6,26,0.25),histtype='step',color='k')
    #py.xlabel('K')
    #py.ylabel('N')
    #py.savefig(r'%s%s/plots/mag_planted%s%s.png' % \
    #           (outdir, adir, suffix, asuffix))
    #py.savefig(r'%s%s/plots/mag_planted%s%s.eps' % \
    #           (outdir, adir, suffix, asuffix))
    #py.close(5)

    usetexFalse()

    # temp for Keck data files (for overplotting)
    #for mm in range(len(magBins)):
    #    print magBins[mm], rmsR[mm]*4000.
    #for mm in range(len(magBins)):
    #    print magBins[mm], comp[mm]
    #out = open('%s/%s/allmags.txt' % (outdir, adir), 'w')
    #for ii in range(len(mag)):
    #    out.write('%8.3f\n' % mag[ii])
    #out.close()
        


def errVsEpoch(alnDir='KeckSims/align_05062013/align/',
               align='align/align_d_rms_1000_abs_t',refEpoch=7):
    """
    plot median alignment error and centroiding error vs. epoch.
    """

    title = alnDir.split('/')[1][6:]
    print 'align: %s' % title

    # Get the star set info
    s = starset.StarSet(root+alnDir+align)

    _epFile = root + alnDir + 'scripts/epochsInfo.txt'
    table = asciidata.open(_epFile)
    epoch = table[0].tonumpy()
    refdate = epoch[refEpoch]
    isAO     = [table[2][ss] == 1 for ss in range(table.nrows)]
    numEpochs = len(epoch)
    dfile = open(root + alnDir + 'align/align_d_rms.date')
    dates = dfile.readline().split()
    date = np.zeros(len(epoch),dtype=float)
    for dd in range(numEpochs):
        date[dd] = dates[dd]

    mag = s.getArray('mag')
    names = s.getArray('name')
    names = np.array(names)

    eps = []
    dates = []
    rerrA_all = []
    rerrC_all = []

    usetexTrue()
    for ee in range(numEpochs):
        if ee == refEpoch:
            print 'ref epoch: %f' % date[ee]
            continue
        xerrA = s.getArrayFromEpoch(ee, 'xerr_a') * 1.e6 # align errors (uas)
        yerrA = s.getArrayFromEpoch(ee, 'yerr_a') * 1.e6
        rerrA = (xerrA + yerrA) / 2.0
        
        xerrC = s.getArrayFromEpoch(ee, 'xerr_p') * 1.e6 # centroid errors (uas)
        yerrC = s.getArrayFromEpoch(ee, 'yerr_p') * 1.e6
        rerrC = (xerrC + yerrC) / 2.0

        x = s.getArrayFromEpoch(ee, 'x')
        y = s.getArrayFromEpoch(ee, 'y')

        # Only get non-zero values
        nz = np.where(xerrA > 0.0)[0]

        if len(nz) > 0:
            name = np.array([names[nn] for nn in nz])
            mm = mag[nz]
            x = x[nz]
            y = y[nz]
            xerrA = xerrA[nz] 
            yerrA = yerrA[nz]
            xerrC = xerrC[nz]
            yerrC = yerrC[nz]
            rerrA = rerrA[nz]
            rerrC = rerrC[nz]
            r = np.hypot(x, y)
        else:
            # Must be the reference epoch
            print 'No stars found in this epoch!!'
            continue

        rerrA_all = np.concatenate([rerrA_all, [np.median(rerrA)]]) # align errors
        eps = np.concatenate([eps, [ee]])
        dates = np.concatenate([dates, [date[ee]]])

        # Get the stars brighter than K=15, r<4"
        kp = np.where((r < 4.0) & (mm < 15))[0]
        rerrC_all = np.concatenate([rerrC_all, [np.median(rerrC[kp])]]) # centroid errors

    def fitfuncLine(p, fjac=None, x=None, y=None, err=None):
        # linear fit
        fun = p[0] + p[1]*x
        
        # deviations from the model
        deviates = (y - fun)
        return [0, deviates]

    def fitLine(p0=None,data=None,quiet=0):
        """Fits line using mpfit.

        Inputs should be given as p0=[A, alpha] and
        data=[date,alignment err]. Returns object of class
        mpfit.
        """
    
        print 'Initial Guess:'
        print '   A     = %6.2f' % p0[0]
        print '   alpha = %5.3f' % p0[1]
    
        # Set data to be passed to fit
        functargs = {'x':data[0],'y':data[1]}
    
        # Set initial values and limits (no limits on parameters)
        pinfo = [{'value':0,'fixed':0,'limited':[0,0],
    	      'limits':[0,0]}]*len(p0)
        for ii in range(len(p0)):
            pinfo[ii]['value'] = p0[ii]
    
        m = nmpfit_sy.mpfit(fitfuncLine, p0, functkw=functargs, parinfo=pinfo,
    		    quiet=quiet)
        if (m.status <= 0):
            print 'Error message = ', m.errmsg
    
        p = m.params                # Best-fit parameters
        perr = m.perror             # Error in parameter fits
                                    # from covariance matrix
    
        m.dof = len(data[0])-len(p) # Number of degrees of freedom
        Rchi2 = m.fnorm/m.dof       # Reduced Chi^2 statistic

        print 'Final Solution:'
        #print '   A      = %6.2f +/- %5.2f' % (p[0],perr[0])
        #print '   alpha  = %5.3f +/- %5.3f' % (p[1],perr[1])
        print '   A      = %6.2f ' % p[0]
        print '   alpha  = %5.3f ' % p[1]
        print '   chi^2  = %5.2f' % m.fnorm
        print '   Rchi^2 = %5.2f' % Rchi2
    
        return m
        
    print
    print '  Average centroid = %5.3f uas' % (rerrC_all.mean())
    print '  Median centroid = %5.3f uas' % (np.median(rerrC_all))
    print '  Average align = %5.3f uas' % (rerrA_all.mean())
    print '  Median align = %5.3f uas' % (np.median(rerrA_all))
    print

    # Fit a line to the alignment errors from ref epoch and back
    # set the ref epoch to be 0.0 and fit data to that
    p0 = [0,10]
    #fitdate1 = np.arange(0,(epoch.max()-refEpoch))
    pfit1 = fitLine(p0,[dates[0:refEpoch], rerrA_all[0:refEpoch]],1)
    fparams1 = pfit1.params
    int1 = fparams1[0]
    slp1 = fparams1[1]
    # Fit a line to the alignment errors from ref epoch forward
    p0 = [0,10]
    #fitdate2 = np.arange(0,(epoch.max()-refdate)-1)
    pfit2 = fitLine(p0,[dates[refEpoch:], rerrA_all[refEpoch:]],1)
    fparams2 = pfit2.params
    int2 = fparams2[0]
    slp2 = fparams2[1]

    #ave_int = (int1 + int2) / 2.0
    ave_slp = (-slp1 + slp2) / 2.0
    print
    print 'average slope: %f' % ave_slp
    print 'range of slopes: %f' % np.abs((np.abs(slp1) - np.abs(slp2)))
    #print 'average intercept: %f' % ave_int
    # plot the average slopes (with opposite signs)
    fitLine1 = int1 - slp1*dates[0:refEpoch]
    fitLine2 = int2 + slp2*dates[refEpoch:]

    # Write out to a file also
    outfile = root + alnDir + 'plots/align_cntrd.txt'
    out = open(outfile,'w')
    hdr = '#%6s  %6s  %7s  %5s  %6s  %7s\n'
    out.write(hdr % ('<cnt>','sd_cnt','med_cnt','<aln>','sd_aln','med_aln'))
    fmt = '%6.2f  %6.2f  %7.2f  %5.2f  %6.2f  %7.2f\n'
    out.write(fmt % (rerrC_all.mean(),rerrC_all.std(ddof=1),np.median(rerrC_all),
                     rerrA_all.mean(),rerrA_all.std(ddof=1),np.median(rerrA_all)))
    out.close()

    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=14)

    # Plot align errors vs. epoch
    py.close('all')
    py.figure(figsize=(6,6))
    py.subplots_adjust(left=0.14, right=0.95, top=0.9)
    py.clf()
    py.plot(dates, rerrA_all, 'k.',ms=9,label=r'$\sigma_{aln}$')
    py.plot(dates, rerrC_all, 'rs',mfc='None',mec='r',mew=1.5,ms=5,label=r'$\sigma_{cnt}$')
    #py.plot(dates[0:refEpoch],fitLine1,'k--')
    #py.plot(dates[refEpoch:],fitLine,'k--')
    #py.plot(fitdate2+refEpoch,fitLine2,'k--')
    #py.plot(dates[refEpoch:],fitLine,'k--')
    py.title(r'%s' % title)
    #py.axis([2005,2019,10,300])
    #py.axis([2005,2019,10,1.e3])
    py.legend(numpoints=1,fancybox=True,loc=2)
    py.xlabel('Date (year)')
    py.ylabel(r'Positional Error in Individual Maps ($\mu$as)')
    py.savefig(root + alnDir + 'plots/align_cntrd_vs_epoch.png')
    py.close()

    usetexFalse()



def convert_to_arcsec(infile='align_16112012/input_gc_positions.lis',scale=0.004):
    """
    Read in a star list with positions in pixels, convert to arcsec, and
    write out to a new star list.
    """

    infile = root + infile
    pixData = asciidata.open(infile)
    name = pixData[0]._data
    mag = pixData[1].tonumpy()
    year = pixData[2].tonumpy()
    xPix = pixData[3].tonumpy()
    yPix = pixData[4].tonumpy()
    snr = pixData[5].tonumpy()
    corr = pixData[6].tonumpy()
    nf = pixData[7].tonumpy()
    ff = pixData[8].tonumpy()

    xhalf = 2048.0
    yhalf = 2048.0
    x = (xPix - xhalf) * scale * -1
    y = (yPix - yhalf) * scale
   

    outfile = infile.split('.')[0] + '_arcsec.lis'
    oo = open(outfile,'w')
 
    fmt = '%12s  %6.3f  %8.3f  %10.5f  %10.5f  %4.2f  %4.2f  %1i  %5.3f\n'
    for ii in range(len(name)):
        oo.write(fmt % (name[ii],mag[ii],year[ii],x[ii],y[ii],
                            snr[ii],corr[ii],nf[ii],ff[ii]))
    
    oo.close()
    
def convert_to_pixels(infile='align_16112012/gc_database_arcsec.lis',keck=False,
                      restrictFOV=True):
    """
    Read in a star list with positions in arcsec, convert to pixels in IRIS
    detector coordinates, with Sgr A* at the center. Write out to a new star list.

    infile must have 'arcsec' in name.
    """

    infile = root + infile
    asecData = asciidata.open(infile)
    name = asecData[0]._data
    mag = asecData[1].tonumpy()
    year = asecData[2].tonumpy()
    xAsec = asecData[3].tonumpy()
    yAsec = asecData[4].tonumpy()
    snr = asecData[5].tonumpy()
    corr = asecData[6].tonumpy()
    nf = asecData[7].tonumpy()
    ff = asecData[8].tonumpy()

    if keck == True:
        scale = 0.01
        sgra = [0,-1]
        size = [1024, 1024]
        x = (xAsec-sgra[0])/scale + size[0]/2.
        y = (yAsec-sgra[1])/scale + size[1]/2.
    else: # TMT
        scale = 0.004
        size = [4096, 4096]
        x = -xAsec/scale + size[0]/2.
        y = yAsec/scale + size[1]/2.

    if restrictFOV == True:
        goodPos = np.where((x > 0) & (x < size[0]) & (y > 0) & (y < size[1]))[0]
        name = [name[gg] for gg in goodPos]
        mag = mag[goodPos]
        year = year[goodPos]
        x = x[goodPos]
        y = y[goodPos]
        snr = snr[goodPos]
        corr = corr[goodPos]
        nf = nf[goodPos]
        ff = ff[goodPos]
    

    outfile = infile.split('.')[0] + '_pixels.lis'
    #outfile = infile.replace('arcsec','pixels')
    oo = open(outfile,'w')
 
    fmt = '%12s  %6.3f  %8.3f  %10.5f  %10.5f  %4.2f  %4.2f  %1i  %5.3f\n'
    for ii in range(len(name)):
        oo.write(fmt % (name[ii],mag[ii],year[ii],x[ii],y[ii],
                            snr[ii],corr[ii],nf[ii],ff[ii]))
    
    oo.close()
    



def single_frames(alnDir='align_18092012/'):
    """
    Reads in the align files that have all 10 individual frames aligned
    to one another and determines the average position and error on the mean.

    Uses align10_t.* files, which have been trimmed to include only those
    stars in all 10 frames.
    """

    # plot the input vs. output aligned positions (must have run align already)
    rootdir = root + alnDir
    f = open(rootdir + 'starfinder/align/align10_t.pos')
    #s = starset.StarSet(rootdir + 'starfinder/align/align10_t')

    for line in f.readlines():
        
        pdb.set_trace()



def make_lis_from_database(outfile='',epoch=None,pixels=False):
    """
    Reads in stars table from GC database and creates a .lis file of
    positions. Set pixels to True to put into IRIS pixel
    coordinates (4 mas pixels).

    Input:
    outfile (str):   Name of output file; give full path after the TMT root directory.
    epoch (float):   Epoch to propagate positions to based on velocities in database.
    		     If epoch=None, the database positions are written out to outfile.
    """
    # Get the input positions
    dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'

    scale = 0.004 # arcsec per pixel
    size = [4096.0, 4096.0]
    radius = size[0]/2.0*scale*np.sqrt(2.0)  # actual simulated radius

    # Create a connection to the database file
    connection = sqlite.connect(dbfile)

    # Create a cursor object
    cur = connection.cursor()
    cur.execute('SELECT name,kp,x,y,r2d,ddate from stars where r2d < ?', [radius])

    name = []
    kp = []
    x0 = []
    y0 = []
    r2d = []
    t0 = []
    fmt = '%12s  %6.3f  %8.3f  %10.5f  %10.5f  %5.2f  %4.2f  %1i  %8.3f\n'
    for row in cur:
        name = np.concatenate([name, [row[0]]])
        kp = np.concatenate([kp, [row[1]]])
        x0 = np.concatenate([x0, [row[2]]]) # arcsec
        y0 = np.concatenate([y0, [row[3]]]) # arcsec
        r2d = np.concatenate([r2d, [row[4]]]) # arcsec
        t0 = np.concatenate([t0, [row[5]]])
        vx = np.concatenate([vx, [row[6]]]) # mas/yr
        vy = np.concatenate([vy, [row[7]]]) # mas/yr

    if epoch != None:
        for ss in range(len(name)):
            dt = epoch - t0[ss]
            x = x0[ss] + (vx[ss]/1.e3 * dt)
            y = y0[ss] + (vy[ss]/1.e3 * dt)

            t0sgr = epoch
    else:
        t0sgr = 2012.0
    
   
    if pixels == True:
        x = -x/scale+size[0]/2.0
        y = y/scale+size[1]/2.0
        sgr = 2048.0
    else:
        sgr = 0.0

    # Write pixel positions to a file
    _outfile = root + outfile
    out = open(_outfile,'w')
    dum = 1.0
    # Include Sgr A* at the top
    for ii in range(len(name)):
        if ii == 0:
            out.write(fmt % ('SgrA', 16.7, t0sgr, sgr, sgr,
                                dum, dum, dum, dum))
            out.write(fmt % (name[ii],kp[ii],t0[ii],x[ii],y[ii],
                                dum, dum, dum, dum))
        else:
            out.write(fmt % (name[ii],kp[ii],t0[ii],x[ii],y[ii],
                                dum, dum, dum, dum))
    out.close()




def combine_simlist_knownlist(alnDir, simlist, knownlist, convert2pixels=True,
                              restrictFOV=True):
    """
    Combines the list of known GC stars with the simulated_stars.txt file
    created when the fainter, theorized population is included in the sims.
    The latter needs to also be converted to pixels, and given dummy star names.
    Also places 16C, 16NW, and 16CC at the top of the list.

    Input files are assumed to be in the TMT directory + alnDir

    """
    wdir = root + alnDir
    
    scale = 0.004
    size = [4096.0, 4096.0]
    radius = size[0]/2.0*scale*np.sqrt(2.0)  # actual simulated radius
    
    # Simulated stars
    sim = asciidata.open(wdir + simlist)
    xSim = sim[0].tonumpy() #arcsec
    ySim = sim[1].tonumpy() #arcsec
    magSim = sim[2].tonumpy()

    # Known stars
    obs = asciidata.open(wdir + knownlist)
    nameObs = [obs[0][ss].strip() for ss in range(obs.nrows)]
    magObs = obs[1].tonumpy()
    xObs = obs[3].tonumpy() 
    yObs = obs[4].tonumpy() 

    suffix = '_arcsec'

    if convert2pixels == True:
        # Convert to pixels (simulated stars)
        xSim = -xSim/scale+size[0]/2.
        ySim = ySim/scale+size[1]/2.

        if restrictFOV == True:
            goodPos = np.where((xSim > 0) & (xSim < size[0]) & (ySim > 0) & (ySim < size[1]))[0]
        else:
            goodPos = np.arange(len(xSim))
        xSim = xSim[goodPos]
        ySim = ySim[goodPos]
        magSim = magSim[goodPos]

        # Convert to pixels (known stars)
        xObs = -xObs/scale+size[0]/2.
        yObs = yObs/scale+size[1]/2.

        if restrictFOV == True:
            goodPos = np.where((xObs > 0) & (xObs < size[0]) & (yObs > 0) & (yObs < size[1]))[0]
        else:
            goodPos = np.arange(len(xObs))
        nameObs = [nameObs[gg] for gg in goodPos]
        xObs = xObs[goodPos]
        yObs = yObs[goodPos]
        magObs = magObs[goodPos]

        suffix = '_pixels'

    # Give all stars a common year, otherwise align may do some funny stuff
    yr = 2012.0
    dum = 1.0

    top3GC = np.array(['irs16C', 'irs16NW', 'irs16CC'])
    top3 = []

    #out = open(wdir + 'known_and_simulated_%s%s.lis' % (str(yr), suffix),'w')
    fmt = '%15s  %6.3f  %8.3f  %10.5f  %10.5f  %4.2f  %4.2f  %1i  %5.3f\n'
    #for jj in range(len(top3GC)):
    #    idx = nameObs.index(top3GC[jj])
    #    out.write(fmt % (nameObs[idx].replace('irs',''), magObs[idx], yr,
    #                     xObs[idx], yObs[idx], dum, dum, dum, dum))
    #    top3 = np.concatenate([top3, [idx]])
#
#    top3 = [int(kk) for kk in top3]
#
#    for ii in range(len(nameObs)):
#        if ii in top3:
#            continue
#        else:
#            out.write(fmt % (nameObs[ii], magObs[ii], yr, xObs[ii], yObs[ii],
#                             dum, dum, dum, dum))
#
#    for jj in range(len(xSim)):
#        nameSim = 'starsim_%s' % str(jj)
#        out.write(fmt % (nameSim, magSim[jj], yr, xSim[jj], ySim[jj],
#                         dum, dum, dum, dum))
#
#
#    out.close()

    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(size=14)
    usetexTrue()
    py.clf()
    py.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.12)
    py.figure(figsize=(7,7))
    #py.hist(magObs,bins=200,histtype='step',color='b',label='Known')
    #nn,bb,pp = py.hist(magSim,bins=200,histtype='step',color='r',label='Simulated')
    py.hist(magObs,bins=np.arange(7,26,0.3),histtype='step',color='b',label='Known')
    nn,bb,pp = py.hist(magSim,bins=np.arange(7,26,0.3),histtype='step',color='r',label='Simulated')
    py.xlabel('K')
    py.ylabel('N')
    py.legend(numpoints=1,loc=2,fancybox=True,prop=prop)
    py.axis([magObs.min()-0.2,magSim.max()+0.2,0,nn.max()+50])
    py.savefig(wdir + 'known_and_simulated_mag_0.3bins.png')
    py.savefig(wdir + 'known_and_simulated_mag_0.3bins.eps')
    usetexFalse()


    

#tmt
#cd align_05102012
#idl
#.r find_stf_tmt
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_05102012/galcenter_stars_noise_bulge_20.fits', '0.8', /makeRes, /makeStars, /makePsf
#
## IDL lines used to test centroiding of PSF:
#xref = 511.0 ; zero-based array
#yref = 511.0 ; zero-based array
#lx = xref - 5.0
#ly = yref - 5.0
#ux = xref + 5.0
#uy = yref + 5.0
#orig1 = readfits('gcpsf_60_20_-4_-4.fits')
#cen1 = centroid(orig1[lx:ux,ly:uy])
#dx = xref - lx - cen1[0]
#dy = yref - ly - cen1[1]
#new = image_shift(orig1, dx, dy)
#writefits,'tmtpsf_centroid_1_4arcsec.fits', new

#orig2 = readfits('gcpsf_60_20_8_0.fits')
#cen2 = centroid(orig2[lx:ux,ly:uy])
#
#orig3 = readfits('gcpsf_1_20_0_8.fits')
#cen3 = centroid(orig3[lx:ux,ly:uy])

#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_1_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_10_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_20_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_30_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_40_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_50_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_60_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_70_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_80_0.8_stf.lis
#calibrate_new -S 16C,16NW,16CC -f 1 -M 6 -c 12 -R -T 0.0 galcenter_stars_noise_bulge_comboKLF_90_0.8_stf.lis



#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_1.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_10.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_20.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_30.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_40.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_50.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_60.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_70.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_80.fits','0.5',/makeRes,/makeStars,/makePsf
#find_stf_tmt,'/u/leo/TMTsims/iris_code_20120130/code/align_30042013/galcenter_stars_noise_bulge_90.fits','0.5',/makeRes,/makeStars,/makePsf


