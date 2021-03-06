##################################################
#
#  Functions and variables to do everything in
#  the doit executable. 
#
##################################################

import asciidata
import os, math, random, shutil
from gcwork import starset
from gcwork import objects
from gcwork import young
import applyLocalDist
import pylab as py
import numpy as np
import numarray as na
from gcreduce import gcutil
import sc_accel_class as acc #####WILL NEED TO CHANGE THIS!!!!!
import pdb

##########
#
# Global variable setup
#
##########

# Root directory
rootDir = os.path.realpath(os.curdir)
if (rootDir.endswith('scripts')): rootDir = rootDir[:-len('scripts')]

# On import, load up epochsInfor.txt
table = asciidata.open(rootDir + '/scripts/epochsInfo.txt')

# List of columns in the table. Make and array for each one.
epoch    = [table[0][ss].strip() for ss in range(table.nrows)]
dirs     = [table[1][ss].strip() for ss in range(table.nrows)]
isAO     = [table[2][ss] == 1 for ss in range(table.nrows)]
doAlign  = [int(table[3][ss]) for ss in range(table.nrows)]
angle    = [float(table[4][ss]) for ss in range(table.nrows)]
wave     = [int(table[5][ss]) for ss in range(table.nrows)]
errScale = [float(table[6][ss]) for ss in range(table.nrows)]
numSubs  = [int(table[7][ss]) for ss in range(table.nrows)]
is06     = [int(table[8][ss]) for ss in range(table.nrows)]

numEpochs = len(epoch)

# Trim out all the non-aligned stuff... we don't do anything with it anyhow
idx = (na.where(na.array(doAlign) == 1))[0]
epoch    = [epoch[i] for i in idx]
dirs     = [dirs[i] for i in idx]
isAO     = [isAO[i] for i in idx]
doAlign  = [doAlign[i] for i in idx]
angle    = [angle[i] for i in idx]
wave     = [wave[i] for i in idx]
errScale = [errScale[i] for i in idx]
numSubs  = [numSubs[i] for i in idx]

# Reference Epoch
refEpoch = '06junlgs_kp'

# Other globals
mcorr = '0.8'  # correlation threshold for main image
scorr = '0.6'  # correlation threshold for subset images

# Set -C and -S flags for calibrate to set which stars
# get used in calibration. This puts 16C and 16SW-E as
# the top two sources.
#calSrcs = '-S 16C,16NE,S2-16,S1-23,S1-3,S1-5,S2-22,S2-5,S1-68,S0-13,S1-25'
calSrcs = '-S 16C,16NW,16CC'


# Number of subset images required in align_rms
num_subs_defs = 3
num_subs_poss = 3

# Number of epochs required for trimming
ecut_d = 7
ecut_p = 4

# When doing trim_align, keep these stars even though they would not
# normally pass the number of epochs cuts. 
keepStars_d = 'S0-102'

##########
#
# Functions
#
##########

def copy_defs():
    """
    Copy uncalibrated star lists from their gc_test location.
    /u/ghezgroup/data/gc_test/lis_stf_new_on/

    These are usually produced using find_stf_new for each epoch.
    """
    for e in range(numEpochs):
        print "*** copying ", epoch[e]

        fromdir = dirs[e] + '/combo/starfinder'
        todir = rootDir + '/lis/'

        # Main
        fromfile = '%s/mag%s_%s_stf.lis' % (fromdir, epoch[e], mcorr)
        os.system('cp %s %s' % (fromfile, todir))

        for s in range(1, numSubs[e]+1):
            fromfile = '%s/m%s_%d_%s_stf.lis' % (fromdir, epoch[e], s, scorr)
            os.system('cp %s %s' % (fromfile, todir))

    os.chdir('../scripts')


def calibrate_defs():
    """
    Run calibrate on each star list.
    """
    os.chdir(rootDir + '/lis')
    
    for e in range(numEpochs):
        print '*** calibrating ', epoch[e]

        # Speckle
        if not isAO[e]:
            cmd = 'calibrate_new %s -f 1 -M %i -c 1 -R ' % (calSrcs, wave[e])
            print 'calibrate_new %s -f 1 -M %i -c 1 -R ' % (calSrcs, wave[e])
        else:
            cmd = 'calibrate_new %s -f 1 -M %i -c 4 -R -T %.1f  ' % \
                  (calSrcs, wave[e], angle[e])
            print 'calibrate_new %s -f 1 -M %i -c 4 -R -T %.1f ' % \
                  (calSrcs, wave[e], angle[e])

        # Main Map
        main = 'mag' + epoch[e] + '_' + mcorr + '_stf.lis'
        os.system(cmd + main)

        # Sub Maps
        for s in range(1, numSubs[e]+1):
            sub = 'm%s_%d_%s_stf.lis' % (epoch[e], s, scorr)
            os.system(cmd + sub)

    os.chdir(rootDir + '/scripts')


def makeLabel():
    """
    Make restricted label.dat file where young stars are removed.
    """
    
    os.chdir(rootDir + '/source_list')
    _in = asciidata.open('ao2006setup_velocities_noAccelSources.dat')
    names = _in[0]._data
    mag = _in[1].tonumpy()
    x = _in[2].tonumpy()
    y = _in[3].tonumpy()
    xerr = _in[4].tonumpy()
    yerr = _in[5].tonumpy()
    vx = _in[6].tonumpy()
    vy = _in[7].tonumpy()
    vxerr = _in[8].tonumpy()
    vyerr = _in[9].tonumpy()
    t0 = _in[10].tonumpy()
    use = _in[11].tonumpy()
    r2d = _in[12].tonumpy()

    # Turn the 'use' column into a string
    use = [str(uu) for uu in use]

    # Do not include young stars in alignment b/c of known net rotation
    yng = young.youngStarNames()
    
    outFile = 'ao2006setup_velocities_noAccelSources_noYng.dat'
    gcutil.rmall([outFile])
    _out = open(outFile,'w')
    _out.write('%-13s %5s  ' % ('#Name', 'K'))
    _out.write('%8s %8s %7s %7s  ' % ('x', 'y', 'xerr', 'yerr'))
    _out.write('%8s %8s %8s %8s  ' % ('vx', 'vy', 'vxerr', 'vyerr'))
    _out.write('%8s %5s %6s\n' % ('t0', 'use?', 'r2d'))

    _out.write('%-13s %5s  ' % ('#()', '(mag)'))
    _out.write('%8s %8s %7s %7s  ' % ('(asec)', '(asec)', '(asec)', '(asec)'))
    _out.write('%8s %8s %8s %8s  ' %
              ('(mas/yr)', '(mas/yr)', '(mas/yr)', '(mas/yr)'))
    _out.write('%8s %5s %6s\n' % ('(year)', '()', '(asec)'))


    for i in range(len(x)):
        if (names[i] in yng):
            use[i] = '0'

        _out.write('%-13s %5.1f  ' % (names[i], mag[i]))
        _out.write('%8.4f %8.4f %7.4f %7.4f  ' %
                  (x[i], y[i], xerr[i], yerr[i]))
        _out.write('%8.3f %8.3f %8.3f %8.3f  ' %
                  (vx[i], vy[i], vxerr[i], vyerr[i]))
        _out.write('%8.3f %5s %6.3f\n' % (t0[i], use[i], r2d[i]))

    _out.close()

    os.chdir(rootDir + '/scripts')
                                                                                                        
def align_defs():
    """
    Call align on the *stf_cal.lis definite star lists.
    """
    os.chdir(rootDir + '/align')
    align_starlists('align_d', '%s_stf_cal' % (mcorr), withErrors=False)
    os.chdir(rootDir + '/scripts')

        
def make_defs_rms():
    """
    Run align_rms on each star list.
    """
    os.chdir(rootDir + '/lis')

    for e in range(numEpochs):
        print '*** align_rms on ', epoch[e]

        # Make align.list
        _alignlist = open('align.list', 'w')

        # Main -- add to list
        main = 'mag' + epoch[e] + '_' + mcorr + '_stf_cal.lis'
        _alignlist.write('%s 2\n' % (main))

        # Subs -- add to list
        for s in range(1, numSubs[e]+1):
            sub = 'm%s_%d_%s_stf_cal.lis' % (epoch[e], s, scorr)
            _alignlist.write('%s 2\n' % (sub))

        _alignlist.close()

        # Call java align
        cmd = 'java align -a 0 -R 5 -p -v '
        cmd += 'align.list'
        os.system(cmd)

        # Call align_rms
        os.system('align_rms -m align %d %d' % (numSubs[e], numSubs[e]))

        # Copy over output to new starlist BUT modify errors first.
        stars = asciidata.open('align_rms.out')

        # Do some cleanup of formatting
        stars[0].reformat('%-13s') # name
        stars[1].reformat('%6.3f') # mag
        stars[2].reformat('%8.3f') # epoch
        stars[3].reformat('%8.3f') # x
        stars[4].reformat('%8.3f') # y
        stars[5].reformat('%5.3f') # xerr
        stars[6].reformat('%5.3f') # yerr
        stars[7].reformat('%7.2f') # snr
        stars[8].reformat('%5.2f') # corr
        stars[9].reformat('%5d')   # numFrames
        stars[10].reformat('%8.3f') # ??

        # Replace with new errors (quad-summed) if needed
        if (errScale[e] != 0):
            err2 = errScale[e]**2
            
            for s in range(stars.nrows):
                stars[5][s] = '%5.3f' % py.sqrt(stars[5][s]**2 + err2)
                stars[6][s] = '%5.3f' % py.sqrt(stars[6][s]**2 + err2)

        stars.writeto('mag%s_rms.lis' % (epoch[e]))

    os.system('rm align*')
    os.chdir(rootDir + '/scripts')


def apply_local_dist():
    """
    Applies local distortion correction to the non-2006 LGSAO setups 
    Creates a new modified starlist with the name mag<epoch>_<filt>_rms.lis
    The original rms.lis file is moved to mag<epoch>_<filt>_rms.orig.lis
    """
    non06 = ['04jullgs', '05junlgs', '05jullgs']

    for ee in range(len(non06)):
        applyLocalDist.go(non06[ee], 'kp', rootDir)
 
    

def align_defs_rms():
    """
    Call align on the *rms.lis definite star lists.
    Don't forget that the 4 non-2006 epochs have slightly
    different filenames after local distortion correction applied.
    """
    os.chdir(rootDir + '/align')
    align_starlists('align_d_rms', 'rms')
    os.chdir(rootDir + '/scripts')


def trim_defs_rms():
    os.chdir(rootDir + '/align')

    cmd = 'java -Xmx512m trim_align '
    cmd += '-r align_d_rms_abs_t -e %s -keep %s ' % (ecut_d, keepStars_d)
    cmd += '-p -f ../points_d/ align_d_rms_abs'

    os.system(cmd)

    os.chdir(rootDir + '/scripts')


def polyfit_defs_rms():
    cmd = 'polyfit -d 2 -linear -i ../align/align_d_rms_abs_t '
    cmd += '-tex -plot -points ../points_d/ -o ../polyfit_d/fit'

    os.system(cmd)


def align_defs_rms_1000():
    """
    Call align on the *rms.lis definite star lists.
    """
    os.chdir(rootDir + '/align')
    align_starlists('align_d_rms_1000', 'rms', bootstrap=True)
    os.chdir(rootDir + '/scripts')


def trim_defs_rms_1000():
    os.chdir(rootDir + '/align')

    cmd = 'java -Xmx512m trim_align '
    cmd += '-r align_d_rms_1000_abs_t -e %s -keep %s ' % (ecut_d, keepStars_d)
    cmd += '-p -f ../points_1000/ align_d_rms_1000_abs'

    os.system(cmd)

    os.chdir(rootDir + '/scripts')


def trim_poss_rms_1000():
    os.chdir(rootDir + '/align')

    cmd = 'java -Xmx512m trim_align '
    cmd += '-r align_p_rms_1000_abs_t -e %s -keep %s ' % (ecut_p, keepStars_d)
    cmd += '-p -f ../points_p/ align_d_rms_1000_abs'

    os.system(cmd)

    os.chdir(rootDir + '/scripts')


def polyfit_defs_rms_1000():
    cmd = 'polyfit -d 2 -linear -i ../align/align_d_rms_1000_abs_t '
    cmd += '-tex -plot -points ../points_1000/ -o ../polyfit_1000/fit'

    os.system(cmd)


def polyfit_poss_rms_1000():
    cmd = 'polyfit -d 2 -linear -i ../align/align_p_rms_1000_abs_t '
    cmd += '-tex -plot -points ../points_p/ -o ../polyfit_p/fit'

    os.system(cmd)


def align_starlists(alignRoot, listSuffix, withErrors=True, bootstrap=False):
                    
    # Make align_d_rms.lis
    _alignlist = open(alignRoot + '.list', 'w')
    for e in range(numEpochs):
        # Pick out the reference epoch
        isRef = ''
        if (epoch[e] == refEpoch):
            isRef = 'ref'

        # Map between datatype and isAO variable
        if (isAO[e] == 0 and withErrors): dataType = 3 # speckle
        if (isAO[e] == 1 and withErrors): dataType = 9 # NIRC2
        if (isAO[e] == 0 and not withErrors): dataType = 2 # speckle
        if (isAO[e] == 1 and not withErrors): dataType = 8 # NIRC2

        if ((is06[e] == 1) and (withErrors == True)):
            _alignlist.write('../lis/mag%s_%s.lis %d %s\n' % \
                             (epoch[e], listSuffix, dataType, isRef))
        elif ((is06[e] == 0) and (withErrors == True)):
            _alignlist.write('../lis/mag%s_%s_ld.lis %d %s\n' % \
                             (epoch[e], listSuffix, dataType, isRef))
        else:
            _alignlist.write('../lis/mag%s_%s.lis %d %s\n' % \
                             (epoch[e], listSuffix, dataType, isRef))
            
    _alignlist.close()

    # Align
    print '*** align ' + alignRoot

    alignType = 4    # second order transformation
    
    cmd = 'java align '
    cmd += '-a %d -r %s -d 1.0 -p -s 1.0 -i -restrict -inputVelAlign ' % (alignType, alignRoot)
    if (bootstrap == True): cmd += '-n 100 '
    cmd += '-N ../source_list/label_restrict_noAccel.dat '
    cmd += '-o ../source_list/orbits.dat -O 2 '
    cmd += alignRoot + '.list'

    os.system(cmd)


    # Align Absolute
    print '*** align absolute ' + alignRoot

    cmd = 'java align_absolute '
    cmd += '-abs ../source_list/absolute_refs_noAccel.dat '
    cmd += alignRoot + ' ' + alignRoot + '_abs'

    os.system(cmd)


def remove_speckle_edge():
    """
    Removes speckle epochs for stars that are only in a fraction
    of all the speckle frames in a given epoch. This should help
    remove problem stars on the edge of speckle field of view.
    Runs polyfit afterward.

    """
    os.chdir(rootDir)

    oldPts = rootDir + 'points_1000/'
    newPts = rootDir + 'points_s/'
    newPoly = rootDir + 'polyfit_s/'

    if os.path.isdir(newPts):
        print 'removing exisiting points_s directory: ' + newPts
        shutil.rmtree(newPts)

    if os.path.isdir(newPoly):
        print 'removing exisiting polyfit_s directory: ' + newPoly
        shutil.rmtree(newPoly)

    gcutil.mkdir('points_s')
    gcutil.mkdir('polyfit_s')

    # copy the files into the new directory
    print 'copying files to :', newPts
    cmdStr = 'cp -rf %s %s' % (oldPts, newPts)
    print 'entering command: '+cmdStr
    os.system(cmdStr)
    print 'done copying files'

    s = starset.StarSet(rootDir + 'align/align_d_rms_1000_abs_t')
    s.loadPolyfit(rootDir + 'polyfit_1000/fit', accel=0, arcsec=0)
    s.loadPolyfit(rootDir + 'polyfit_1000/fit', accel=1, arcsec=0)

    names = s.getArray('name')
    r2d = s.getArray('r2d')
    mag = s.getArray('mag')
    bnames = np.zeros((numEpochs,len(names)))
    years = s.getArray('years')[0]

    ####################
    # Speckle edge removal:
    ####################
    
    # Read the align*param file to figure out fraction of speckle
    # frames a star was detected in for a given epoch
    f_par = open(rootDir + 'align/align_d_rms_1000_abs_t.param', 'r')
    tab = asciidata.open(rootDir + '/scripts/epochsInfo.txt')
    imtype = tab[2].tonumpy()

    # The first line is for 16C, which is in every speckle frame
    # Every 4th column gives the number of frames for that epoch
    par = f_par.readline()
    n16c = [par.split()[cc] for cc in range(3,len(par.split()),4)]
    n16c = np.array([float(nn) for nn in n16c])

    # Start the loop at index 1 b/c we've read in first line already
    for i in range(1, len(names)):
        star = names[i]
        #if (r2d[i] < 1.5):
        #    par = f_par.readline()
        #    continue

        par = f_par.readline()
        # If a star is in < 80% of speckle frames, remove that epoch for that star
        nframes = [par.split()[cc] for cc in range(3,len(par.split()),4)]
        nframes = np.array([float(nn) for nn in nframes])
        frac = nframes / n16c

        # Star is on edge of speckle FOV
        edge = np.where((frac < 0.8) & (frac >= 0.0) & (imtype == 0))[0]

        # The index 'edge' tells us the epoch that we want to throw out
        # We must tie that to the dates in self.allEpochs
        spec_remove = years[edge]

        if len(edge) > 0:
            # Open this star's points file and remove the bad epochs
            f_pts = asciidata.open(rootDir + 'points_s/' + star + '.points')
            f_phot = asciidata.open(rootDir + 'points_s/' + star + '.phot')

            for ee in spec_remove:
                epochs = f_pts[0].tonumpy()
                ep2remove = np.where(epochs == ee)[0]
                if len(ep2remove) > 0:
                    f_pts.delete(ep2remove[0])
                    f_phot.delete(ep2remove[0])

            f_pts.flush()
            f_phot.flush()

    os.chdir(rootDir)

    print 'now running polyfit in ' + rootDir + 'polyfit_s'
    cmd = 'polyfit -d 2 -linear -i '+ rootDir + 'align/align_d_rms_1000_abs_t'
    cmd += ' -tex -plot -points '+ rootDir + '/points_s/ -o '+ rootDir + '/polyfit_s/fit'
    os.system(cmd)

    #py.clf()
    #py.hist(frac, bins=np.arange(0,1.0,0.05), histtype='step', color='b')
    #py.xlabel('Fraction')
    #py.ylabel('N')
    #py.savefig(rootDir + 'plots/speckle_fraction_%s.png' % ep)

    os.chdir(rootDir + '/scripts')
    

def remove_confused():
    """
    Calls accel_class.py to find nearest neighbors and removes
    confused epochs. Creates new points files in points_c/
    and runs polyfit again on these points files.
    """
    os.chdir(rootDir)

    gcutil.mkdir('points_c')
    gcutil.mkdir('polyfit_c')

    data = acc.accelClass(rootDir=rootDir,align='align/align_d_rms_1000_abs_t',
                          poly='polyfit_s/fit',points='points_s/')
    data.findNearestNeighbors()
    data.removeConfusedEpochs(mkPointsFiles=True, runPolyfit=True)
    
    os.chdir(rootDir + '/scripts')

def manual_confusion_fix():
    """
    S1-8 and S1-14 are clearly confused by something, but
    the automatic confusion removal (above) does not catch these
    cases.  Used compare_pos to determine which epochs should be
    removed. Saves a copy of the original points files.
    Then re-runs polyfit one last time.
    """

    stars = ['S1-8', 'S1-14']
    # S1-8 is confused by S0-36 during:
    # 	-- 3 epochs in 1998 and one in 1999,
    #      when S0-36 is not detected
    s18_conf = [1998.251, 1998.366, 1998.771, 1999.333]
    # S1-14 is confused by 28star_345 during:
    #	-- all 2006-2008 epochs
    s114_conf = [ 2006.336, 2006.470, 2006.541, 2007.374,
                  2007.382, 2007.612, 2008.371, 2008.562]

    ptsDir = rootdir + 'points_c/'

    for ss in stars:
        pts = ptDir + ss + '.points'
        phot = ptDir + ss + '.phot'

        # Open up the latest points file and remove confused epochs
        # First copy points file
        print 'copying points files for :', ss
        cmdStr = 'cp -rf %s/%s.points %s' % (ptsDir, ss, newDir)
        print 'entering command: '+cmdStr
        os.system(cmdStr)
        print 'done copying files'
    


def new_zero():
    """
    Calls accel_class.py to update zero point taken from 
    efit of S0-2. Creates new points files in points_nz/
    and runs polyfit again on these points files.
    """
    os.chdir(rootDir)

    gcutil.mkdir('points_nz')
    gcutil.mkdir('polyfit_nz')

    # WILL NEED TO CHANGE WHICH ACCEL CLASS IS CALLED !!!!
    data = acc.accelClass(rootDir=rootDir,align='align/align_d_rms_1000_abs_t',
                          poly='polyfit_c/fit',points='points_c/')
    data.newZero(mkPointsFiles=False, runPolyfit=False)
    
    os.chdir(rootDir + '/scripts')



def plotStar(star, align = 'align/align_d_rms_1000_abs_t', poly='polyfit_1000/fit',
             points='points_1000/',fit='linear',suffix='',LGSonly=False):
    """
    Plot positions and best fit velocity and/or acceleration for star. Possible
    arguments for fit are: linear, accel, both.
    """


    pointsFile = rootDir + points + star + '.points'
    if fit != None:
        if fit == 'linear':
            fitFile = rootDir + poly + '.' + star + '.lfit'
        elif fit == 'accel':
            fitFile = rootDir + poly + '.' + star + '.pfit'
        elif fit =='both':
            fitFile1 = rootDir + poly + '.' + star + '.lfit'
            fitFile2 = rootDir + poly + '.' + star + '.pfit'

    _tabPoints = asciidata.open(pointsFile)
    date = _tabPoints[0].tonumpy()
    x = _tabPoints[1].tonumpy() * -1.0
    y = _tabPoints[2].tonumpy()
    xerr = _tabPoints[3].tonumpy()
    yerr = _tabPoints[4].tonumpy()

    if ((fit != 'both') and (fit != None)):
        _fitPoints = asciidata.open(fitFile)
        date_f = _fitPoints[0].tonumpy()
        x_f = _fitPoints[1].tonumpy() * -1.0
        y_f = _fitPoints[2].tonumpy()
        xerr_f = _fitPoints[3].tonumpy()
        yerr_f = _fitPoints[4].tonumpy()
        py.close(2)
        py.figure(2, figsize=(6.2,5.4))
        py.clf()
        py.plot(x_f, y_f, 'k-')
        py.plot(x_f + xerr_f, y_f + yerr_f, 'k--')
        py.plot(x_f - xerr_f, y_f - yerr_f, 'k--')
    elif (fit == 'both'):
        _fitPoints = asciidata.open(fitFile1)
        date_f1 = _fitPoints[0].tonumpy()
        x_f1 = _fitPoints[1].tonumpy() * -1.0
        y_f1 = _fitPoints[2].tonumpy()
        xerr_f1 = _fitPoints[3].tonumpy()
        yerr_f1 = _fitPoints[4].tonumpy()
        py.close(2)
        py.figure(2, figsize=(6.2,5.4))
        py.clf()
        py.plot(x_f1, y_f1, 'k-')
        py.plot(x_f1 + xerr_f1, y_f1 + yerr_f1, 'k--')
        py.plot(x_f1 - xerr_f1, y_f1 - yerr_f1, 'k--')

        _fitPoints = asciidata.open(fitFile2)
        date_f2 = _fitPoints[0].tonumpy()
        x_f2 = _fitPoints[1].tonumpy() * -1.0
        y_f2 = _fitPoints[2].tonumpy()
        xerr_f2 = _fitPoints[3].tonumpy()
        yerr_f2 = _fitPoints[4].tonumpy()
        py.plot(x_f2, y_f2, 'k-')
        py.plot(x_f2 + xerr_f2, y_f2 + yerr_f2, 'k--')
        py.plot(x_f2 - xerr_f2, y_f2 - yerr_f2, 'k--')

    # Find range of plot
    halfRange = max([ abs(x.max() - x.min()), abs(y.max() - y.min()) ]) / 2.0

    padd = 0.0005
    xmax = x.min() + ((x.max() - x.min())/2.0) + halfRange + padd
    ymax = y.min() + ((y.max() - y.min())/2.0) + halfRange + padd
    xmin = x.min() + ((x.max() - x.min())/2.0) - halfRange - 2*padd
    ymin = y.min() + ((y.max() - y.min())/2.0) - halfRange - padd

    if LGSonly:
        legend_items = ['2006', '2007', '2008', '2009', '2010']
        legend_colors = ['gold','lightsalmon','orange','darkorange',
                         'orangered', 'red', 'maroon']
    else:
        legend_items = ['1995', '1996', '1997', '1998', '1999',
                        '2000', '2001', '2002', '2003', '2004',
                        '2005', '2006', '2007', '2008', '2009',
                        '2010']
        legend_colors = ['olive', 'brown', 'purple', 'darkslateblue',
                         'mediumblue', 'steelblue', 'teal', 'green',
                         'greenyellow', 'gold', 'lightsalmon', 'orange',
                         'darkorange', 'orangered', 'red', 'maroon']

    # Assign color to each epoch
    year = [str(int(na.floor( d ))) for d in date]
    color_arr = []


    for i in range(len(year)):
        # find out which year
        try:
            idx = legend_items.index( year[i] )
            color_arr.append( legend_colors[idx] )
        except ValueError:
            color_arr.append('black')

        py.errorbar([x[i]], [y[i]], fmt='ko', xerr=xerr[i], yerr=yerr[i],
                    ms=5,
                    color=color_arr[i], mec=color_arr[i], mfc=color_arr[i])

    # Set axis ranges
    py.title(star)
    ax = py.axis([xmax+0.005, xmin-0.005, ymin-0.005, ymax+0.005])
    py.gca().set_aspect('equal')
    py.xlabel('X (arcsec)')
    py.ylabel('Y (arcsec)')

    # Draw legend
#    py.legend(legend_items, numpoints=1, loc=4)
#    ltext = py.gca().get_legend().get_texts()
#    lgdLines = py.gca().get_legend().get_lines()
#    for j in range(len(ltext)):
#        py.setp(ltext[j], color=legend_colors[j])
#        lgdLines[j].set_marker('o')
#        lgdLines[j].set_mec(legend_colors[j])
#        lgdLines[j].set_mfc(legend_colors[j])
#        lgdLines[j].set_ms(5)


    # Retrieve information on the acceleration fit from the .accelFormal file
    fitFile = rootDir + poly + '.accelFormal'
    _fit = open(fitFile, 'r')

    _vel = asciidata.open(rootDir + align + '.vel')
    name = _vel[0]._data
    vcnt = _vel[1].tonumpy()
    idx = name.index(star)
    
    for line in _fit:
        entries = line.split()

        if (entries[0] == star):
            chiSqX = float( entries[7] )
            chiSqY = float( entries[15] )
            chiSqXr = chiSqX / (vcnt[idx] - 3.0)
            chiSqYr = chiSqY / (vcnt[idx] - 3.0)
            break

    print 'Reduced Chi^2 (X): %6.3f' % chiSqXr
    print 'Reduced Chi^2 (Y): %6.3f' % chiSqYr

    _fit.close()

    py.savefig(rootDir+'plots/plotStar_'+star+'_orbit'+suffix+'.png')    
    py.close(2)

