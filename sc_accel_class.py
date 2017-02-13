import os
import shutil
import sys
import math
import numpy as np
import pylab as py
from pylab import *
import objects
import util
import gcutil
from gcwork import starset
from gcwork import objects
from gcwork.polyfit import accel
from scipy import stats
import nmpfit_sy2 as nmpfit_sy
from realmedian import *
from scipy import stats
import random
import asciidata
import cPickle
import pdb

class accelClass(starset.StarSet):
    """Object containing acceleration analysis. This class will
    inherit methods from the StarSet class 

    CLASS METHODS:
        computeSpeckleAOVel(): compute the speckle and AO velocities individually and
        check the difference in the velocity vectors.

        epochDistance() - find the range in epochs when two stars are within some
        threshold distance. This function is helpful for trimming epochs when sources
        may be confused
        
        findMisMatch(starName = 'irs16NE', bootStrap=False) - check to
        see whether there are points that are mismatched in the star
        fit.

        findNearestNeighbors() - find the nearest neighbors
        
        findNearestStar() - find the minimum distance to every star
        
        findNonPhysical() - find non-physical accelerations

        fitAccel(time, x, xerr) - fit for acceleration and returns mpfit
        object.
        
        fitVelocity(time, x, xerr) - fit for a velocity and returns an mpfit
        object.

        plotchi2(removeConfused = False) - plot the chi2 distribution
        of the acceleration and velocity fits.
        
        removeConfusedEpochs(mkPointsFiles=False, runPolyfit=False) -
        make a copy of the points files into points_c/ and then trim
        out all epochs from each confused star. Can be set to produce
        new points files and run polyfit on those

        testChiSqFitXY - look what additive error is necessary to add to
        the x and y coordinates
        
        updateAccel() - update the acceleration calculations of radial
        and tangential acceleration

        newZero(mkPointsFiles=False, runPolyfit=False) -
        make a copy of the points files into points_nz/ and then updates
        zero in x and y positions and x an y errors according to zero in
        efit of S0-2. Can be set to produce new points files and run
        polyfit on those
        
        
    HISTORY: 2009-03-11 - T. Do
             2010-04-05 - T. Do - added function to compute the additive
             error necessary to get the Chi2 to be as expected.
    """

    # Notes X positions: the positions from polyfit are 'correct'
    # in the sense of RA increasing to the east, however, the
    # positions from starSet are NOT in the right direction. Must
    # multiply any starSet positions by -1.0 to get into the right
    # orientation.  Accelerations computed by polyfit are also in the
    # correct direction.

    # initialize physical and non-physical acceleration variables
    # these will be updated with real names when nonPhysical() is run
    #nonPhysical = ''
    #physical = ''
    nonPhysical = ''  # names of stars with sig. non-physical accel
    physical = ''   # names of stars with sig. physical accelerations
    nearest = np.zeros(1)   # the nearest star to every star
    nearestNonPhysical = np.zeros(1)  # the nearest star to unphysical stars
    bad = np.zeros(1)  # indices of non-physical stars
    goodAccel = np.zeros(1) # indices of sig. physical acceleration
    chiThreshold = 0     # reduced chi-square threshold

    prefix = '' # default prefix for file names
    maxEpochInd =  [] # index of stars detected in all epochs
    speckInd = [] # indices of speckle epochs
    aoInd = [] # indicies of AO epochs
    images = None # images corresponding to each aligned epoch

    rmax = 4.0 # maximum radius to compute the chi-square values for additive error tests
    rmin = 0.0 # minium radius to compute the chi-square values for additive error tests
    magmax = 14.5  # set the maximum magnitude, faintest end
    magmin = 0.0   # set the minimum magnitude, brightest end

    addErrStarsInd = 0 # the stars that are used in the additive error calculations
    
    # a python dictionary type that will hold the closest stars. The
    # keys will be the star names and it will contain an array of
    # names of closest stars. In addition, there are keys with
    # starname_dist which contain the corresponding closest approach
    # of that neighboring star. starname_ind will contain the index
    neighbors = None
    confusionThreshold = 0.06  # distance within which a star is said to be confused
    confusionMagThreshold = 5  # < delta magnitude required to be considered confused
    confusedSources = [] # names of confused sources
    confusedSourcesInd = [] # index corresponding to the confused sources
    
    # F test probabilities for justifying the use of acceleration for fits
    xFProb = []
    yFProb = []
#    pdb.set_trace()
    
    # rootDir='/u/syelda/research/gc/aligndir/06setup_only/10_03_02/'
    # rootDir='/u/tdo/nirc2/aligndir/10_03_02/',
    def __init__(self, rootDir='/u/syelda/research/gc/aligndir/11_07_05/',
                 align='align/align_d_rms_1000_abs_t',
                 poly='polyfit_1000/fit', points='points_1000/', sigma=5,
                 namesFile = False, findAddErr = False, rmax=0.0,
                 verbose=True):

        # Load the align files into the class so other methods can use
        # the info.
        self.rootDir = rootDir
        self.align = align
        self.points = points
        self.poly = poly
        self.sigma = sigma

        # prefix for plots
        #self.plotPrefix = 'plots/'+os.path.split(os.path.dirname(rootDir))[1]+'_'
        self.plotPrefix = '/u/schappell/plots/'
        self.prefix = os.path.split(os.path.dirname(rootDir))[1]
        # Load GC constants
        cc = objects.Constants()
        self.cc = cc
        
        # Load up positional information from align. 
        s = starset.StarSet(rootDir + align)
        s.loadPolyfit(rootDir + poly, arcsec=1)
        s.loadPolyfit(rootDir + poly, arcsec=1, accel=1)

        # load the points files corresponding to the polyfit. Note
        # that this points file might not be the same as the align
        # file if some trimming of confused epochs is done.
        s.loadPoints(rootDir+points)
        
        names = s.getArray('name')
        mag = s.getArray('mag')*1.0

        # load an alternative align with names. This align is to
        # provide names for the larger align and should be aligned to
        # the reference epoch. DOESN'T WORK RIGHT NOW - trimmed list
        # different than original reference list.
        if namesFile:
            print 'getting names from a different align: ', namesFile
            nameAlign = starset.StarSet(namesFile)
            refNames = nameAlign.getArray('name')

            print refNames[20:30]
            print names[20:30]
            # take the names here and replace the names in the input
            # alignment
            names[0:len(refNames)-1] = refNames

        # In mas/yr^2
        x = s.getArray('fitXa.p')
        y = s.getArray('fitYa.p')
        xerr = s.getArray('fitXa.perr')
        yerr = s.getArray('fitYa.perr')
                        
        vx = s.getArray('fitXa.v')
        vy = s.getArray('fitYa.v')
        vxe = s.getArray('fitXa.verr')
        vye = s.getArray('fitYa.verr')
        ax = s.getArray('fitXa.a')
        ay = s.getArray('fitYa.a')
        axe = s.getArray('fitXa.aerr')
        aye = s.getArray('fitYa.aerr')
        r2d = s.getArray('r2d')

        #    chi2x = s.getArray('fitXa.chi2red')
        #    chi2y = s.getArray('fitYa.chi2red')

        chi2x = s.getArray('fitXa.chi2')
        chi2y = s.getArray('fitYa.chi2')

        chi2xv = s.getArray('fitXv.chi2')
        chi2yv = s.getArray('fitYv.chi2')

        self.nEpochs = s.getArray('velCnt')

        # set the chi-sq threshold to be three times the DOF for velocity fits
        self.chiThreshold = (np.max(self.nEpochs)- 2.0)*3.0
        # index of stars with the maximum number of epochs detected
        self.maxEpochInd = np.where(self.nEpochs == np.max(self.nEpochs))[0] 
        
        # T0 for each of the acceleration fits
        self.epoch = s.getArray('fitXa.t0')

        # All epochs sampled
        self.allEpochs = np.array(s.stars[0].years)


        self.rmax = rmax  # maximum radius for chi-sq in additive error
        
        # figure out which  epochs are speckle and AO
        epochFile = rootDir+'scripts/epochsInfo.txt'

        if os.path.isfile(epochFile):
            epochTab = asciidata.open(epochFile)
            aoFlag = np.array(epochTab[2])
            alignFlag = np.array(epochTab[3])
            epochList = np.array(epochTab[0])
            epochDir = np.array(epochTab[1])
            
            # remove points not used in the align
            bad = where(alignFlag == 0)[0]
            aoFlag = np.delete(aoFlag,bad)
            epochList = np.delete(epochList,bad)
            epochDir = np.delete(epochDir,bad)
            
            # separate the speckle points from the AO points
            speckInd = np.where((aoFlag == 0))[0]
            aoInd = np.where((aoFlag == 1))[0]

            if verbose == True:
                print 'speckle epochs: ', self.allEpochs[speckInd]
                print 'AO epochs: ', self.allEpochs[aoInd]
            self.speckInd = speckInd
            self.aoInd = aoInd

            # get the names of the image files corresponding to the
            # different epochs
            
            # self.images = epochDir+'combo/'+epochList+'.fits'
        else:
            print 'epochsInfo file missing: '+epochFile


        self.names = np.array(names)  # star names
        self.mag = mag  # magnitudes
        self.cc = cc    # constants
        self.x, self.y = x, y  # fited positions at T0
        self.xerr, self.yerr = xerr, yerr # positional errors
        self.vx, self.vy = vx, vy  # velocities
        self.vxe, self.vye = vxe, vye  # velocity errors
        self.ax, self.ay = ax, ay  # accelerations
        self.axe, self.aye = axe, aye  # acceleration errors

        self.chi2x, self.chi2y = chi2x, chi2y   # acceleartion fit chi-square
        self.chi2xv, self.chi2yv = chi2xv, chi2yv # velocity fit chi-square
        
        self.r2d = r2d    # radial distance from Sgr A*
        self.starSet = s    # the original starset object

        self.nearestStarName = zeros(len(x), 'str')
        self.nearestStarDist = zeros(len(x))

        # scale factor for the errors in computing additive factor
        self.errScaleFactor = 1.0

        # initialize the acceleration information
        self.updateAccel()
        self.computeFTest()  # compute the F test 
        self.computeJerk()  # compute the jerk -- not functional yet
        self.computeFTestJerk()  # compute the F test 
        self.findNearestStar()
        self.findNonPhysical(verbose=verbose)
        self.findNearestNeighbors()
        self.removeConfusedEpochs()
        self.newZero()
        
        # self.findMismatch(bootStrap=True)
        #if (len(self.speckInd) > 3):
        #    self.computeSpeckleAOVel(requireAllEpochs=True)
        #self.saveClass()
        #self.testChi2Fit()

        if len(self.speckInd) > 0:
            data = 'speckle'
        else:
            data = 'ao'
        if findAddErr:
            # testErr should be in arcseconds
            self.testChi2Fit(scaleFactor = self.errScaleFactor, data = data,
                             testErr = np.arange(0.00005,0.0005,0.00001))

    def updateAccel(self):
        """ Updates the radial and tangential acceleration arrays.
        """
        cc = self.cc
        x, y = self.x, self.y
        vx, vy = self.vx, self.vy
        ax, ay = self.ax, self.ay
        axe, aye = self.axe, self.aye
        
        # Lets do radial/tangential
        r = np.sqrt(x**2 + y**2) 
        if ('Radial' in self.poly):
            at = ax
            ar = ay
            ate = axe
            are = aye
        else:
            ar = ((ax*x) + (ay*y)) / r
            at = ((ax*y) - (ay*x)) / r
            are = np.sqrt((axe*x)**2 + (aye*y)**2) / r
            ate = np.sqrt((axe*y)**2 + (aye*x)**2) / r

        # Total acceleration
        atot = py.hypot(ax, ay)
        atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

        # Calculate the acceleration limit set by the projected radius
        # Convert into cm
        r2d = r * cc.dist * cc.cm_in_au
        rarc =np.arange(0.01,10.0,0.01)
        rsim = rarc * cc.dist * cc.cm_in_au

        # acc1 in cm/s^2
        a2d = -cc.G * cc.mass * cc.msun / r2d**2
        a2dsim = -cc.G * cc.mass * cc.msun / rsim**2
        # acc1 in km/s/yr

        a2d *= cc.sec_in_yr / 1.0e5
        # a2d *= 1000.0 / cc.asy_to_kms

        # convert between arcsec/yr^2 to km/s/yr
        ar *= cc.asy_to_kms
        are *= cc.asy_to_kms
        
        at *= cc.asy_to_kms
        ate *= cc.asy_to_kms

        a2dsim *= cc.sec_in_yr / 1.0e5
        # a2dsim *= 1000.0 / cc.asy_to_kms

        self.ar = ar      # radial acceleration
        self.are = are    # radial acceleration error
        self.at = at      # tangential acceleration
        self.ate = ate    # tangential acceleration error
        self.atot = atot  # total acceleration
        self.a2d = a2d    # maximum acceleration for each star at the x, y position
        self.rarc, self.a2dsim = rarc, a2dsim # radial distance vs. accel (z = 0)

    def minDistance(self, xfit1,yfit1,xfit2,yfit2,t1,t2,trange,pause=0):
        """
        From the coefficients of the fits to two stars, find the minimum
        distance between them by using the acceleration fits. 
        xfit1 = [x0, vx, ax]
        yfit1 = [y0, vy, ay]
        
        HISTORY: 2010-01-13 - T. Do
        """
        time = np.arange(trange[0],trange[1]+1,0.01)
        time1 = time-t1
        time2 = time-t2
        xpos1 = self.poly2(xfit1, time1)
        ypos1 = self.poly2(yfit1, time1)
        xpos2 = self.poly2(xfit2, time2)
        ypos2 = self.poly2(yfit2, time2)
        
        distance = np.sqrt((xpos1-xpos2)**2 + (ypos1 - ypos2)**2)
        if pause:
            clf()
            plot(xpos1,ypos1)
            plot(xpos2,ypos2)
        #print time
        #pdb.set_trace()
        return np.amin(distance)

    def epochDistance(self, xfit1,yfit1,xfit2,yfit2,t1,t2,trange,dist,pause=False):
        """ Given two sets of acceleration fits, and a threshold distance,
        compute the range of epochs in which the two stars are within
        that distance. Returns None if the two stars never gets that close. 
        """
        
        time = np.arange(trange[0],trange[1]+1,0.01)
        time1 = time-t1
        time2 = time-t2
        xpos1 = self.poly2(xfit1, time1)
        ypos1 = self.poly2(yfit1, time1)
        xpos2 = self.poly2(xfit2, time2)
        ypos2 = self.poly2(yfit2, time2)
        
        distance = np.sqrt((xpos1-xpos2)**2 + (ypos1 - ypos2)**2)

        good = np.where(distance <= dist)[0]
        if pause:
            clf()
            plot(xpos1,ypos1)
            plot(xpos2,ypos2)
            plot(xpos1[good],ypos1[good],'bo')
            plot(xpos2[good],ypos2[good],'go')
        #print time
        #pdb.set_trace()

        if (len(good) > 0):
            return np.array([np.min(time[good]), np.max(time[good])])
        else:
            return None

    def epochDistanceTest(self, index = 0):
        """ test the epochDistance procedure. Should have run
        findNearestNeighbors() funcition.
        """
        star1 = self.confusedSources[index]
        star2 = (self.neighbors[star1])[0]
        print star1, star2
        names = list(self.names)
        s1 = names.index(star1)
        s2 = names.index(star2)

        x,y = self.x,self.y
        vx, vy = self.vx, self.vy
        ax, ay = self.ax, self.ay
        
        xfit1 = [x[s1],vx[s1],ax[s1]]
        xfit2 = [x[s2],vx[s2],ax[s2]]
        yfit1 = [y[s1],vy[s1],ay[s1]]
        yfit2 = [y[s2],vy[s2],ay[s2]]

        epoch1 = self.epoch[s1]
        epoch2 = self.epoch[s2]
        minEpoch = np.amin(self.allEpochs)
        maxEpoch = np.amax(self.allEpochs)
        epochRange = self.epochDistance(xfit1, yfit1, xfit2, yfit2, epoch1, epoch2,
                                        [minEpoch, maxEpoch],self.confusionThreshold, 
                                        pause = True)
        text(x[s1],y[s1],star1)
        text(x[s2],y[s2],star2)
        print epochRange
    
    def fitLin(self, p, fjac=None, x=None, y=None, err=None):
        # linear fit
        fun = p[0] + p[1]*x
        
        # deviations from the model
        deviates = (y - fun)/err
        return [0, deviates]

    def line(self, p, x):
        return p[0]+p[1]*x
    
    def fitfunPoly2(self, p, fjac=None, x=None, y=None, err=None):
        # second order polynomial
        fun = p[0] + p[1]*x + 0.5*p[2]*x**2
        
        # deviations from the model
        deviates = (y - fun)/err
        return [0, deviates]

    def poly2(self, p, x):
        return p[0] + p[1]*x + 0.5*p[2]*x**2


    def fitfunPoly3(self, p, fjac=None, x=None, y=None, err=None):
        # third order polynomial
        fun = p[0] + p[1]*x + 0.5*p[2]*x**2 + (1./6.)*p[3]*x**3
        
        # deviations from the model
        deviates = (y - fun)/err
        return [0, deviates]

    def poly3(self, p, x):
        return p[0] + p[1]*x + 0.5*p[2]*x**2 + (1./6.)*p[3]*x**3


    def isPhysical(self, x, y, ax, axe, ay, aye, arcsec = True):
        """
        Return a True or False depending on whether the acceleration
        measurement is physical. Unphysical accelerations are defined as
        either: 1. Significant tangential acceleration 2. Significant
        positive radial acceleration 3. Negative acceleration greater than
        the maximum allowed at the 2D position.

        RETURN: status array where 1 is true [tangential, pos. radial, > max radial]
        """

        status = np.zeros(3)  # return array
        
        # Lets do radial/tangential
        r = np.sqrt(x**2 + y**2)
        ar = ((ax*x) + (ay*y)) / r
        at = ((ax*y) - (ay*x)) / r
        are = np.sqrt((axe*x)**2 + (aye*y)**2) / r
        ate = np.sqrt((axe*y)**2 + (aye*x)**2) / r

        # Total acceleration
        atot = py.hypot(ax, ay)
        atoterr = np.sqrt((ax*axe)**2 + (ay*aye)**2) / atot

        # Calculate the acceleration limit set by the projected radius

        if arcsec:
            #convert to mks
            cc = objects.Constants()

            # Convert into cm
            r2d = r * cc.dist * cc.cm_in_au

            rarc =np.arange(0.01,10.0,0.01)
            rsim = rarc * cc.dist * cc.cm_in_au

            # acc1 in cm/s^2
            a2d = -cc.G * cc.mass * cc.msun / r2d**2
            a2dsim = -cc.G * cc.mass * cc.msun / rsim**2

            # acc1 in km/s/yr

            a2d *= cc.sec_in_yr / 1.0e5

            #a2d *= 1000.0 / cc.asy_to_kms

            # convert between arcsec/yr^2 to km/s/yr
            ar *= cc.asy_to_kms
            are *= cc.asy_to_kms

            at *= cc.asy_to_kms
            ate *= cc.asy_to_kms

            a2dsim *= cc.sec_in_yr / 1.0e5
            #a2dsim *= 1000.0 / cc.asy_to_kms

        # tests to see if the accelerations are physical
        if (abs(at/ate) > self.sigma):
            #print 'significant tangential acceleration'
            status[0] = 1

        #print 'radial acceleration %f +- %f' % (ar, are)
        #print 'tangential acceleration %f +- %f' % (at, ate)
        if ((ar - (self.sigma*are)) > 0):
            #print 'positive radial acceleration'
            status[1] = 1

        if (ar + (self.sigma*are) < a2d):
            #print 'too large radial acceleration'
            status[2] = 1

        # print some diagnostic info
        print 'isPhsyical: ar, are, at, ate: ', ar, are, at, ate
        return status


    def findNearestStar(self):
        """Use the acceleration fits to find the nearest stars that
        gets within a certain distance. Will store in the dictionary
        self.neighbors
        """
        
        # get some variables in to local scope
        names = np.array(self.names)
        sigma = self.sigma
        x, y = self.x, self.y
        vx, vy = self.vx, self.vy
        ax, ay = self.ax, self.ay
        mag = self.mag
        cc = self.cc
        chi2x, chi2y = self.chi2x, self.chi2y
        chi2xv, chi2yv = self.chi2xv, self.chi2yv
        r2d = self.r2d
        rarc, a2dsim = self.rarc, self.a2dsim
        epoch, allEpochs = self.epoch, self.allEpochs
        
        nearestStar = np.zeros(len(x))  # distance to nearest star
        nearestStarInd = np.zeros(len(x),'int32') # index of nearest star
        
        minEpoch = np.amin(allEpochs)
        maxEpoch = np.amax(allEpochs)

        # look for the minium distance between all stars according to the fit
        for ii in np.arange(0,len(x)):
            distances = sqrt((x[ii] - x)**2 + (y[ii] - y)**2)
            srt = distances.argsort()
            distances = distances[srt]
            #nearestStar[ii] = distances[1]

            #loop over the 5 closest sources to see which one is closest
            #print names[ii]
            fitMin = np.zeros(5)
            
            for rr in np.arange(1,6):
                #print names[srt[rr]]

                xfit1 = [x[ii],vx[ii],ax[ii]]
                xfit2 = [x[srt[rr]],vx[srt[rr]],ax[srt[rr]]]
                yfit1 = [y[ii],vy[ii],ay[ii]]
                yfit2 = [y[srt[rr]],vy[srt[rr]],ay[srt[rr]]]
                fitMin[rr-1] = self.minDistance(xfit1,yfit1,xfit2,yfit2,
                                                epoch[ii],epoch[srt[rr]],[minEpoch,maxEpoch])

            minInd = np.argmin(fitMin)
            nearestStar[ii] = fitMin[minInd]
            # add 1 to the index because the first star will be itself
            nearestStarInd[ii] = srt[minInd+1] 
            
        self.nearestStarName = names[nearestStarInd]
        self.nearestStarDist = nearestStar

    def findNearestNeighbors(self):
        """Use the acceleration fits to find the nearest stars that
        gets within a certain distance. Will store in the dictionary
        self.neighbors
        """
        threshold = self.confusionThreshold
        magThreshold = self.confusionMagThreshold
        # get some variables in to local scope
        names = np.array(self.names)
        sigma = self.sigma
        x, y = self.x, self.y
        vx, vy = self.vx, self.vy
        ax, ay = self.ax, self.ay
        mag = self.mag
        cc = self.cc
        chi2x, chi2y = self.chi2x, self.chi2y
        chi2xv, chi2yv = self.chi2xv, self.chi2yv
        r2d = self.r2d
        rarc, a2dsim = self.rarc, self.a2dsim
        epoch, allEpochs = self.epoch, self.allEpochs

        # reset the confusion arrays
        self.neighbors = None
        self.confusedSources = []
        self.confusedSourcesInd = []
        
        nearestStar = np.zeros(len(x))  # distance to nearest star
        nearestStarInd = np.zeros(len(x),'int32') # index of nearest star
        
        minEpoch = np.amin(allEpochs)
        maxEpoch = np.amax(allEpochs)

        # look for the minium distance between all stars according to the fit
        for ii in np.arange(0,len(x)):
            distances = sqrt((x[ii] - x)**2 + (y[ii] - y)**2)
            srt = distances.argsort()
            distances = distances[srt]
            closeStarNames = names[srt[0:10]]
            closeStarInds = srt[0:10]
            # find the differences in magnitude between primary and neighbors
            magDiff = mag[closeStarInds] - mag[ii]
            #nearestStar[ii] = distances[1]

            #loop over the 6 closest sources to see which one is closest

            fitMin = np.zeros(len(closeStarNames))
            
            for rr in np.arange(1,len(closeStarNames)):
                #print names[srt[rr]]
                xfit1 = [x[ii],vx[ii],ax[ii]]
                xfit2 = [x[srt[rr]],vx[srt[rr]],ax[srt[rr]]]
                yfit1 = [y[ii],vy[ii],ay[ii]]
                yfit2 = [y[srt[rr]],vy[srt[rr]],ay[srt[rr]]]

                fitMin[rr] = self.minDistance(xfit1,yfit1,xfit2,yfit2,
                                             epoch[ii],epoch[srt[rr]],[minEpoch,maxEpoch])

                
            # find all the stars that will come within the threshold
            # (and exclude the star itself)
            
            # include only stars that are brighter than the magnitude
            # threshold difference. Note that the magnitude filter is
            # non-symmetric. If magThreshold = 5, A 10th magnitude
            # star will not be considered confused next to a 15th
            # magnitude star, but the 15th magnitude star will be
            # deemd confused by the 10th magnitude star
            
            good = np.where((fitMin < threshold) & (fitMin > 0) & (magDiff < magThreshold))[0]
            if (len(good) > 0):
                # record the name of the confused source
                self.confusedSources = np.concatenate((self.confusedSources, [names[ii]]))
                self.confusedSourcesInd = np.concatenate((self.confusedSourcesInd, [ii]))
                
                minInd = np.argmin(fitMin)

                    
                # fill in the arrays
                if self.neighbors is None:
                    self.neighbors = {names[ii]:closeStarNames[good],
                                      names[ii]+'_dist':fitMin[good],
                                      names[ii]+'_ind': closeStarInds[good]}
                else:
                    self.neighbors[names[ii]] = closeStarNames[good]
                    self.neighbors[names[ii]+'_dist'] = fitMin[good]
                    self.neighbors[names[ii]+'_ind'] = closeStarInds[good]
            else:
                if self.neighbors is None:
                    self.neighbors = {names[ii]:None,
                                      names[ii]+'_dist':None,
                                      names[ii]+'_ind':None}
                else:
                    self.neighbors[names[ii]] = None
                    self.neighbors[names[ii]+'_dist'] = None
                    self.neighbors[names[ii]+'_ind'] = None

    def removeConfusedEpochs(self, mkPointsFiles = False, debug = False, runPolyfit=False):
        """ Take a look at the closest neighbors of all stars to
        figure out what epochs are possibly confused. The epochs that
        are found in one star but not the other are assumed to be
        confused and should be removed. This function requires that
        findNearestNeighbors() has been run.

        Will create new points files in points_c/ without the confused
        epochs. The file points_c/removedEpochs.txt will contain a
        list of stars, number of epochs removed, and the epochs
        removed
        
        Keywords: mkPointsFiles - create new points files that have
        the confused sources removed and put them in points_c/

        runPolyfit - run poly fit on the points_c directory and put
        polyfit results in polyfit_c/
        """
        names = list(self.names)
        x = self.x
        y = self.y
        vx, vy, ax, ay = self.vx, self.vy, self.ax, self.ay
        r2d = self.r2d
        mag = self.mag
        neighbors = self.neighbors
        epochsRemoved = zeros(len(names))
        oldDir = os.path.split(self.rootDir+self.points)[0]
        newDir = self.rootDir+'points_c/'
        epochLimits = np.array([np.min(self.allEpochs), np.max(self.allEpochs)])
        epochsRemoved = None
        nEpochsRemoved = np.zeros(len(names))
        culprits = None
        
        if mkPointsFiles:
            # make a new directory and then copy all the old points files there
            if os.path.isdir(newDir):
                print 'removing exisiting points_c directory: '+newDir
                shutil.rmtree(newDir)

            # copy the files into the new directory
            print 'copying files to :', newDir
            cmdStr = 'cp -rf %s %s' % (oldDir, newDir)
            print 'entering command: '+cmdStr
            os.system(cmdStr)
            print 'done copying files'
            workDir = newDir+'/'
        else:
            # if not making new files then use the old directory
            workDir = oldDir+'/'

        if debug:
            nRange = arange(30)
        else:
            nRange = np.arange(len(names))


        # loop through stars and go over the ones that are greater
        # than 0.5" away because those stars could potentially have a
        # lot of problems and are likely not able to be fit by
        # acceleration alone

        for i in nRange:
            if (r2d[i] > 0.5) & (neighbors[names[i]] is not None):
                star1 = names[i]
                closeStars = neighbors[names[i]]
                closeStarsInd = neighbors[names[i]+'_ind']

                # open the points files
                ptFile = workDir + names[i]+'.points'
                phFile = workDir + names[i]+'.phot'

                tab1 = asciidata.open(ptFile)
                phot1 = asciidata.open(phFile)

                if tab1.ncols > 0:
                    epochs1 = tab1[0].tonumpy()

                    # loop through the stars close to this one
                    for j in arange(len(closeStars)):
                        star2 = closeStars[j]
                        ptFile2 = workDir+closeStars[j]+'.points'
                        phFile2 = workDir+closeStars[j]+'.phot'

                        tab2 = asciidata.open(ptFile2)
                        phot2 = asciidata.open(phFile2)

                        if tab2.ncols > 0:
                            epochs2 = tab2[0].tonumpy()

                            # get the range of epochs when they were confused
                            s1 = i  
                            s2 = names.index(star2)  # index of close neighbor

                            xfit1 = [x[s1],vx[s1],ax[s1]]
                            xfit2 = [x[s2],vx[s2],ax[s2]]
                            yfit1 = [y[s1],vy[s1],ay[s1]]
                            yfit2 = [y[s2],vy[s2],ay[s2]]
                            epochRange = self.epochDistance(xfit1, yfit1, xfit2, yfit2,
                                                            self.epoch[s1], self.epoch[s2],
                                                            epochLimits, self.confusionThreshold)

                            good1 = np.where((epochs1 >= epochRange[0]) & (epochs1 <= epochRange[1]))[0]
                            good2 = np.where((epochs2 >= epochRange[0]) & (epochs2 <= epochRange[1]))[0]

                            # find the epochs that are not overlapping (the
                            # star was detected in one epoch for one of the
                            # stars, but not the other)
                            missingEpochs = np.setxor1d(epochs1[good1],epochs2[good2])


                            # record the epochs that will be removed in the star
                            missingEpochs1 = np.intersect1d(missingEpochs, epochs1)

                            if epochsRemoved is None:
                                epochsRemoved = {star1:missingEpochs1}
                            else:
                                if star1 in epochsRemoved:
                                    epochsRemoved[star1] = \
                                                np.concatenate([epochsRemoved[star1], missingEpochs1])
                                    epochsRemoved[star1] = np.unique(epochsRemoved[star1])

                                else:
                                    epochsRemoved[star1] = missingEpochs1
  
                            nEpochsRemoved[s1] = len(epochsRemoved[star1])

                            if debug:
                                print star1, self.mag[s1], star2, self.mag[s2]
                                print star1, epochs1
                                print star2, epochs2
                                print 'missing: ', missingEpochs

                            #print star1
                            #if star1 == 'S3-8':
                            #    pdb.set_trace()
                            # Keep track of the culprit star to print to text file
                            #if culprits is None:
                            #    culprits = {star1:star2}
                            #else:
                                #if star1 in culprits:
                                #    culprits[star1] = \
                                #                np.concatenate([culprits[star1], star2])
                                #    culprits[star1] = np.unique(culprits[star1])
                                #else:
                                    #culprits[star1] = star2
  
                            # look for the epochs with bad points and remove them
                            if mkPointsFiles:
                                for ee in missingEpochs:
                                    epochs1 = tab1[0].tonumpy()
                                    epochs2 = tab2[0].tonumpy()
                                    bad1 = np.where(epochs1 == ee)[0]
                                    bad2 = np.where(epochs2 == ee)[0]
                                    if debug:
                                        print 'bad1 ', bad1
                                        print 'bad2 ', bad2
                                    if len(bad1) > 0:
                                        tab1.delete(bad1[0])
                                        phot1.delete(bad1[0])
                                    if len(bad2) > 0:
                                        tab2.delete(bad2[0])
                                        phot2.delete(bad2[0])

                                # write out the points files
                                tab1.flush()
                                phot1.flush()

                                tab2.flush()
                                phot2.flush()
                            #if star1 == 'S1-2':
                            #    pdb.set_trace()
                                 
        # output a file with the number of epochs removed
        if mkPointsFiles:
            print 'writing: '+ workDir+'epochsRemoved.txt'
            output = open(workDir+'epochsRemoved.txt','w')
            for rr in np.arange(len(nEpochsRemoved)):
                if nEpochsRemoved[rr] > 0:
                    epochStr = ' '.join(np.array(epochsRemoved[names[rr]],dtype='str'))+'\n'
                else:
                    epochStr = '\n'
                outStr = names[rr]+'\t'+str(int(nEpochsRemoved[rr]))+'\t'+epochStr
                output.write(outStr)
            output.close()

        # output a file with the stars that caused another star to be confused 
        #print 'writing: '+ workDir+'confusingSources.txt'
        #output = open(workDir+'confusingSources.txt','w')
        #for cc in np.arange(len(names)):
        #    if culprits[cc] > 0:
        #        epochStr = ' '.join(np.array(culprits[rr],dtype='str'))+'\n'
        #    else:
        #        epochStr = '\n'
        #    output.write(epochStr)
        #output.close()

        if runPolyfit:
            print 'now running polyfit in '+self.rootDir+'polyfit_c'
            #gcutil.mkdir(self.rootDir+'polyfit_c')
            #cmd = 'polyfit -d 2 -linear -i '+self.rootDir+'/align/align_d_rms_1000_abs_t '
            cmd = 'polyfit -d 2 -linear -i '+self.rootDir+self.align
            cmd += ' -tex -plot -points '+self.rootDir+'/points_c/ -o '+self.rootDir+'/polyfit_c/fit'
            os.system(cmd)
                
            
    def findNonPhysical(self, plotDist = False, epochsRequired=0.0, verbose=True):
        """
        Print out a list of stars that have non-physical accelerations.
        
        Return: returns a list of star names that are unphysical 
        """

        # get some variables in to local scope
        names = np.array(self.names)
        at, ate = self.at, self.ate
        ar, are = self.ar, self.are
        sigma = self.sigma
        x, y = -self.x, self.y
        vx, vy = self.vx, self.vy
        ax, ay = self.ax, self.ay
        epoch, allEpochs = self.epoch, self.allEpochs
        nEpochs = self.nEpochs  # number of epochs each star was detected
        a2d = self.a2d
        mag = self.mag
        cc = self.cc
        chi2x, chi2y = self.chi2x, self.chi2y
        chi2xv, chi2yv = self.chi2xv, self.chi2yv
        r2d = self.r2d
        rarc, a2dsim = self.rarc, self.a2dsim
        
        ##########
        #
        # Non-physical accelerations.
        #
        ##########

        # 1. Look for signficant non-zero tangential accelerations.
        bad1 = (np.where( abs(at/ate) > sigma))[0]
        if verbose == True:
            print 'Found %d stars with tangential accelerations > %4.1f sigma' % \
                  (len(bad1), sigma)

            print names[bad1]
        # 2. Look for significant positive radial accelerations
        bad2 = (np.where( ((ar - (sigma*are)) > 0) ))[0]
        if verbose == True:
            print 'Found %d stars with positive radial accelerations > %4.1f sigma' % \
                  (len(bad2), sigma)

            print names[bad2]
        
        # 3. Look for negative radial accelerations that are higher
        # than physically allowed.
        bad3 = (np.where( ar + (sigma*are) < a2d))[0]
       
        if verbose == True:
            print 'Found %d stars with too large negative radial accelerations' % \
              (len(bad3))
            print names[bad3]
        bad = np.unique( np.concatenate((bad1, bad2, bad3)) )

        # chose only bright stars
        nbad = len(bad)
        
        bright = (np.where((mag[bad] < 16.0) & (r2d[bad] < 4.0)))[0]

        bad = bad[bright]

    #    bad = bad2
        # Sort
        rdx = r2d[bad].argsort()
        bad = bad[rdx]

        ######### look for physical accelerations
        goodAccel = np.where((ar + sigma*are < 0) & (mag < 16.0) & (r2d < 4.0) & (nEpochs > epochsRequired))[0]

        # remove the ones that also have unphysical accelerations
        goodAccel = np.setdiff1d(goodAccel,bad)

        nearestStar = self.nearestStarDist
        if verbose == True:
            print 'total number of stars detected: ', len(names)
            print '%d faint stars with non-physical accelerations' % (nbad - len(bright))
            print '*** Found %d stars with significant physical accelerations that do not have an unphysical acceleration ***' % len(goodAccel)
            print '|Name     |Mag     |x      |y    |chi2x     |chi2y      |ar     |are     |at    |ate     |nearest| X F-prob | Y F-prob|'
        
            for bb in goodAccel:
                # returnNames = np.concatenate([returnNames,[names[bb]]])
                print '|  %6s  |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f| %7.5f | %7.5f |' % (names[bb],  mag[bb], x[bb], y[bb],chi2x[bb],chi2y[bb],ar[bb],are[bb],at[bb],ate[bb],nearestStar[bb],self.xFProb[bb],self.yFProb[bb])

        returnNames = names[bad]
        if verbose == True:
            print ''
            print '*** Found %d stars with non-physical accelerations ***' % len(bad)

            print '|Name     |Mag     |x      |y    |chi2x     |chi2y      |ar     |are     |at    |ate     |nearest| X F-prob | Y F-prob|'

        
            for bb in bad:
                # returnNames = np.concatenate([returnNames,[names[bb]]])
                print '|  %6s  |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f| %7.5f | %7.5f |' % (names[bb],  mag[bb], x[bb], y[bb],chi2x[bb],chi2y[bb],ar[bb],are[bb],at[bb],ate[bb],nearestStar[bb],self.xFProb[bb],self.yFProb[bb])
    


        if plotDist:
            # plot the histograms
            clf()
            subplot(221)
            nAll, allbins, patches = hist(nearestStar,bins=50)

            nBad, badbins, patches2 = hist(nearestStar[bad],bins=10)
            setp(patches2, 'facecolor','r')
            xlabel('Nearest Neighbor (arcsec)')
            ylabel('N Stars')


            # do a ks test of the two distributions
            ksTest = stats.ks_2samp(nAll, nBad)
            print 'K-S test of the distribution of nearest star in whole sample vs. in non-physical sources'
            print ksTest

            subplot(222)
            plot(rarc,a2dsim)
            plot([0,10],[0,0],hold=True)

            #plot(r2d[bad],ar[bad],are[bb],'ro')
            py.errorbar(r2d[bad],ar[bad],are[bad],fmt='ro')
            xlim(0,5)
            ylim(np.max(ar[bad]),np.min(ar[bad]))
            xlabel('Distance from Sgr A* (arcsec)')
            ylabel('Acceleration (km/s/yr)')


            subplot(223)
            # look at the distribution of chi-squares for stars brighter than 16 magnitude
            good = np.where(mag < 16.0)

            # use mean degree of freedom:
            dof = np.floor(np.mean(self.nEpochs)-3)

            print 'Mean degree of freedom %f' % dof

            n, bins, patches1 = hist(chi2x[good], bins = np.arange(0, 20, 0.1),
                                     normed=1,label='x',alpha=0.6,color='blue')
            n, bins, patches2 = hist(chi2y[good], bins = bins, normed=1,
                                     label='y',alpha=0.6,color='green')
            xlabel('Reduced Chi^2')

            chiInd = np.arange(0.01,100,0.01)
            chi2Theory = stats.chi2.pdf(chiInd,dof)
            plot(chiInd, chi2Theory)

            legend()
            subplot(224)
            n, bins, patches1 = hist(chi2x[bad], bins = np.arange(0,20,0.5), normed=1,label='x',alpha=0.6)
            n, bins, patches2 = hist(chi2y[bad], bins = bins,normed=1,label='y',alpha=0.6)
            xlabel('Reduced Chi^2')
            plot(chiInd, chi2Theory)

            legend()

        # update list of non-physical stars
        self.nonPhysical = returnNames
        self.physical = names[goodAccel]
        self.nearest = nearestStar
        self.nearestNonPhysical = nearestStar[bad]
        self.bad = bad  # indices of bad points
        self.goodAccel = goodAccel 

    def findMismatch(self, starName = 'irs16NE', mkplots = True, bootStrap = False, iterate = False, useAlign = False):
        """Check the list of unphysical stars to find which of them
        are mismatches. This will be done by looking at fits with
        large chisquare values in the acceleration and removing points
        that are causing the deviation to see if the chi-square
        improves.

        KEYWORDS: useAlign - use the points from the align file (by
        default uses the points from the points directory)
        """

         # threshold in chi-square to check for mismatches
        chiThreshold = self.chiThreshold
        
        # find the star and get its position
        names = list(self.names)
        ind = names.index(starName)

        if useAlign:
            x = -self.starSet.stars[ind].getArrayAllEpochs('x')
            y = self.starSet.stars[ind].getArrayAllEpochs('y')
            # positional errors
            xerr_p = self.starSet.stars[ind].getArrayAllEpochs('xerr_p')
            yerr_p = self.starSet.stars[ind].getArrayAllEpochs('yerr_p')
            
            # alignment errors
            xerr_a = self.starSet.stars[ind].getArrayAllEpochs('xerr_a')
            yerr_a = self.starSet.stars[ind].getArrayAllEpochs('yerr_a')

            # add the two errors in quadrature
            xerr = np.sqrt(xerr_p**2 + xerr_a**2)
            yerr = np.sqrt(yerr_p**2 + yerr_a**2)
        else:
            x = -self.starSet.stars[ind].getArrayAllEpochs('pnt_x')
            y = self.starSet.stars[ind].getArrayAllEpochs('pnt_y')
            # positional+alignment errors
            xerr = self.starSet.stars[ind].getArrayAllEpochs('pnt_xe')
            yerr = self.starSet.stars[ind].getArrayAllEpochs('pnt_ye')

        
        t0 = self.starSet.stars[ind].fitXa.t0    # use the same T0 as polyfit
        years = np.array(self.starSet.stars[ind].years)

        good = (np.where(abs(x) < 500))[0]
        if (len(good) > 1):
            x = x[good]
            y = y[good]
            xerr = xerr[good]
            yerr = yerr[good]
            years = years[good]

        print 'x ', x
        print 'xerr ', xerr
        print 'y ', y
        print 'y err', yerr
        functargsX = {'x':years-t0, 'y':x, 'err': xerr}
        functargsY = {'x':years-t0, 'y':y, 'err': yerr}

        p0x = [self.x[ind],(np.amax(x) - np.amin(x))/(np.amax(years)-np.amin(years)), 0]
        p0y = [self.y[ind],(np.amax(y) - np.amin(y))/(np.amax(years)-np.amin(years)), 0]
        
        #clf()
        #plot(x,y,'o')
        #py.errorbar(x,y,xerr,yerr,fmt='o')
        #do the fitting
        xfit = nmpfit_sy.mpfit(self.fitfunPoly2, p0x, functkw=functargsX,quiet=1)
        yfit = nmpfit_sy.mpfit(self.fitfunPoly2, p0y, functkw=functargsY,quiet=1)
        dof = len(x) - len(xfit.params)
        xredChiSq = np.sqrt(xfit.fnorm/dof)
        yredChiSq = np.sqrt(yfit.fnorm/dof)

        print xredChiSq
        print yredChiSq
        print xfit.fnorm
        print yfit.fnorm
        
        simTime = np.arange(np.amin(years)-t0,np.amax(years)+0.1-t0,0.1)
        simX = self.poly2(xfit.params,simTime)
        simY = self.poly2(yfit.params,simTime)
        #plot(simX,simY)

        modelAccelX = self.poly2(xfit.params, years - t0)
        modelAccelY = self.poly2(yfit.params, years - t0)

        xfitChi2 = np.sum((modelAccelX - x)**2/xerr**2.0)
        yfitChi2 = np.sum((modelAccelY - y)**2/yerr**2.0)
        print 'xfitchi2 %f' % xfitChi2
        print 'yfitchi2 %f' % yfitChi2
        # fit a line
        guessX =  [self.x[ind],(np.amax(x) - np.amin(x))/(np.amax(years)-np.amin(years))]
        guessY =  [self.y[ind],(np.amax(y) - np.amin(y))/(np.amax(years)-np.amin(years))]
        xfitLin = nmpfit_sy.mpfit(self.fitLin, guessX, functkw=functargsX,quiet=1)
        yfitLin = nmpfit_sy.mpfit(self.fitLin, guessY, functkw=functargsY,quiet=1)

        
        simXLin = self.line(xfitLin.params,simTime)
        simYLin = self.line(yfitLin.params,simTime)
        #plot(simXLin,simYLin)

        # if the accleration fit has too large a chi-square, go
        # through and remove points that are bad. then refit until the
        # chi-square goes below the threshold

        # initialize the chi square values
        loopChiX = xfit.fnorm/dof
        loopChiY = yfit.fnorm/dof

        loopLinChiX = xfitLin.fnorm/(dof + 1.0)
        loopLinChiY = yfitLin.fnorm/(dof + 1.0)
        xPoints = x
        xPointsErr = xerr
        yPoints = y
        yPointsErr = yerr
        subTime = years - t0


        print xfit.params
        print yfit.params
        x0 = self.poly2(xfit.params, 0.0)
        y0 = self.poly2(yfit.params, 0.0)
        print 'running isPhyiscal on polyfit results'
        print self.isPhysical(self.x[ind],self.y[ind],self.ax[ind],self.axe[ind],self.ay[ind],self.aye[ind])
        
        print 'running isPhysical on fit: '
        accelStatus = self.isPhysical(x0, y0,xfit.params[2],xfit.perror[2],
                                      yfit.params[2], yfit.perror[2], arcsec = True)
        print accelStatus
        print 't0: ', t0
        print 'from poly fit'
        print 'x0, vx, ax:'
        print [self.x[ind],self.vx[ind],self.ax[ind]]
        print [self.xerr[ind],self.vxe[ind],self.axe[ind]]
        print 'y0, yx, yx:'
        print [self.y[ind],self.vy[ind],self.ay[ind]]
        print [self.yerr[ind],self.vye[ind],self.aye[ind]]
##         print '[x, y, ax, ay]'
##         print [self.x[ind],self.y[ind],self.ax[ind],self.ay[ind]]
##         print '[xerr, yerr, axerr, ayerr]'
##         print [self.xerr[ind],self.yerr[ind],self.axe[ind],self.aye[ind]]
        print 'from current accel fit'
        print 'xfit: '
        print xfit.params
        print xfit.perror
        print 'yfit: '
        print yfit.params
        print yfit.perror
                
        print [x0, y0, xfit.params[2], yfit.params[2]]
        #yfit.perror = yfit.perror*np.sqrt(yfit.fnorm/dof)
        #xfit.perror = xfit.perror*np.sqrt(xfit.fnorm/dof)
        print [xfit.perror[0], yfit.perror[0], xfit.perror[2], yfit.perror[2]]
        print sqrt(xfit.perror[2]**2+yfit.perror[2]**2)
        print accelStatus

        print 'Start chi2x: %7.3f chi2y: %7.3f' % (loopChiX, loopChiY)
        print 'Start linear chi2x: %7.3f chi2y: %7.3f' % (loopLinChiX, loopLinChiY)

        # loop until a good chi2 is obtained or if the number of points goes below five
        while (((loopChiX > chiThreshold) | (loopChiY > chiThreshold)) & (len(xPoints) > 4)):
            
            # look for the difference between model and data
            xModel = self.poly2(xfit.params, subTime)
            yModel = self.poly2(yfit.params, subTime)
            xDiff = (xModel - xPoints)/xPointsErr
            yDiff = (yModel - yPoints)/yPointsErr
            totalDiff = sqrt(xDiff**2 + yDiff**2)

            xLinModel = self.line(xfitLin.params, subTime)
            yLinModel = self.line(yfitLin.params, subTime)
            xLinDiff = (xModel - xPoints)/xPointsErr
            yLinDiff = (yModel - yPoints)/yPointsErr
            totalLinDiff = sqrt(xLinDiff**2 + yLinDiff**2)

            # find the maximum difference
            maxInd = np.argmax(totalDiff)
            print 'Epoch %7.3f is off by %7.3f sigma' % (subTime[maxInd]+t0, totalDiff[maxInd])

            # look at how the fit imporves if we remove either end point
            xPointsEnd = xPoints[0:-1]
            xPointsErrEnd = xPointsErr[0:-1]
            yPointsEnd = yPoints[0:-1]
            yPointsErrEnd = yPointsErr[0:-1]
            subTimeEnd = subTime[0:-1]
            
            functargsXEnd = {'x':subTimeEnd, 'y':xPointsEnd, 'err': xPointsErrEnd}
            functargsYEnd = {'x':subTimeEnd, 'y':yPointsEnd, 'err': yPointsErrEnd}

            p0xEnd = [xPointsEnd[0],(np.amax(xPointsEnd) - np.amin(xPointsEnd))/(np.amax(subTimeEnd)-np.amin(subTimeEnd)), 0]
            p0yEnd = [yPointsEnd[0],(np.amax(yPointsEnd) - np.amin(yPointsEnd))/(np.amax(subTimeEnd)-np.amin(subTimeEnd)), 0]
        
            # do the fitting
            xfitEnd = nmpfit_sy.mpfit(self.fitfunPoly2, p0xEnd, functkw=functargsXEnd,quiet=1)
            yfitEnd = nmpfit_sy.mpfit(self.fitfunPoly2, p0yEnd, functkw=functargsYEnd,quiet=1)
            dofEnd = len(xPointsEnd) - len(xfitEnd.params)
            loopChiXEnd = xfitEnd.fnorm/dofEnd
            loopChiYEnd = yfitEnd.fnorm/dofEnd

            
            # remove max difference and refit
            xPoints = np.delete(xPoints, maxInd)
            xPointsErr = np.delete(xPointsErr, maxInd)
            yPoints = np.delete(yPoints, maxInd)
            yPointsErr = np.delete(yPointsErr, maxInd)
            subTime = np.delete(subTime,maxInd)
            
            functargsX = {'x':subTime, 'y':xPoints, 'err': xPointsErr}
            functargsY = {'x':subTime, 'y':yPoints, 'err': yPointsErr}

            p0x = [xPoints[0],(np.amax(xPoints) - np.amin(xPoints))/(np.amax(subTime)-np.amin(subTime)), 0]
            p0y = [yPoints[0],(np.amax(yPoints) - np.amin(yPoints))/(np.amax(years)-np.amin(years)), 0]
        
            # do the fitting
            xfit = nmpfit_sy.mpfit(self.fitfunPoly2, p0x, functkw=functargsX,quiet=1)
            yfit = nmpfit_sy.mpfit(self.fitfunPoly2, p0y, functkw=functargsY,quiet=1)
            dof = len(xPoints) - len(xfit.params)
            loopChiX = xfit.fnorm/dof
            loopChiY = yfit.fnorm/dof

            # fit a line
            guessX =  [self.x[ind],(np.amax(x) - np.amin(x))/(np.amax(years)-np.amin(years))]
            guessY =  [self.y[ind],(np.amax(y) - np.amin(y))/(np.amax(years)-np.amin(years))]
            xfitLin = nmpfit_sy.mpfit(self.fitLin, guessX, functkw=functargsX,quiet=1)
            yfitLin = nmpfit_sy.mpfit(self.fitLin, guessY, functkw=functargsY,quiet=1)

            loopLinChiX = xfitLin.fnorm/(dof + 1.0)
            loopLinChiY = yfitLin.fnorm/(dof + 1.0)

            print 'New chi2x: %7.3f chi2y: %7.3f' % (loopChiX, loopChiY)
            print 'New linear chi2x: %7.3f chi2y: %7.3f' % (loopLinChiX, loopLinChiY)
            print 'New chi2x without end point: %7.3f chi2y: %7.3f' % (loopChiXEnd, loopChiYEnd)
                        


        x0 = self.poly2(xfit.params, 0.0)
        y0 = self.poly2(yfit.params, 0.0)
        accelStatus = self.isPhysical(x0, y0,xfit.params[2],xfit.perror[2],
                                      yfit.params[2], yfit.perror[2], arcsec = True)
        print accelStatus
        
        if mkplots:
            mTime = arange(np.min(subTime), np.max(subTime), 0.05)
            xModel = self.poly2(xfit.params, mTime)
            yModel = self.poly2(yfit.params, mTime)

            xLinModel = self.line(xfitLin.params, mTime)
            yLinModel = self.line(yfitLin.params, mTime)
            clf()
            py.errorbar(x, y, xerr, yerr, fmt = 'o',color='red')
            py.errorbar(xPoints, yPoints, xPointsErr, yPointsErr, fmt = 'o',color='blue')
            plot(xModel, yModel,color='red')
            plot(xLinModel, yLinModel, color = 'green')
            xlim(np.max(xPoints+xPointsErr), np.min(xPoints-xPointsErr))
            title(starName)
            
            # plot the original fit
            
            for ii in range(len(x)):
                text(x[ii],y[ii], years[ii])

            # find the median point and find the velocity between each point    
            d1 = np.sqrt(x**2 + y**2)
            d1Err = np.sqrt((x**2*xerr**2 + y**2*yerr**2)/d1**2)
            medianInd = realMedian(d1,index=True)
            plot(x[medianInd],y[medianInd],'co')

            # find distance from the median point
            d1 = np.sqrt((x-x[medianInd])**2+(y-y[medianInd])**2)
            
            d2 = np.delete(d1, medianInd)
            d2Err = np.delete(d1, medianInd)
            d2Year = np.delete(years, medianInd)

            d1 = d1[medianInd]
            d1Err = d1Err[medianInd]
            d1Year = years[medianInd]
            for i in range(len(d2)):
                t = array([d1Year, d2Year[i]]) - t0
                d3 = array([d1, d2[i]])
                d3err = array([d1Err, d2Err[i]])

                # sort the time so that all points will fit in one direction
                vFit = self.fitVelocity(t, d3, d3err)
                #print d2Year[i]
                #print np.abs(vFit.params[1])

        if bootStrap:    
            self.bootstrapVel(x, y, xerr, yerr, years, t0, iterate=iterate)
            
    def fitVelocity(self, x, y,  yerr):
        """ takes in an array of x and y position as well as time and
        fit for the velocity. Will return an mpfit object
        """
        guessY =  [y[0],(np.amax(y) - np.amin(y))/(np.amax(x)-np.amin(x))]
        functargsY = {'x':x, 'y':y, 'err': yerr}
        yfitLin = nmpfit_sy.mpfit(self.fitLin, guessY, functkw=functargsY,quiet=1)

        return yfitLin

    def fitAccel(self, x, y,  yerr):
        """ takes in an array of x and y position as well as time and
        fit for the acceleration. Will return an mpfit object
        """
        guessY =  [y[0],(np.amax(y) - np.amin(y))/(np.amax(x)-np.amin(x)), 0]
        functargsY = {'x':x, 'y':y, 'err': yerr}
        yfitLin = nmpfit_sy.mpfit(self.fitfunPoly2, guessY, functkw=functargsY,quiet=1)

        return yfitLin

    def fitJerk(self, x, y,  yerr):
        """ takes in arrays for x and y position and time and
        fit for jerk. Will return an mpfit object
        """
        guessY =  [y[0],(np.amax(y) - np.amin(y))/(np.amax(x)-np.amin(x)), 0, 0]
        functargsY = {'x':x, 'y':y, 'err': yerr}
        yfitLin = nmpfit_sy.mpfit(self.fitfunPoly3, guessY, functkw=functargsY,quiet=1)

        return yfitLin

    def bootstrapVel(self, xInput, yInput, xerrInput, yerrInput, yearsInput,
                     t0, n = 500, fraction = 0.05, iterate = False):
        """ Do a half-sample bootstrap of the velocity fits to
        determine if there are incorrect points.

        KEYWORDS: n = 500 - number of bootstraps to do
                  fraction = 0.05 - fraction of times a point is detected
                                    to be kept
        """
        xStack = np.zeros(n)
        yStack = np.zeros(n)
        velx = np.zeros(n)
        velxErr = np.zeros(n)
        velxChi2 = np.zeros(n)
        vely = np.zeros(n)
        velyErr = np.zeros(n)
        velyChi2 = np.zeros(n)
        drawStack = np.zeros((len(xInput)/2, n))

        # degree of freedom for each half sample boot strap
        bootStrapDOF = len(xInput)/2.0 - 2.0
        chiThreshold = bootStrapDOF*3.0
        print 'boot strap chi2 threshold: ', chiThreshold
        xMedianErr = realMedian(xerrInput)
        yMedianErr = realMedian(yerrInput)

        for i in range(n):
            # draw N/2 indices without replacement
            draw = random.sample(arange(len(xInput)), len(xInput)/2)
            drawStack[:,i] = draw
            x = xInput[draw]
            y = yInput[draw]

            # make all the errors the median error (won't really work
            # because this will make points that are badly measured
            # have low errors)
            
            # xerr = np.zeros(len(xInput)/2)+xMedianErr
            # yerr = np.zeros(len(xInput)/2)+yMedianErr
            
            xerr = xerrInput[draw]
            yerr = yerrInput[draw]

            years = yearsInput[draw]-t0
            #print '%f xerr' % i
            #print x
            #print xerr
            #print years
            xfit = self.fitVelocity(years, x, xerr)
            yfit = self.fitVelocity(years, y, yerr)
            
            dof = len(x) - len(xfit.params)
            xredChiSq = xfit.fnorm
            yredChiSq = yfit.fnorm

            #print xfit.params
            xStack[i] = xfit.params[0]
            yStack[i] = yfit.params[0]
            velx[i] = xfit.params[1]
            #velxErr[i] = xfit.perror[1]
            
            velxChi2[i] = xredChiSq
            vely[i] = yfit.params[1]
            #velyErr[i] = yfit.perror[1]
            velyChi2[i] = yredChiSq

        # compute the tangential and radial directions to the velocity
        xMed = realMedian(xStack)
        yMed = realMedian(yStack)
        r = np.sqrt(xMed**2 + yMed**2)
        vr = ((velx*xMed) + (vely*yMed))/r
        vt = ((velx*yMed) - (vely*xMed))/r
        clf()
        
        subplot(231)
        mVr = np.mean(vr)
        sVr = np.std(vr)
        mVt = np.mean(vt)
        sVt = np.std(vt)

        binWidth = (np.max(vr) - np.min(vr))/20.0  # make sure there are at least 10 bins
        
        nbins, bins, patches1 = hist(vr, bins = arange(np.min(vr), np.max(vr), binWidth))
        xlabel('Vr')
        ylabel('N stars')
        subplot(232)
        binWidth = (np.max(vt) - np.min(vt))/20.0
        #nbins, bins, patches1 = hist(vt, bins = arange(mVt - sVt*10.0, mVt+sVt*10.0, 0.0001))
        nbins, bins, patches1 = hist(vt, bins = arange(np.min(vt), np.max(vt), binWidth))
        xlabel('Vt')
        ylabel('N stars')
        subplot(233)
        print 'min chi2 x: %6.3f, max chi2 x: %6.3f' % (np.min(velxChi2), np.max(velxChi2))
        nbins, bins, patches1 = hist(velxChi2, bins = arange(0, self.chiThreshold, self.chiThreshold/30.0))
        xlabel('X Vel. Chi-Sq')
        ylabel('N Stars')
        subplot(234)
        nbins, bins, patches1 = hist(velyChi2, bins = arange(0, self.chiThreshold, self.chiThreshold/30.0))
        print 'min chi2 y: %6.3f, max chi2 y: %6.3f' % (np.min(velyChi2), np.max(velyChi2))
        xlabel('Y Vel. Chi-Sq')
        ylabel('N Stars')
        subplot(235)
        good = np.where((velxChi2 < chiThreshold) & (velyChi2 < chiThreshold))[0]
        subStack = drawStack[:,good]
        nbins, bins, patches = hist(subStack.flatten(),bins = np.arange(len(xInput)+1))
        xlabel('Epoch')
        ylabel('Times Used')
        print bins
        py.xticks(arange(len(xInput)), yearsInput, rotation = 45)
        #nbins, bins, patches1 = hist(velyChi2, bins = arange(0, 100, 0.5))
        nbins = np.array(nbins,dtype='float')
        nbins = nbins/np.sum(nbins)
        print nbins

        # select out all the bad points that aren't used often. Take
        # out the points that are found in fewer than a fraction of
        # the maximum bin. This is done so that if stars are found in
        # all epochs, there could potentially be a case where all the
        # points will be dropped.
        bad = np.where(nbins < fraction*np.max(nbins))[0]

        # recompute the fits without the bad points
        print 'dropping ', yearsInput[bad]
        x = np.delete(xInput, bad)
        xerr = np.delete(xerrInput, bad)
        y = np.delete(yInput, bad)
        yerr = np.delete(yerrInput, bad)                
        years = np.delete(yearsInput, bad)-t0
        print years
        print x
        xfit = self.fitVelocity(years, x, xerr)
        yfit = self.fitVelocity(years, y, yerr)
        subplot(236)
        py.errorbar(xInput, yInput, xerrInput, yerrInput, fmt='ro')
        py.errorbar(x,y,xerr,yerr,fmt='o')
        plot(self.line(xfit.params,years),self.line(yfit.params,years))
        xlabel('RA offset (arcsec)')
        ylabel('DEC offset (arcsec)')
        velChi2x = xfit.fnorm
        velChi2y = yfit.fnorm
        print 'xfit', xfit.params
        print 'xfit error: ', xfit.perror
        print 'xfit chi2: ', xfit.fnorm
        print 'yfit', yfit.params
        print 'yfit error: ', yfit.perror
        print 'yfit chi2: ', yfit.fnorm
        xfit = self.fitAccel(years, x, xerr)
        yfit = self.fitAccel(years, y, yerr)
        print 'acceleration fits:'
        print 'xfit', xfit.params
        print 'xfit error: ', xfit.perror
        print 'xfit chi2: ', xfit.fnorm
        print 'yfit', yfit.params
        print 'yfit error: ', yfit.perror
        print 'yfit chi2: ', yfit.fnorm
        accelChi2x = xfit.fnorm
        accelChi2y = yfit.fnorm

        # the F-test for which model to prefer
        fValue = ((velChi2x - accelChi2x)/1.0)/(accelChi2x/(len(x)-3.0))
        print 'F-test value X fit', fValue
        print 'F-test probability: ', stats.f.sf(fValue, 1, len(x) -3.0)
        fValue = ((velChi2y - accelChi2y)/1.0)/(accelChi2x/(len(x)-3.0))
        print 'F-test value Y fit', fValue
        print 'F-test probability: ', stats.f.sf(fValue, 1, len(x) -3.0)

        if iterate:
            self.bootstrapVel(x, y, xerr, yerr, years+t0,
                         t0, n = 500, fraction = 0.05, iterate = False)

        
        
    def computeAccel(self):
        """Recomputes the velocity and acceleration of all sources and
        updates those fields as well as the chi2 values. This will use
        the internal fitter instead of the one from polyfit
        """
        # get variables into local scope
        cc = self.cc
        names = list(self.names)
        x, y = self.x, self.y
        xerr, yerr = self.xerr, self.yerr
        vx, vy = self.vx, self.vy

        # loop through the names and get the points for the acceleration fits
        for i in range(len(names)):
            ind = names.index(names[i])

            x = self.starSet.stars[ind].getArrayAllEpochs('x')
            y = self.starSet.stars[ind].getArrayAllEpochs('y')
            xerr_p = self.starSet.stars[ind].getArrayAllEpochs('xerr_p')
            yerr_p = self.starSet.stars[ind].getArrayAllEpochs('yerr_p')

            xerr_a = self.starSet.stars[ind].getArrayAllEpochs('xerr_a')
            yerr_a = self.starSet.stars[ind].getArrayAllEpochs('yerr_a')

            xerr = np.sqrt(xerr_p**2 + xerr_a**2)
            yerr = np.sqrt(yerr_p**2 + yerr_a**2)

            t0 = self.starSet.stars[ind].fitXa.t0    # use the same T0 as polyfit
            years = np.array(self.starSet.stars[ind].years)

            good = (np.where(abs(x) < 500))[0]

            if (len(good) > 3):
                x = x[good]
                y = y[good]
                xerr = xerr[good]
                yerr = yerr[good]
                years = years[good]

                #functargsX = {'x':years-t0, 'y':x, 'err': xerr}
                #functargsY = {'x':years-t0, 'y':y, 'err': yerr}

                #p0x = [self.x[ind],(np.amax(x) - np.amin(x))/(np.amax(years)-np.amin(years)), 0]
                #p0y = [self.y[ind],np.amax(y) - np.amin(y)/(np.amax(years)-np.amin(years)), 0]
                #xfit = nmpfit_sy.mpfit(self.fitfunPoly2, p0x, functkw=functargsX,quiet=1)
                #yfit = nmpfit_sy.mpfit(self.fitfunPoly2, p0y, functkw=functargsY,quiet=1)
                xfit = self.fitAccel(years - t0, x, xerr)
                yfit = self.fitAccel(years - t0, y, yerr)
                dof = len(x) - len(xfit.params)
                if names[i] in ['irs29N']:
                    print 'yerr_a: ', yerr_a
                    print 'yerr_p: ', yerr_p
                    print 'yerr: ', yerr
                # put the results in the correct star
                if names[i] in self.nonPhysical:
                    print '%s    \t xAccel: %10.8f xPolyFit: %10.8f Diff: %10.8f  xErr: %10.8f polyErr: %10.8f ErrDiff: %10.8f' % (names[i], xfit.params[2], self.ax[i],
                                                                                               xfit.params[2] - self.ax[i], xfit.perror[2], self.axe[i],
                                                                                                                                self.vxe[i]/xfit.perror[1])
                self.x[i] = xfit.params[0]
                self.ax[i] = xfit.params[2]
                self.axe[i] = xfit.perror[2]
                self.vx[i] = xfit.params[1]
                self.vxe[i] = xfit.perror[1]

                self.y[i] = yfit.params[0]
                self.ay[i] = yfit.params[2]
                self.aye[i] = yfit.perror[2]
                self.vy[i] = yfit.params[1]
                self.vye[i] = yfit.perror[1]
        self.updateAccel() # update the radial and tangential acceleration calculations


    def computeFTest(self):
        """ Go through the acceleration and velocity chi-squares to
        compute the F test values.
        """
        print 'computing F-Test for all valid stars '
        chi2xAccel = self.chi2x
        chi2yAccel = self.chi2y
        chi2xVel = self.chi2xv
        chi2yVel = self.chi2yv
        nEpochs = self.nEpochs
        
        xFValue = zeros(len(chi2xAccel))
        yFValue = zeros(len(chi2xAccel))
        xFProb = zeros(len(chi2xAccel))
        yFProb = zeros(len(chi2xAccel))
        dof = zeros(len(chi2xAccel))
        
        good = np.where(chi2xAccel > 0)[0]
        dof[good] = nEpochs[good] - 3.0  # dof freedom for acceleration
        
        xFValue[good] = ((chi2xVel[good] - chi2xAccel[good])/1.0)/(chi2xAccel[good]/dof[good])
        yFValue[good] = ((chi2yVel[good] - chi2yAccel[good])/1.0)/(chi2yAccel[good]/dof[good])

        xFProb[good] = stats.f.sf(xFValue[good], 1, dof[good])
        yFProb[good] = stats.f.sf(yFValue[good], 1, dof[good])

        self.xFProb = xFProb
        self.yFProb = yFProb
        

    def computeFTestJerk(self):
        """ Go through the jerk and acceleration chi-squares to
        compute the F test values.
        """
        print 'computing Jerk F-Test for all valid stars '
#        self.computeJerk()
        chi2xAccel = self.chi2x
        chi2yAccel = self.chi2y
        chi2xJerk = self.chi2xj
        chi2yJerk = self.chi2yj
        nEpochs = self.nEpochs
        
        xFValue = zeros(len(chi2xAccel))
        yFValue = zeros(len(chi2xAccel))
        xFProb = zeros(len(chi2xAccel))
        yFProb = zeros(len(chi2xAccel))
        dof = zeros(len(chi2xAccel))
        
        good = np.where(chi2xJerk > 0)[0]
        dof[good] = nEpochs[good] - 4.0  # dof freedom for jerk 
        
        xFValue[good] = ((chi2xAccel[good] - chi2xJerk[good])/1.0)/(chi2xJerk[good]/dof[good])
        yFValue[good] = ((chi2yAccel[good] - chi2yJerk[good])/1.0)/(chi2yJerk[good]/dof[good])

        xFProb[good] = stats.f.sf(xFValue[good], 1, dof[good])
        yFProb[good] = stats.f.sf(yFValue[good], 1, dof[good])

        self.xFProbJrk = xFProb
        self.yFProbJrk = yFProb
        

    def computeSpeckleAOVel(self, requireAllEpochs = False, verbose=True):
        """ compute the speckle and AO velocities individually and
        check the difference in the velocity vectors.

        Keywords: requireAllEpochs - require the stars to be detected
                  in all speckle or all AO epochs to be computed
        """
        names = self.names
        starsX = self.x
        starsY = self.y
        
        # look up the epoch for speckle and AO
        epochFile = self.rootDir+'scripts/epochsInfo.txt'

        if os.path.isfile(epochFile):
            epochTab = asciidata.open(epochFile)
            aoFlag = np.array(epochTab[2])
            alignFlag = np.array(epochTab[3])

            bad = where(alignFlag == 0)[0]
            aoFlag = np.delete(aoFlag,bad)
            
            # separate the speckle points from the AO points
            speckInd = np.where((aoFlag == 0))[0]
            aoInd = np.where((aoFlag == 1))[0]

            if verbose == True:
                print 'speckle epochs: ', self.allEpochs[speckInd]
                print 'AO epochs: ', self.allEpochs[aoInd]
        else:
            return False

        # define some new arrays for velocity calculations
        n = len(names)
        speckVx = zeros(n)
        speckVy = zeros(n)
        speckVxErr = zeros(n)
        speckVyErr = zeros(n)
        speckT0 = zeros(n)
        nSpeckle = zeros(n)  # number of speckle epochs
        speckChi2x = zeros(n) # chi2 for speckle velocities
        speckChi2y = zeros(n) # chi2 for speckle velocities
        
        aoVx = zeros(n)
        aoVy = zeros(n)
        aoVxErr = zeros(n)
        aoVyErr = zeros(n)
        aoT0 = zeros(n)
        aoChi2x = zeros(n) # chi2 for AO velocities
        aoChi2y = zeros(n) # chi2 for AO velocities
        nAO = zeros(n) # number of AO epochs
        
        velDiffX = zeros(n)
        velDiffY = zeros(n)
        
        # loop through names
        #for ind in arange(n):
        for ind in arange(n):
            
            x = -self.starSet.stars[ind].getArrayAllEpochs('x')
            y = self.starSet.stars[ind].getArrayAllEpochs('y')
            # positional errors
            xerr_p = self.starSet.stars[ind].getArrayAllEpochs('xerr_p')
            yerr_p = self.starSet.stars[ind].getArrayAllEpochs('yerr_p')

            # alignment errors
            xerr_a = self.starSet.stars[ind].getArrayAllEpochs('xerr_a')
            yerr_a = self.starSet.stars[ind].getArrayAllEpochs('yerr_a')

            # add the two errors in quadrature
            xerr = np.sqrt(xerr_p**2 + xerr_a**2)
            yerr = np.sqrt(yerr_p**2 + yerr_a**2)

            t0 = self.starSet.stars[ind].fitXa.t0    # use the same T0 as polyfit
            years = np.array(self.starSet.stars[ind].years)
            
            speckX = x[speckInd]
            speckY = y[speckInd]
            speckXErr = xerr[speckInd]
            speckYErr = yerr[speckInd]
            speckYears = years[speckInd]

            aoX = x[aoInd]
            aoY = y[aoInd]
            aoXErr = xerr[aoInd]
            aoYErr = yerr[aoInd]
            aoYears = years[aoInd]
            
            goodS = (np.where(abs(speckX) < 500))[0]
            nSpeckle[ind] = len(goodS)
            if (len(goodS) > 3):
                speckX = speckX[goodS]
                speckY = speckY[goodS]
                speckXErr = speckXErr[goodS]
                speckYErr = speckYErr[goodS]
                speckYears = speckYears[goodS]
                speckT0[ind] = np.sum(speckYears/speckXErr**2)/np.sum(1.0/speckXErr**2)
                speckYears = speckYears - speckT0[ind]
                
            goodAO = (np.where(abs(aoX) < 500))[0]
            nAO[ind] = len(goodAO)
            if (len(goodAO) > 3):
                aoX = aoX[goodAO]
                aoY = aoY[goodAO]
                aoXErr = aoXErr[goodAO]
                aoYErr = aoYErr[goodAO]
                aoYears = aoYears[goodAO]
                aoT0[ind] = np.sum(aoYears/aoXErr**2)/np.sum(1.0/aoXErr**2)
                aoYears = aoYears - aoT0[ind]

            if requireAllEpochs:
                # only do fitting if both speckle and AO are found in all epochs
                # print 'requiring stars to be detected in all speckle and AO epochs'
                filterTest = (len(goodS) == len(speckInd)) & (len(goodAO) == len(aoInd))
            else:
                # do fitting if stars are found in more than three epochs
                filterTest = (len(goodAO) > 3) & (len(goodS) > 3)

            if filterTest:
                functargsX = {'x':speckYears, 'y':speckX, 'err': speckXErr}
                functargsY = {'x':speckYears, 'y':speckY, 'err': speckYErr}
                guessX =  [mean(speckX),(np.amax(speckX) - np.amin(speckX))/(np.amax(speckYears)-np.amin(speckYears))]
                guessY =  [mean(speckY),(np.amax(speckY) - np.amin(speckY))/(np.amax(speckYears)-np.amin(speckYears))]
                xfitLin = nmpfit_sy.mpfit(self.fitLin, guessX, functkw=functargsX,quiet=1)
                yfitLin = nmpfit_sy.mpfit(self.fitLin, guessY, functkw=functargsY,quiet=1)
                print 'speckle fit x: ', xfitLin.params, xfitLin.perror
                print 'speckle fit y : ', yfitLin.params, yfitLin.perror
                speckVx[ind] = xfitLin.params[1]
                if xfitLin.perror == None:
                    speckVxErr[ind] = 0
                else:
                    speckVxErr[ind] = xfitLin.perror[1]
                    speckChi2x[ind] = xfitLin.fnorm
                    
                speckVy[ind] = yfitLin.params[1]
                if yfitLin.perror == None:
                    speckVyErr[ind] = 0
                else:
                    speckVyErr[ind] = yfitLin.perror[1]
                    speckChi2y[ind] = yfitLin.fnorm
                
                

                # fit the AO portion
                functargsX = {'x':aoYears, 'y':aoX, 'err': aoXErr}
                functargsY = {'x':aoYears, 'y':aoY, 'err': aoYErr}
                guessX =  [mean(aoX),(np.amax(aoX) - np.amin(aoX))/(np.amax(aoYears)-np.amin(aoYears))]
                guessY =  [mean(aoY),(np.amax(aoY) - np.amin(aoY))/(np.amax(aoYears)-np.amin(aoYears))]
                xfitLin = nmpfit_sy.mpfit(self.fitLin, guessX, functkw=functargsX,quiet=1)
                yfitLin = nmpfit_sy.mpfit(self.fitLin, guessY, functkw=functargsY,quiet=1)
                print 'ao fit: ', xfitLin.params, xfitLin.perror
                aoVx[ind] = xfitLin.params[1]
                if xfitLin.perror == None:
                    aoVxErr[ind] = 0
                else:
                    aoVxErr[ind] = xfitLin.perror[1]
                    aoChi2x[ind] = xfitLin.fnorm
                    
                aoVy[ind] = yfitLin.params[1]
                if yfitLin.perror == None:
                    aoVyErr[ind] = 0
                else:
                    aoVyErr[ind] = yfitLin.perror[1]
                    aoChi2y[ind] = yfitLin.fnorm


                # take the difference between the two velocities
                velDiffX[ind] = aoVx[ind] - speckVx[ind]
                velDiffY[ind] = aoVy[ind] - speckVy[ind]

        self.speckVx = speckVx
        self.speckVy = speckVy
        self.speckVxErr = speckVxErr
        self.speckVyErr = speckVyErr
        self.nSpeckle = nSpeckle
        self.speckChi2x = speckChi2x
        self.speckChi2y = speckChi2y
        
        self.aoVx = aoVx
        self.aoVy = aoVy
        self.aoVxErr = aoVxErr
        self.aoVyErr = aoVyErr
        self.nAO = nAO
        self.aoChi2x = aoChi2x
        self.aoChi2y = aoChi2y
        self.velDiffX = velDiffX
        self.velDiffY = velDiffY
        self.plotSpeckleAOVel()

    def chi2DistFit(self, p, direction = 'x', returnChi2 = False, maxChi2 = None,
                    scaleFactor = 1.0, data = 'both', removeConfused = True):
        """ function to fit the chi2 distribution of velocities to
        figure out the additive error that needs to be included to
        bring the chi2 distribution in line with the expected. Will only
        consider stars that are detected in all epochs. 

        Input: p - additive error parameter
               x - bin locations
               y - theoretical chi2 curve

        Keywords: direction - 'x' or 'y' to calculate the chi2 in
                  either direction.

                  returnChi2 - return the chi2
                  values instead of the difference between chi2 and
                  model

                  maxChi2 - the maximum chi2 value to create the
                  histogram. Defaults to 3 times the number of epochs

                  removeConfused - look at the self.confusedStars array and not
                  include them in the fit

                  self.rmax - the maximum radius to compute the additive
                  error. Default rmax = 0 for the entire range.
        """
        names = self.names
        maxEpochInd = self.maxEpochInd
        mag = self.mag
        
        # remove the points that are confused
        if removeConfused:
            # remove points from maxEpochInd that are in the confusedSourcesInd
            maxEpochInd2 = np.array(np.setdiff1d(maxEpochInd, self.confusedSourcesInd),dtype='int32')
            print 'removing %f confused sources out of %f total sources' % (len(maxEpochInd)-len(maxEpochInd2),len(maxEpochInd))
            maxEpochInd = maxEpochInd2

        # limit the additive error within a certain radial range
        if (self.rmax > 0):
            goodR = np.where((self.r2d < self.rmax) & (self.r2d > self.rmin))[0]
            maxEpochInd = np.intersect1d(goodR, maxEpochInd)
            print 'using only stars within %f arcsec, faintest magnitude: %f' % (self.rmax,np.amax(mag[maxEpochInd]))

        if (self.magmax > 0):
            goodMag = np.where((mag <= self.magmax) & (mag > self.magmin))[0]
            maxEpochInd = np.intersect1d(goodMag, maxEpochInd)
            print 'limiting magnitude to those brighter than %f, n stars: %f' % (self.magmax,len(maxEpochInd))
            
        n = len(maxEpochInd)

        # keep track of which stars were used to compute the additive error
        self.addErrStarsInd = maxEpochInd
        
        chi2x = zeros(n)
        chi2y = zeros(n)

        # allow for the possibility of subtracting some error
        if p[0] < 0:
            sign = -1.0
        else:
            sign = 1.0

        # make sure if there are no speckle data that we set the flag to AO only
        if len(self.speckInd) == 0:
            data = 'ao'
        
        if (data == 'both'):
            epochInd = np.union1d(self.speckInd, self.aoInd)
        elif (data == 'speckle'):
            epochInd = self.speckInd
        elif (data == 'ao'):
            epochInd = self.aoInd

        print 'using error scale factor: ', scaleFactor

        # remove confused sources

        
        # go through good stars and compute the chi2 values
        for i in arange(n):
            ind = maxEpochInd[i]

            if (direction == 'x'):
                xpts = -self.starSet.stars[ind].getArrayAllEpochs('x')[epochInd]
                # positional errors - divide by sqrt(2) to get error in mean
                xerr_p = self.starSet.stars[ind].getArrayAllEpochs('xerr_p')[epochInd]/scaleFactor
                xerr_a = self.starSet.stars[ind].getArrayAllEpochs('xerr_a')[epochInd]
            elif (direction == 'y'):
                xpts = self.starSet.stars[ind].getArrayAllEpochs('y')[epochInd]
                xerr_p = self.starSet.stars[ind].getArrayAllEpochs('yerr_p')[epochInd]/scaleFactor
                xerr_a = self.starSet.stars[ind].getArrayAllEpochs('yerr_a')[epochInd]
            elif (direction == 'both'):
                xpts = -self.starSet.stars[ind].getArrayAllEpochs('x')[epochInd]
                # positional errors - divide by sqrt(2) to get error in mean
                xerr_p = self.starSet.stars[ind].getArrayAllEpochs('xerr_p')[epochInd]/scaleFactor
                xerr_a = self.starSet.stars[ind].getArrayAllEpochs('xerr_a')[epochInd]
                
                ypts = self.starSet.stars[ind].getArrayAllEpochs('y')[epochInd]
                yerr_p = self.starSet.stars[ind].getArrayAllEpochs('yerr_p')[epochInd]/scaleFactor
                yerr_a = self.starSet.stars[ind].getArrayAllEpochs('yerr_a')[epochInd]


            # add the two errors in quadrature and add the new error
            xerr = np.sqrt(xerr_p**2 + xerr_a**2 + sign*p[0]**2)
            if (direction == 'both'):
                yerr = np.sqrt(yerr_p**2 + yerr_a**2 + sign*p[0]**2)
            

            years = np.array(self.starSet.stars[ind].years)[epochInd]
            if (data == 'both'):
                t0 = self.starSet.stars[ind].fitXa.t0    # use the same T0 as polyfit
            else:
                # compute TO with weighted mean epoch by x positional error
                t0 = np.sum(years/xerr**2)/np.sum(1.0/xerr**2)

            # go through each star and compute the velocity as well as the errors
            functargsX = {'x':years - t0, 'y':xpts, 'err': xerr}
            guessX = [mean(xpts),(np.amax(xpts) - np.amin(xpts))/(np.amax(years) - np.amin(years))]
            xfitLin = nmpfit_sy.mpfit(self.fitLin, guessX, functkw=functargsX,quiet=1)
            chi2x[i] = xfitLin.fnorm
            
            if (direction == 'both'):
                functargsY = {'x':years - t0, 'y':ypts, 'err': yerr}
                guessY = [mean(xpts),(np.amax(xpts) - np.amin(xpts))/(np.amax(years) - np.amin(years))]
                yfitLin = nmpfit_sy.mpfit(self.fitLin, guessY, functkw=functargsY,quiet=1)
                chi2y[i] = yfitLin.fnorm


        # put both together if doing both directions
        if (direction == 'both'):
            chi2x = np.concatenate((chi2x,chi2y))

        if returnChi2:
            return chi2x
        else:
            if maxChi2 is None:
                maxChi2 = len(years)*3
            width = maxChi2/30.0
            nStars, bins, patches = hist(chi2x, bins = arange(0,maxChi2, width),color='b')

            chiInd = bins
            # scale the theoretical chi-sq curve by the total number
            # of stars so we can use the poisson error for each bin
            dof = len(years) - 2.0
            print dof
            chi2Theory = stats.chi2.pdf(chiInd,dof)*np.sum(nStars)
            plot(chiInd, chi2Theory)
            plot(bins[0:-1],nStars,'go')
            good = np.where(nStars != 0)[0]
            bad = np.where(nStars == 0)[0]
            err = sqrt(nStars)
            err[bad] = 1.0
            diff = zeros(len(nStars))
            diff = (nStars - chi2Theory[0:-1])/err
            return diff

    def testChi2Fit(self, testErr = None, scaleFactor = 1.0, data = 'both'):
        """ test the chi2DistFit function by going through a series
        ofadditive factors to minimize the chi2
        """
        #testErr = np.arange(0.00001,0.0004,0.00001)
        if testErr is None:
            testErr = np.arange(0.00001,0.001,0.00001)
        chi2x = zeros(len(testErr))
        chi2y = zeros(len(testErr))

        # record the number of degree of freedom in the test
        if (data == 'speckle'):
            self.errAddDof = len(self.speckInd)-2.0
        elif data == ('ao'):
            self.errAddDof = len(self.aoInd)-2.0
        else:
            self.errAddDof = len(self.allEpochs)-2.0

        for ii in arange(len(testErr)):
            fit = self.chi2DistFit(np.array([testErr[ii]]), scaleFactor = scaleFactor,
                                   data = data)
            chi2x[ii] = np.sum(fit**2)
            print 'additive error factor: %f chi2: %f' % (testErr[ii],chi2x[ii])
            
        for ii in arange(len(testErr)):
            fit = self.chi2DistFit(np.array([testErr[ii]]), direction = 'y',
                                   scaleFactor = scaleFactor, data = data)
            chi2y[ii] = np.sum(fit**2)
            print 'additive error factor: %f chi2: %f' % (testErr[ii],chi2y[ii])

        clf()
        plot(testErr, chi2x, 'bo', label ='x')
        plot(testErr, chi2y, 'go', label = 'y')
        xlabel('Additive Factor (arcsec)')
        ylabel('Chi Sq. Difference')
        legend()
        self.errAdd = testErr
        self.errAddChi2x = chi2x
        self.errAddChi2y = chi2y
        savefig(self.plotPrefix+'additive_error_test.png')

    def testChi2FitXY(self, testErr = None, scaleFactor = None, data = 'both'):
        """ test the chi2DistFit function by going through a series
        ofadditive factors to minimize the chi2 over BOTH directions
        """

        if scaleFactor is None:
            scaleFactor = self.errScaleFactor
            
        if testErr is None:
            testErr = np.arange(0.00001,0.001,0.00001)

        chi2x = zeros(len(testErr))
        chi2y = zeros(len(testErr))

        # record the number of degree of freedom in the test
        if (data == 'speckle'):
            self.errAddDof = len(self.speckInd)-2.0
        elif data == ('ao'):
            self.errAddDof = len(self.aoInd)-2.0
        else:
            self.errAddDof = len(self.allEpochs)-2.0

        for ii in arange(len(testErr)):
            fit = self.chi2DistFit(np.array([testErr[ii]]),direction = 'both',
                                   scaleFactor = scaleFactor,data = data)
            chi2x[ii] = np.sum(fit**2)
            print 'additive error factor: %f chi2: %f' % (testErr[ii],chi2x[ii])
            
        clf()
        plot(testErr, chi2x, 'bo', label ='x')
        xlabel('Additive Factor (arcsec)')
        ylabel('Chi Sq. Difference')
        legend()
        self.errAdd = testErr
        self.errAddChi2xy = chi2x
        savefig(self.plotPrefix+'additive_error_test_xy.png')


    def computeJerk(self):
        """Computes the jerk for all sources. Uses the internal fitter for jerk.
              ***NOT FUNCTIONAL YET***
        """
        # get variables into local scope
        cc = self.cc
        names = list(self.names)
        x, y = self.x, self.y
        xerr, yerr = self.xerr, self.yerr
        vx, vy = self.vx, self.vy
        axe, aye = self.axe, self.aye

        #jx = np.zeros(len(names), dtype=float)
        #jxe = np.zeros(len(names), dtype=float)
        #jy = np.zeros(len(names), dtype=float)
        #jye = np.zeros(len(names), dtype=float)
        jx = []
        jxe = []
        jy = []
        jye = []
        tmp_chi2xj=np.zeros(len(names))
        tmp_chi2yj=np.zeros(len(names))

        # loop through the names and get the points for the jerk fits
        for i in range(len(names)):
            ind = names.index(names[i])

            x = self.starSet.stars[ind].getArrayAllEpochs('x')
            y = self.starSet.stars[ind].getArrayAllEpochs('y')
            xerr_p = self.starSet.stars[ind].getArrayAllEpochs('xerr_p')
            yerr_p = self.starSet.stars[ind].getArrayAllEpochs('yerr_p')

            xerr_a = self.starSet.stars[ind].getArrayAllEpochs('xerr_a')
            yerr_a = self.starSet.stars[ind].getArrayAllEpochs('yerr_a')

            xerr = np.sqrt(xerr_p**2 + xerr_a**2)
            yerr = np.sqrt(yerr_p**2 + yerr_a**2)

            try:
                t0 = self.starSet.stars[ind].fitXa.t0    # use the same T0 as polyfit
            except AttributeError:
                continue
            years = np.array(self.starSet.stars[ind].years)

            good = (np.where(abs(x) < 500))[0]

            if (len(good) > 4):
                x = x[good]
                y = y[good]
                xerr = xerr[good]
                yerr = yerr[good]
                years = years[good]

                xfit = self.fitJerk(years - t0, x, xerr)
                yfit = self.fitJerk(years - t0, y, yerr)
                dof = len(x) - len(xfit.params)
                # put the results in the correct star
                try:
                    self.x[i] = xfit.params[0]
                    jx = np.concatenate([jx, [xfit.params[3]]])
                    jxe = np.concatenate([jxe, [xfit.perror[3]]])
                except TypeError:
                    continue
                self.ax[i] = xfit.params[2]
                self.axe[i] = xfit.perror[2]
                self.vx[i] = xfit.params[1]
                self.vxe[i] = xfit.perror[1]

                self.y[i] = yfit.params[0]
                jy = np.concatenate([jy, [yfit.params[3]]])
                jye = np.concatenate([jye, [yfit.perror[3]]])
                self.ay[i] = yfit.params[2]
                self.aye[i] = yfit.perror[2]
                self.vy[i] = yfit.params[1]
                self.vye[i] = yfit.perror[1]
#                pdb.set_trace()

                tmp_chi2xj[i] = xfit.fnorm
                tmp_chi2yj[i] = yfit.fnorm


        self.jx = jx
        self.jxe = jxe
        self.jy = jy
        self.jye = jye
        self.chi2xj = tmp_chi2xj
        self.chi2yj = tmp_chi2yj


    def plotSpeckleAOVel(self):
        """ plot the AO - speckle velocities
        """
        # plot the differences in velocity
        velDiffX = self.velDiffX*1e3
        velDiffY = self.velDiffY*1e3
        starsX = self.x
        starsY = self.y
        
        good = np.where((velDiffX != 0) & (velDiffX < 5) & (velDiffY < 5))[0]
        py.figure(1,figsize=(8,5))
        py.clf()
        py.subplots_adjust(wspace=0.25, right=0.9, left=0.1, top=0.85, bottom=0.15)
        py.subplot(1,2,1)
        qvr = py.quiver(starsX[good],starsY[good],velDiffX[good],velDiffY[good],width=0.03,units='x',angles='xy',scale=1.0)
        py.quiverkey(qvr,5,3,1,'1 mas/yr',coordinates='data')
        py.xlim(6,-6)
        py.ylim(-4,4)
        py.subplot(1,2,2)
        binsIn = np.arange(-20, 20, 0.5)
        (nx, bx, ptx) = py.hist(velDiffX[good],bins=binsIn,color='r',histtype='step',linewidth=1,label='X')
        (ny, by, pty) = py.hist(velDiffY[good],bins=binsIn,color='b',histtype='step',linewidth=1,label='Y')
        py.plot([0,0],[0,120],'k--')
        py.xlabel('Velocity Differences (mas/yr)',fontsize=12)
        py.legend(numpoints=1)
        print 'Median (AO - speckle) velocity (X,Y) = (%6.3f, %6.3f) mas/yr' % \
              (np.median(velDiffX[good]),np.median(velDiffY[good]))
        
    def plotAccelMap(self,save = False):
        """Plot the locations of the unphysical accelerations along
        with the acceleration vector
        """
        bad = self.bad
        goodChi = np.where((self.chi2x < self.chiThreshold*10.0) | (self.chi2y < self.chiThreshold*10))[0]
        bad = np.intersect1d(bad,goodChi)
        
        names = np.array(self.names)[bad]
        x = self.x[bad]
        y = self.y[bad]
        ax = self.ax[bad]*500.0
        ay = self.ay[bad]*500.0
        clf()
        plot(x,y,'o')
        xlim(6,-6)
        xlabel('RA Offset (arcsec)')
        ylabel('DEC Offset (arcsec)')
        quiver(x,y,ax,ay,angles='xy',scale=.7,width=.03,units='x')
        for ii in arange(len(names)):
            text(x[ii],y[ii],names[ii])
        if save:
            savefig(self.plotPrefix+'accelMap.png')
        
    def plotchi2(self, removeConfused = False):
        """Plot the chi-square distribution of the velocity and acceleration fits

        KEYWORDS: removeConfused - remove the confused sources from the analysis
        """
        # get some variables in to local scope
        names = np.array(self.names)
        at, ate = self.at, self.ate
        ar, are = self.ar, self.are
        sigma = self.sigma
        x, y = self.x, self.y
        vx, vy = self.vx, self.vy
        ax, ay = self.ax, self.ay
        epoch, allEpochs = self.epoch, self.allEpochs
        nEpochs = self.nEpochs
        a2d = self.a2d
        mag = self.mag
        cc = self.cc
        chi2x, chi2y = self.chi2x, self.chi2y
        chi2xv, chi2yv = self.chi2xv, self.chi2yv
        r2d = self.r2d
        rarc, a2dsim = self.rarc, self.a2dsim
        bad, goodAccel = self.bad, self.goodAccel
        
        # plot some chi-square info
        clf()
        subplot(231)
        # look at the distribution of chi-squares for stars brighter
        # than 16 magnitude and have maximum degree of freedom
        maxEpoch = np.amax(self.nEpochs)
        dof = maxEpoch - 3
        good = np.where((mag < 16.0) & (nEpochs == maxEpoch))[0]

        # remove confused sources if specified
        if removeConfused:
            good = np.setdiff1d(good, self.confusedSourcesInd)
        
        # filter out the accelerations that do not have the max number of epochs measured
        print shape(good)
        print shape(bad)
        bad = np.intersect1d(good, bad)
        goodAccel = np.intersect1d(goodAccel, good)

        # use mean degree of freedom:
        #dof = np.floor(np.mean(nEpochs)-3)

        print 'Degree of freedom %f' % dof
        maxChi=dof*4.0
        binWidth = (maxChi-0)/30.0
        n, bins, patches1 = hist(chi2x[good], bins = np.arange(0, maxChi, binWidth),normed=1,label='x',alpha=0.6,color='blue')
        n, bins, patches2 = hist(chi2y[good], bins = bins, normed=1, label='y',alpha=0.6,color='green')
        xlabel('Chi-Sq Acc')
        title('All stars')
        chiInd = np.arange(0.01,100,0.01)
        chi2Theory = stats.chi2.pdf(chiInd,dof)
        plot(chiInd, chi2Theory)

        legend()
        subplot(232)
        n, bins, patches1 = hist(chi2x[bad], bins = np.arange(0,maxChi,binWidth), normed=1,label='x',alpha=0.6)
        n, bins, patches2 = hist(chi2y[bad], bins = bins,normed=1,label='y',alpha=0.6)
        title('Non-physical Accelerations')
        xlabel('Chi-Sq Acc')
        plot(chiInd, chi2Theory)

        legend()

        # chi-sq distribution of physical accelerations
        subplot(233)
        n, bins, patches1 = hist(chi2x[goodAccel], bins = np.arange(0,maxChi,binWidth), normed=1,label='x',alpha=0.6)
        n, bins, patches2 = hist(chi2y[goodAccel], bins = bins,normed=1,label='y',alpha=0.6)
        title('Physical Accelerations')
        xlabel('Chi-Sq Acc')
        plot(chiInd, chi2Theory)

        legend()

        #savefig('./chiSq_dist_accel.png')

        #clf()
        # plot the velocities
        velDof = dof + 1
        subplot(234)
        n, bins, patches1 = hist(chi2xv[good], bins = np.arange(0, maxChi, binWidth),normed=1,label='x',alpha=0.6,color='blue')
        n, bins, patches2 = hist(chi2yv[good], bins = bins, normed=1, label='y',alpha=0.6,color='green')
        xlabel('Chi-Sq Vel')
        title('All stars')
        chiInd = np.arange(0.01,100,0.01)
        chi2Theory = stats.chi2.pdf(chiInd,velDof)
        plot(chiInd, chi2Theory)

        legend()
        subplot(235)
        n, bins, patches1 = hist(chi2xv[bad], bins = np.arange(0,maxChi,binWidth), normed=1,label='x',alpha=0.6)
        n, bins, patches2 = hist(chi2yv[bad], bins = bins,normed=1,label='y',alpha=0.6)
        title('Non-physical Accelerations')
        xlabel('Chi-Sq Vel')
        plot(chiInd, chi2Theory)

        legend()

        subplot(236)
        n, bins, patches1 = hist(chi2xv[goodAccel], bins = np.arange(0,maxChi,binWidth), normed=1,label='x',alpha=0.6)
        n, bins, patches2 = hist(chi2yv[goodAccel], bins = bins,normed=1,label='y',alpha=0.6)
        title('Physical Accelerations')
        xlabel('Chi-Sq Vel')
        plot(chiInd, chi2Theory)

        legend()

        savefig(self.plotPrefix+'chiSq_dist_accel_vel.png')

        clf()
        subplot(131)
        loglog(chi2xv[good]/(dof+1.0),chi2x[good]/dof,'bo',label='X',alpha=0.6)
        plot(chi2yv[good]/(dof+1.0),chi2y[good]/dof,'go',label='Y',alpha=0.6)
        plot([0.001,100],[0.001,100])
        xlim(0.01,maxChi)
        ylim(0.01,maxChi)
        xlabel('Vel. Fit Reduced Chi-Sq')
        ylabel('Acc. Fit Reduced Chi-Sq')
        title('All Stars')
        legend(loc=2)

        subplot(132)
        loglog(chi2xv[bad]/(nEpochs[bad]-2.0),chi2x[bad]/(nEpochs[bad]-3.0),'bo',label='X',alpha=0.6)
        plot(chi2yv[bad]/(nEpochs[bad]-2.0),chi2y[bad]/(nEpochs[bad]-3.0),'go',label='Y',alpha=0.6)
        plot([0.001,100],[0.001,100])
        xlim(0.01,maxChi)
        ylim(0.01,maxChi)
        xlabel('Vel. Fit Reduced Chi-Sq')
        ylabel('Acc. Fit Reduced Chi-Sq')
        title('Non-physical Accelerations')
        legend(loc=2)

        subplot(133)
        loglog(chi2xv[goodAccel]/(nEpochs[goodAccel]-2.0),chi2x[goodAccel]/(nEpochs[goodAccel]-3.0),'bo',label='X',alpha=0.6)
        plot(chi2yv[goodAccel]/(nEpochs[goodAccel]-2.0),chi2y[goodAccel]/(nEpochs[goodAccel]-3.0),'go',label='Y',alpha=0.6)
        plot([0.001,100],[0.001,100])
        xlim(0.01,maxChi)
        ylim(0.01,maxChi)
        xlabel('Vel. Fit Reduced Chi-Sq')
        ylabel('Acc. Fit Reduced Chi-Sq')
        title('Physical Accelerations')
        legend(loc=2)

        savefig(self.plotPrefix+'chi2_vel_vs_accel.png')

        # plot chi-sq as a function of magnitude
        clf()
        subplot(231)
        semilogy(mag[good],chi2x[good]/(nEpochs[good]-3.0),'bo',label='X',alpha=0.6)
        plot(mag[good],chi2y[good]/(nEpochs[good]-3.0),'go',label='Y',alpha=0.6)
        ylim(0.01,maxChi)
        xlabel('K mag')
        ylabel('Accel. Fit Chi-Sq')
        title('All Stars')
        legend(loc=3)

        subplot(232)
        semilogy(mag[bad],chi2x[bad]/(nEpochs[bad]-3.0),'bo',label='X',alpha=0.6)
        plot(mag[bad],chi2y[bad]/(nEpochs[bad]-3.0),'go',label='Y',alpha=0.6)
        ylim(0.01,maxChi)
        xlabel('K mag')
        ylabel('Accel. Fit Chi-Sq')
        title('Non-physical Accelerations')
        legend(loc=3)

        subplot(233)
        semilogy(mag[goodAccel],chi2x[goodAccel]/(nEpochs[goodAccel]-3.0),'bo',label='X',alpha=0.6)
        plot(mag[goodAccel],chi2y[goodAccel]/(nEpochs[goodAccel]-3.0),'go',label='Y',alpha=0.6)
        ylim(0.01,maxChi)
        xlabel('K mag')
        ylabel('Accel. Fit Chi-Sq')
        title('Physical Accelerations')
        legend(loc=3)

        subplot(234)
        semilogy(mag[good],chi2xv[good]/(nEpochs[good]-2.0),'bo',label='X',alpha=0.6)
        plot(mag[good],chi2yv[good]/(nEpochs[good]-2.0),'go',label='Y',alpha=0.6)
        ylim(0.01,maxChi)
        xlabel('K mag')
        ylabel('Vel. Fit Chi-Sq')
        title('All Stars')
        legend(loc=3)

        subplot(235)
        semilogy(mag[bad],chi2xv[bad]/(nEpochs[bad]-2.0),'bo',label='X',alpha=0.6)
        plot(mag[bad],chi2yv[bad]/(nEpochs[bad]-2.0),'go',label='Y',alpha=0.6)
        ylim(0.01,maxChi)
        xlabel('K mag')
        ylabel('Vel. Fit Chi-Sq')
        title('Non-physical Accelerations')
        legend(loc=3)

        subplot(236)
        semilogy(mag[goodAccel],chi2xv[goodAccel]/(nEpochs[goodAccel]-2.0),'bo',label='X',alpha=0.6)
        plot(mag[goodAccel],chi2yv[goodAccel]/(nEpochs[goodAccel]-2.0),'go',label='Y',alpha=0.6)
        ylim(0.01,maxChi)
        xlabel('K mag')
        ylabel('Vel. Fit Chi-Sq')
        title('Physical Acceleration')
        legend(loc=3)

        savefig(self.plotPrefix+'mag_vs_chi2.png')



    def plotChi2RadialDependence(self):
        """ Plot the chi2 distribution by splitting the stars into
        different distances from Sgr A* to see if there is a
        dependence in the additive error that we need
        """
    def saveClass(self, filename = False):
        """Save the class with a pickle file
        """
        if not(filename):
            filename = self.prefix+'_accelClass.pkl'
        # open a file
        f = open(filename,'wb')
        print 'saving file: '+filename
        cPickle.dump(self,f,-1)
        f.close()


    
    def newZero(self, mkPointsFiles = False, runPolyfit=False, S02Dir = '/g/ghez/align/13_08_21/efit/chains_S0-2_newRV2/', marg_mean = False):
        """ Loads points after remove confused and corrects for the
        moving zero point given by S0-2 efit

        Will create new points files in points_nz/ with new zero. 

        Keywords: mkPointsFiles - create new points files that have
        the confused sources removed and put them in points_nz/

        runPolyfit - run poly fit on the points_nz directory and put
        polyfit results in polyfit_nz/

        S02Dir - set directory within rootDir of the efit_summary.txt file that will be used

        marg_mean - If true, uses the mean values of the marginalized distribrutions of position and velocity
        as opposed to using the values from the maximum likelihood
        """

        #Values for zero from S0-2 orbit
#        ori_x0= -0.398237949233159853E-02
#        ori_y0= -0.852560358066134505E-02
#        ori_vx= 0.184690930460452768E-03
#        ori_vy= 0.548734449435500692E-03
        ori_t0= 2000. #hard coded in
#        ori_xerr= 0.632109315557945173E-03
#        ori_yerr= 0.136635578849187800E-02
#        ori_vxerr= 0.378284674963935368E-04
#        ori_vyerr= 0.378284674963935368E-04
        origin_val = asciidata.open(S02Dir + 'efit_summary.txt')
        print 'Loaded file'
        if marg_mean:
            ori_x0 = origin_val[1][0]
            ori_y0 = origin_val[2][0]
            ori_vx = origin_val[3][0]
            ori_vy = origin_val[4][0]
        else:
            ori_x0 = origin_val[35][0]
            ori_y0 = origin_val[36][0]
            ori_vx = origin_val[37][0]
            ori_vy = origin_val[38][0]
        print 'Loaded values'
        ori_vxe = origin_val[20][0]
        ori_vye = origin_val[21][0]



        oldDir = os.path.split(self.rootDir+self.points)[0]
        newDir = self.rootDir+'points_nz'
        names=list(self.names)
        nRange=np.arange(len(names))

        if mkPointsFiles:
            # make a new directory and then copy all the old points files there
            if os.path.isdir(newDir):
                print 'removing exisiting points_nz directory: '+newDir
                shutil.rmtree(newDir)

            # copy the files into the new directory
            print 'copying files to :', newDir
            cmdStr = 'cp -rf %s %s' % (oldDir, newDir)
            print 'entering command: '+cmdStr
            os.system(cmdStr)
            print 'done copying files'
            workDir = newDir+'/'
        else:
            #if not making new files then use old directory
            workDir = oldDir+'/'

        for i in nRange:
            ptFile = workDir + names[i] + '.points'
            tab=asciidata.open(ptFile)
            
            time = tab[0].tonumpy()
#            x = tab[1].tonumpy()
#            y = tab[2].tonumpy()
#            xerr = tab[3].tonumpy()
#            yerr = tab[4].tonumpy()

            #Cycle through epochs, update X and Y positions with new origin
            for j in range(len(time)):
                tab[1][j] = tab[1][j] - ((tab[0][j]-ori_t0)*ori_vx + ori_x0)
                tab[2][j] = tab[2][j] - ((tab[0][j]-ori_t0)*ori_vy + ori_y0)

                #Update errors JUST WITH ERRORS IN VELOCITY
 #               tab[3][j] = math.sqrt(tab[3][j]**2 + ((tab[0][j]-ori_t0)*ori_vxe)**2)
 #               tab[4][j] = math.sqrt(tab[4][j]**2 + ((tab[0][j]-ori_t0)*ori_vye)**2)


            #Update points file
            if mkPointsFiles:
                tab.flush()
        
        if runPolyfit:
            print 'now running polyfit in '+self.rootDir+'polyfit_nz'
            #gcutil.mkdir(self.rootDir+'polyfit_c')
            #cmd = 'polyfit -d 2 -linear -i '+self.rootDir+'/align/align_d_rms_1000_abs_t '
            cmd = 'polyfit -d 2 -linear -i '+self.rootDir+self.align
            cmd += ' -tex -plot -points '+self.rootDir+'/points_nz/ -o '+self.rootDir+'/polyfit_nz/fit'
            os.system(cmd)
