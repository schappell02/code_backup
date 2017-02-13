import matplotlib as mtpl
mtpl.use('PDF')
import scipy
import pyfits
from scipy import stats
from scipy import special
from scipy import integrate
import pickle
import nmpfit_sy
import asciidata, os, sys, pickle
import nmpfit_sy2 as nmpfit_sy
import numpy as np
import math
import pdb
import time
import scipy.optimize as opter
from scipy.optimize import fsolve
from scipy.optimize import minimize
from matplotlib.ticker import ScalarFormatter 
import datetime
import time
import threading
import pylab as py
# Several functions were taken from:
# /ghezgroup/code/python/gcwork/polyfit/accel.py,
# velocity.py and residuals.py, written by JLu.


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




def forCluster(num=10,start=1,gamma=1.0,alpha=2.0,delta=2.0,r_break=0.5,cfact=0.3,base_label=''):
    for i in range(num):
        tmplabel = base_label + '__' + str(i+start)
        testCMN(gamma=gamma,alpha=alpha,delta=delta,r_break=r_break,cfact=cfact,label=tmplabel)

def testCMN(numStars=3705,min_r=0.0,max_r=5.0,gamma=-1.0,alpha=4.0,delta=4.0,r_break=0.5,perturb=False,R2d_cut=5.0,
            label='test_r2donly',offset_al=0.0,offset_de=0.0,offset_br=0.0,solver=True,situation=3,innerCut=5.0,outerCut=15.0,
            nonRadial=1,resume=False,cfact=0.5):


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

    if ((resume == False)):
        if ((situation != 1) & (situation != 2) & (situation != 3) & (situation != 4)):
            sys.exit("Did not pick available situation, try 1, 2, 3, or 4")
        

        if (solver==False):
    #if not using solver
            numt = 100000

        R2d_cut *= dist / au_in_pc #R2d_cut now in pc
        innerCut *= dist / au_in_pc
        outerCut *= dist / au_in_pc

    #read in actual stars
        actual = np.loadtxt('gcows_for_mock.dat')
        are_act = actual[:,2]
        pOld_act = actual[:,3]
        
        if (nonRadial>0):
            plate_scale=0.00995 #arcsec/pixel
            gcows_zero = np.array([1500.0,1500.0])
            gcows = pyfits.getdata('nirc2_gcows_2010_all_mask.fits')

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
                    #print ii
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
                        #print ii
                        ii += 1

            #else:
            if ((rand_R2d < outerCut) & (rand_R2d > innerCut)):
                if ((situation > 2) & (cval > cfact)):
                    radiusM = np.append(radiusM,rand_r)
                    R2dM = np.append(R2dM,rand_R2d)
                    randFinalM = np.append(randFinalM,rand01)
                    #print ii
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
        np.savetxt('stars_mn'+label+'.dat',np.transpose([R2d,accel,are,pOld]),delimiter=' ')

        if (situation > 2):
            radiusM *= cm_in_pc
            R2dM *= cm_in_pc #in cm
        #read in actual stars
        #actual = np.loadtxt('/u/schappell/code/c/maser_mn_actual.dat')
            actual = np.loadtxt('schodel_for_mock.dat')
            pOldM_act = actual[:,1]
            randex = np.random.randint(len(pOldM_act),size=len(radiusM))
            pOldM = pOldM_act[randex]
            np.savetxt('maser_mn'+label+'.dat',np.transpose([R2dM,pOldM]),delimiter=' ')
            radius = np.append(radius,radiusM)
            radius2d = np.append(R2d,R2dM)

    #pdb.set_trace()
        np.savetxt('stars_r'+label+'.dat',np.transpose([radius2d,radius]),delimiter=' ')
        #np.savetxt('/u/schappell/code/c/xy_pixels'+label+'.dat',np.transpose([xPIX,yPIX]),delimiter=' ')




#def compareMulti(flagArray=['g1_a2_d2_b0.5_c0.3'],legendArray=['Mock 1-2-2-0.5-0.3'],priorArray=['priors.txt']):
    
 #   for i in range(len(flagArray)):
        

        

def plotMockParam(mockTag='g1_a2_d2_b0.5_c0.3',priorFile='/u/schappell/code/cluster/priors.txt',numM=100,Cplot=True,
                  r_gamma=1.0,r_alpha=2.0,r_delta=2.0,r_rbreak=0.5,r_cfact=0.3,mockTag_p='g1_a2_d2_b0.5_c0.3__'):
    priors = np.loadtxt(priorFile)
    gamma = np.array([])
    alpha = np.array([])
    delta = np.array([])
    rbreak = np.array([])
    if (Cplot==True):
        cfact = np.array([])
    weights = np.array([])
    for i in range(100):
        i += 1
        tmpMock = np.loadtxt('/u/schappell/code/cluster/results_'+mockTag+'/'+mockTag_p+str(i)+'.txt')
        tmpgamma = tmpMock[:,2]*(priors[0,1] - priors[0,0]) + priors[0,0]
        tmpalpha = tmpMock[:,3]*(priors[1,1] - priors[1,0]) + priors[1,0]
        tmpdelta = tmpMock[:,4]*(priors[2,1] - priors[2,0]) + priors[2,0]
        tmpbreak = tmpMock[:,5]*(priors[3,1] - priors[3,0]) + priors[3,0]
        if (Cplot==True):
            tmpcfact = tmpMock[:,6]
            cfact = np.append(cfact, tmpcfact)
        tmpweights = tmpMock[:,0]/np.sum(tmpMock[:,0])

        gamma = np.append(gamma, tmpgamma)
        alpha = np.append(alpha, tmpalpha)
        delta = np.append(delta, tmpdelta)
        rbreak = np.append(rbreak, tmpbreak)
        weights = np.append(weights,tmpweights)

    py.clf()
    hist,bins,junk=py.hist(gamma,bins=numM,normed=1,weights=weights)
    py.xlabel('$\gamma$ (Inner Slope)')
    py.ylabel('Weighted Posterior')
    print 'Gamma --------------------------'
    print 'Mean: '+str(np.average(gamma,weights=weights))+' +/- '+str(math.sqrt(np.average((gamma-np.average(gamma,weights=weights))**2,weights=weights)))
    print 'Median: '+str(np.median(gamma))
    tmpdex = np.argmax(hist)
    print 'Max bin: '+str(bins[tmpdex])+' to '+str(bins[tmpdex+1])
    levels = getContourLevels(hist)#,percLevels=[0.999999998027,0.999999426697,0.999936657516])
    tmpdex=np.where(hist > levels[2])[0]
    ax = py.gca()
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
    py.plot([r_gamma,r_gamma],ylim,color='r',linewidth=4.0,linestyle='--')
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+mockTag+'_gamma_post.png')
    py.clf()

    hist,bins,junk=py.hist(alpha,bins=numM,normed=1,weights=weights)
    py.xlabel(r'$\alpha$ (Outer Slope)')
    py.ylabel('Weighted Posterior')
    ax = py.gca()
    ylim = ax.get_ylim()
    py.plot([r_alpha,r_alpha],ylim,color='r',linewidth=4.0,linestyle='--')
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
    py.savefig('/u/schappell/plots/'+mockTag+'_alpha_post.png')
    py.clf()

    hist,bins,junk=py.hist(delta,bins=numM,normed=1,weights=weights)
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel('Weighted Posterior')
    ax = py.gca()
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
    py.plot([r_delta,r_delta],ylim,color='r',linewidth=4.0,linestyle='--')
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+mockTag+'_delta_post.png')
    py.clf()

    hist,bins,junk=py.hist(rbreak,bins=numM,normed=1,weights=weights)
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel('Weighted Posterior')
    ax = py.gca()
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
    py.plot([r_rbreak,r_rbreak],ylim,color='r',linewidth=4.0,linestyle='--')
    py.ylim(ylim)
    py.savefig('/u/schappell/plots/'+mockTag+'_rbreak_post.png')
    py.clf()

    if (Cplot==True):
        hist,bins,junk=py.hist(cfact,bins=4*numM,normed=1,weights=weights)
        py.xlabel('C Factor')
        py.ylabel('Weighted Posterior')
        ax = py.gca()
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
        py.plot([r_cfact,r_cfact],ylim,color='r',linewidth=4.0,linestyle='--')
        py.ylim(ylim)
        py.savefig('/u/schappell/plots/'+mockTag+'_cfactor_post.png')
        py.clf()


    #2D contours
    py.clf()
    hist,ybins,xbins = np.histogram2d(gamma,alpha,bins=numM/4,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
    py.plot(r_alpha,r_gamma,marker='d',ms=15,color='r')
    py.xlabel(r'$\alpha$ (Outer Slope)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/'+mockTag+'_gamma_alpha.png')
    print 'GAMMA AND ALPHA --------------------------'
    tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

    hist,ybins,xbins = np.histogram2d(gamma,delta,bins=numM/4,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
    py.plot(r_delta,r_gamma,marker='d',ms=15,color='r')
    py.xlabel(r'$\delta$ (Sharpness)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/'+mockTag+'_gamma_delta.png')
    print 'GAMMA AND DELTA --------------------------'
    tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

    hist,ybins,xbins = np.histogram2d(gamma,rbreak,bins=numM/4,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.plot(r_rbreak,r_gamma,marker='d',ms=15,color='r')
    py.xlabel(r'$r_{break}$ (pc)')
    py.ylabel('$\gamma$ (Inner Slope)')
    py.savefig('/u/schappell/plots/'+mockTag+'_gamma_rbreak.png')
    print 'GAMMA AND r_BREAK --------------------------'
    tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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
        hist,ybins,xbins = np.histogram2d(gamma,cfact,bins=numM/4,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
        #        py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.plot(r_cfact,r_gamma,marker='d',ms=15,color='r')
        py.xlabel('C Factor')
        py.ylabel('$\gamma$ (Inner Slope)')
        py.savefig('/u/schappell/plots/'+mockTag+'_gamma_cfactor.png')
        print 'GAMMA AND C FACTOR --------------------------'
        tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

        hist,ybins,xbins = np.histogram2d(alpha,cfact,bins=numM/4,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.plot(r_cfact,r_alpha,marker='d',ms=15,color='r')
        py.ylabel(r'$\alpha$ (Outer Slope)')
        py.savefig('/u/schappell/plots/'+mockTag+'_alpha_cfactor.png')
        print 'ALPHA AND C FACTOR --------------------------'
        tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

        hist,ybins,xbins = np.histogram2d(delta,cfact,bins=numM/4,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.plot(r_cfact,r_delta,marker='d',ms=15,color='r')
        py.ylabel(r'$\delta$ (Sharpness)')
        py.savefig('/u/schappell/plots/'+mockTag+'_delta_cfactor.png')
        print 'DELTA AND C FACTOR --------------------------'
        tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

        hist,ybins,xbins = np.histogram2d(rbreak,cfact,bins=numM/4,weights=weights)
        x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
        y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
        levels = getContourLevels(hist)
        X,Y=np.meshgrid(x,y)
        py.clf()
        py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
        py.xlabel('C Factor')
        py.plot(r_cfact,r_rbreak,marker='d',ms=15,color='r')
        py.ylabel(r'$r_{break}$ (prc)')
        py.savefig('/u/schappell/plots/'+mockTag+'_rbreak_cfactor.png')
        print 'r_BREAK AND C FACTOR --------------------------'
        tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

    hist,ybins,xbins = np.histogram2d(alpha,delta,bins=numM/4,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$\delta$ (Sharpness)')
    py.plot(r_delta,r_alpha,marker='d',ms=15,color='r')
    py.ylabel(r'$\alpha$ (Outer Slope)')
    py.savefig('/u/schappell/plots/'+mockTag+'_alpha_delta.png')
    print 'ALPHA AND DELTA --------------------------'
    tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

    hist,ybins,xbins = np.histogram2d(alpha,rbreak,bins=numM/4,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$r_{break}$ (pc)')
    py.plot(r_rbreak,r_alpha,marker='d',ms=15,color='r')
    py.ylabel(r'$\alpha$ (Outer Slope)')
    py.savefig('/u/schappell/plots/'+mockTag+'_alpha_rbreak.png')
    print 'ALPHA AND r_BREAK --------------------------'
    tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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

    hist,ybins,xbins = np.histogram2d(delta,rbreak,bins=numM/4,weights=weights)
    x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
    y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
    levels = getContourLevels(hist)
    X,Y=np.meshgrid(x,y)
    py.clf()
    py.contour(X,Y,hist,levels=levels, colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid',))
#py.contour(X,Y,hist,levels=levels, colors=('limegreen','k','gray'),linestyles=('solid','dashed','dotted'))
    py.xlabel(r'$r_{break}$ (pc)')
    py.plot(r_rbreak,r_delta,marker='d',ms=15,color='r')
    py.ylabel(r'$\delta$ (Sharpness)')
    py.savefig('/u/schappell/plots/'+mockTag+'_delta_rbreak.png')
    print 'DELTA AND r_BREAK --------------------------'
    tmpdex = np.unravel_index(hist.argmax(),hist.shape)
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


    #Make pyramid plot with 1D posteriors
    py.clf()
    #left,bottom,width,height,space = 0.12,0.12,0.2,0.2,0.02
    Xarray = np.zeros((10,numM/4,numM/4))
    Yarray = np.zeros((10,numM/4,numM/4))
    histarray = np.zeros((10,numM/4,numM/4))
    ij = 0
    if (Cplot==True):
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
                hist,ybins,xbins = np.histogram2d(yarray,xarray,bins=numM/4,weights=weights)
                x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
                y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
                X,Y=np.meshgrid(x,y)
                Xarray[ij,:,:] = X
                Yarray[ij,:,:] = Y
                histarray[ij,:,:] = hist
                ij += 1

        ij = 0
        xticks = np.array([2.0,4.0,0.75,0.3,4.0])#,0.1])
        x_real=np.array([r_gamma,r_delta,r_rbreak,r_cfact])
        yticks = np.array([4.0,0.3,1.0,4.0])
        y_real=np.array([r_alpha,r_cfact,r_rbreak,r_delta])
        py.figure(1,figsize=(8.0,8.0))
        for i in range(5):
            for j in range(5-i):
                py.axes([0.07+(j*0.175),0.084+(i*0.18),0.165,0.165])
                if ((i+j) < 4):
                    levels = getContourLevels(histarray[ij,:,:])
                    py.contour(Xarray[ij,:,:],Yarray[ij,:,:],histarray[ij,:,:],levels=levels, 
                               colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid'))
                    py.plot(x_real[j],y_real[i],marker='d',ms=7,color='r')
                    ij += 1
                    ax = py.gca()
                    ax.xaxis.set_major_locator(mtpl.ticker.MultipleLocator(xticks[j]))
                    ax.yaxis.set_major_locator(mtpl.ticker.MultipleLocator(yticks[i]))
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
                    ax = py.gca()
                    yticks1 = np.array([0.1,10.0,1.0,0.06,1.0])
                    ax.xaxis.set_major_locator(mtpl.ticker.MultipleLocator(xticks[j]))
                    ax.yaxis.set_major_locator(mtpl.ticker.MultipleLocator(yticks1[i]))
                    if (j > 0):
                        ax.axes.yaxis.tick_right()
                        ax.yaxis.set_label_position('right')
                    if (j < 4):
                        ax.axes.xaxis.set_ticklabels([])
                    if (j == 0):
                        hist,bins,junk=py.hist(gamma,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_gamma,r_gamma],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(gxlim)
                        py.ylabel('Posterior')
                    elif (j == 1):
                        hist,bins,junk=py.hist(delta,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_delta,r_delta],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(dxlim)
                    elif (j == 2):
                        hist,bins,junk=py.hist(rbreak,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_rbreak,r_rbreak],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(rxlim)
                    elif (j == 3):
                        hist,bins,junk=py.hist(cfact,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_cfact,r_cfact],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(cxlim)
                    elif (j == 4):
                        hist,bins,junk=py.hist(alpha,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        py.xlabel(r'$\alpha$ (Outer Slope)')
                        ylim=ax.get_ylim()
                        py.plot([r_alpha,r_alpha],ylim,color='r')
                        py.ylim(ylim)
                        ytmp = ax.yaxis.get_major_ticks()
                        ytmp[0].label1.set_visible(False)
                        xtmp = ax.yaxis.get_major_ticks()
                        xtmp[0].label1.set_visible(False)



    else:
        for i in range(3):
            if (i==0):
                yarray = alpha * 1.0
            elif (i==1):
                yarray = rbreak * 1.0
            else:
                yarray = delta * 1.0
            for j in range(3-i):
                if (j==0):
                    xarray = gamma * 1.0
                elif (j==1):
                    xarray = delta * 1.0
                else:
                    xarray = rbreak * 1.0
                hist,ybins,xbins = np.histogram2d(yarray,xarray,bins=numM/4,weights=weights)
                x = np.array([(xbins[i]+xbins[i+1])/2.0 for i in range(numM/4)])
                y = np.array([(ybins[i]+ybins[i+1])/2.0 for i in range(numM/4)])
                X,Y=np.meshgrid(x,y)
                Xarray[ij,:,:] = X
                Yarray[ij,:,:] = Y
                histarray[ij,:,:] = hist
                ij += 1

        ij = 0
        xticks = np.array([2.0,4.0,0.75,4.0])#,0.1])
        x_real=np.array([r_gamma,r_delta,r_rbreak])
        yticks = np.array([4.0,1.0,4.0])
        y_real=np.array([r_alpha,r_cfact,r_rbreak,r_delta])
        py.figure(1,figsize=(8.0,8.0))
        for i in range(4):
            for j in range(4-i):
                py.axes([0.06+(j*0.225),0.084+(i*0.225),0.21,0.21])
                if ((i+j) < 3):
                    levels = getContourLevels(histarray[ij,:,:])
                    py.contour(Xarray[ij,:,:],Yarray[ij,:,:],histarray[ij,:,:],levels=levels, 
                               colors=('gray','k','limegreen'),linestyles=('dotted','dashed','solid'))
                    py.plot(x_real[j],y_real[i],marker='d',ms=7,color='r')
                    ij += 1
                    ax = py.gca()
                    ax.xaxis.set_major_locator(mtpl.ticker.MultipleLocator(xticks[j]))
                    ax.yaxis.set_major_locator(mtpl.ticker.MultipleLocator(yticks[i]))
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
                    if (j == 0):
                        if (i==0):
                            py.ylabel(r'$\alpha$ (Outer Slope)')
                        elif (i==1):
                            py.ylabel(r'$r_{break}$ (pc)')
                        elif (i==2):
                            py.ylabel(r'$\delta$ (Sharpness)')

                else:
                    ax = py.gca()
                    yticks1 = np.array([0.1,1.0,0.06,1.0])
                    ax.xaxis.set_major_locator(mtpl.ticker.MultipleLocator(xticks[j]))
                    ax.yaxis.set_major_locator(mtpl.ticker.MultipleLocator(yticks1[i]))
                    if (j > 0):
                        ax.axes.yaxis.tick_right()
                        ax.yaxis.set_label_position('right')
                    if (j < 3):
                        ax.axes.xaxis.set_ticklabels([])
                    if (j == 0):
                        hist,bins,junk=py.hist(gamma,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_gamma,r_gamma],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(gxlim)
                        py.ylabel('Posterior')
                    elif (j == 1):
                        hist,bins,junk=py.hist(delta,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_delta,r_delta],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(dxlim)
                    elif (j == 2):
                        hist,bins,junk=py.hist(rbreak,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        ylim=ax.get_ylim()
                        py.plot([r_rbreak,r_rbreak],ylim,color='r')
                        py.ylim(ylim)
                        py.xlim(rxlim)
                    elif (j == 3):
                        hist,bins,junk=py.hist(alpha,weights=weights,normed=1,bins=100,histtype='step',color='b')
                        py.xlabel(r'$\alpha$ (Outer Slope)')
                        ylim=ax.get_ylim()
                        py.plot([r_alpha,r_alpha],ylim,color='r')
                        py.ylim(ylim)
                        ytmp = ax.yaxis.get_major_ticks()
                        ytmp[0].label1.set_visible(False)
                        xtmp = ax.yaxis.get_major_ticks()
                        xtmp[0].label1.set_visible(False)
    py.savefig('/u/schappell/plots/'+mockTag+'_pyramid.eps',format='eps',dpi=1000)




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
