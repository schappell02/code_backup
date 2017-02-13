from func_ellipse import *
import asciidata
import numpy as np
import os
import pylab
from sim_AAS12 import getContourLevels
import multiprocessing
import time

def plot_param_hist(mnestfile,star='',color='black',star_num=0,
                    a_lim=None,m_lim=None,p_lim=None,to_lim=None,
                    e_lim=None,i_lim=None,o_lim=None,O_lim=None,
                    linestyle='solid',BH_fixed=False,allparamsinput=False,
                    paramfile=None,nbins=40,extmass=False,norange=False,
                    new=True):
    ''' Use this to plot multinest fitting outputs! '''


    def plotProb(var, label, suffix, bestFit, weights, BH_fixed, norange):   # can probably remove the BH_fixed case from this function
        # Compute the probability distribution
        # (not the probability density function)               

        #(prob, bins) = matplotlib.mlab.hist(var, bins=40, normed=False)
        prob, bins, patch = pylab.hist(var, bins=nbins, normed=False, weights=weights, histtype='step', color=color, linestyle=linestyle)
        
        if (BH_fixed and (suffix == 'D' or suffix == 'Vx'  or suffix == 'Vy'  or suffix == 'Vz'  or suffix == 'mass'  or suffix == 'y0'  or suffix == 'x0')) or len(np.unique(var)) <= 1 or norange:
            print 'not calculating range for '+suffix
        else:
            prob = np.array(prob, dtype=float) / prob.sum() # normalize
        
            # Calculate the peak of the probability distribution
            # and the confidence intervals from the 1D Probs.
            sid = (prob.argsort())[::-1]  #  indices for a reverse sort
            probSort = prob[sid]

            peakPix = sid[0]
            peakVal = bins[peakPix]
            peakProb = prob[peakPix]

            # Make a cumulative distribution function starting from the
            # highest pixel value. This way we can find the level above
            # which 68% of the trials will fall.
            cdf = np.cumsum(probSort)
        
            # Determine point at which we reach XX confidence
            idx1 = (np.where(cdf > 0.6827))[0] # 1 sigma
            idx2 = (np.where(cdf > 0.9545))[0] # 2 sigma
            idx3 = (np.where(cdf > 0.9973))[0] # 3 sigma

            if ((len(idx1) < 2) or (len(idx2) < 2) or (len(idx3) < 2)):
                pylab.clf()
                pylab.hist(var,weights=weights)
                print 'Min, Max = ', var.min(), var.max()
                print idx1
                print idx2
                print idx3
        
            level1 = probSort[idx1[0]]
            level2 = probSort[idx2[0]]
            level3 = probSort[idx3[0]]


            # Find the range of values 
            idx1 = (np.where((prob > level1)))[0]
            idx2 = (np.where((prob > level2)))[0]
            idx3 = (np.where((prob > level3)))[0]

            # Parameter Range:
            range1 = np.array([ bins[idx1[0]], bins[idx1[-1]] ])
            range2 = np.array([ bins[idx2[0]], bins[idx2[-1]] ])
            range3 = np.array([ bins[idx3[0]], bins[idx3[-1]] ])

            # Plus/Minus Errors:
            pmErr1 = np.abs(range1 - peakVal)
            pmErr2 = np.abs(range2 - peakVal)
            pmErr3 = np.abs(range3 - peakVal)

            pmErr1_best = np.abs(range1 - bestFit)
            pmErr2_best = np.abs(range2 - bestFit)
            pmErr3_best = np.abs(range3 - bestFit)
        

            # Print it out to the screen:
            print ''
            print 'Best Fit vs. Peak of Prob. Dist. for the %s' % label
            print '   %6s = %f   vs.   %f' % (suffix, bestFit, peakVal)
            print '1, 2, 3 Sigma Confidence Intervals for the %s' % label
            print '   68.27%% = [%10.4f -- %10.4f] or -/+ [%10.4f, %10.4f] [%10.4f, %10.4f]' % \
                (range1[0], range1[1], pmErr1_best[0], pmErr1_best[1], pmErr1[0], pmErr1[1])
            print '   95.45%% = [%10.4f -- %10.4f] or -/+ [%10.4f, %10.4f] [%10.4f, %10.4f]' % \
                (range2[0], range2[1], pmErr2_best[0], pmErr2_best[1], pmErr2[0], pmErr2[1])
            print '   99.73%% = [%10.4f -- %10.4f] or -/+ [%10.4f, %10.4f] [%10.4f, %10.4f]' % \
                (range3[0], range3[1], pmErr3_best[0], pmErr3_best[1], pmErr3[0], pmErr3[1])

            # Write in an output file:
            _out.write('%6s  %10.4f  %10.4f    ' % (suffix, bestFit, peakVal))
            _out.write('%10.4f %10.4f  %10.4f %10.4f    ' % \
                           (range1[0], range1[1], pmErr1[0], pmErr1[1]))
            _out.write('%10.4f %10.4f  %10.4f %10.4f    ' % \
                           (range2[0], range2[1], pmErr2[0], pmErr2[1]))
            _out.write('%10.4f %10.4f  %10.4f %10.4f\n' % \
                           (range3[0], range3[1], pmErr3[0], pmErr3[1]))
       
        #pylab.clf()
        #(pbins, pprob) = histNofill.convertForPlot_new(bins, prob)
        #plot(pbins, pprob, color='black')
        ##plot(bins, prob, color='black')
        pylab.xlabel(label)
        pylab.ylabel('Probability')
    
        # Plot the best-fit value
        #quiver([bestFit], [peakProb * 1.1], [0], [-peakProb*0.1])

        if (suffix == 't0'):
            pylab.gca().get_xaxis().set_major_formatter(pylab.FormatStrFormatter('%.2f'))
            
        #savefig('plotParamsProb_' + suffix + '.png')
        #savefig('plotParamsProb_' + suffix + '.eps')
        
        #return (suffix, label)
        #return (bins, prob)



    # read in all the arrays if the BH is not fixed
    if not BH_fixed:
        print '****************************'
        print 'BH parameters were NOT fixed'
        print '****************************'

        print 'Opening file...'
        inFile = np.loadtxt(mnestfile)

        weights = inFile[:,0]
        print 'Putting samples in arrays...'

        if extmass:
            skip = 1            
        else:
            skip = 0


        mass = inFile[:,2]
        xo = -inFile[:,3]
        yo = inFile[:,4]
        Vx = -inFile[:,5]
        Vy = inFile[:,6]
        Vz = inFile[:,7]
        D = inFile[:,8]
        Omega = inFile[:,9+star_num*6+skip]
        omega = inFile[:,10+star_num*6+skip]
        i = inFile[:,11+star_num*6+skip]
        P = inFile[:,12+star_num*6+skip]
        To = inFile[:,13+star_num*6+skip]
        e = inFile[:,14+star_num*6+skip]

        
    elif BH_fixed and allparamsinput:
        print '*********************************************'
        print 'BH parameters were fixed, allparamsinput=True'
        print '*********************************************'

        print 'Opening file...'
        inFile = np.loadtxt(mnestfile)

        weights = inFile[:,0]
        print 'Putting samples in arrays...'

        if extmass:
            skip = 1
        else:
            skip = 0

        mass = inFile[:,8]
        xo = -inFile[:,9]
        yo = inFile[:,10]
        Vx = -inFile[:,11]
        Vy = inFile[:,12]
        Vz = inFile[:,13]
        D = inFile[:,14]
        Omega = inFile[:,2+star_num*6+skip]
        omega = inFile[:,3+star_num*6+skip]
        i = inFile[:,4+star_num*6+skip]
        P = inFile[:,5+star_num*6+skip]
        To = inFile[:,6+star_num*6+skip]
        e = inFile[:,7+star_num*6+skip]

    # if BH_fixed, then read the weight file (output from master_wrapper)
    # save the rest of the files for later, when they are being plotted
    else:
        print '************************'
        print 'BH parameters were FIXED'
        print '************************'

        print 'Reading weights file...'
        file_weights = open(mnestfile[0:-4]+'_weights.txt')
        if 'new' in mnestfile or new:
            weights = np.array(file_weights.readlines(),dtype=float)
        else:
            weights = np.array(file_weights.readline().split(),dtype=float)
        file_weights.close()


        #print 'Reading mass file...'
        #file_mass = open(mnestfile[0:-4]+'_mass.txt')
        #mass = np.array(file_mass.readline().split(),dtype=float)
        #file_mass.close()

        #print 'Reading xo file...'
        #file_xo = open(mnestfile[0:-4]+'_xo.txt')
        #xo = -np.array(file_xo.readline().split(),dtype=float)
        #file_xo.close()

        #print 'Reading yo file...'
        #file_yo = open(mnestfile[0:-4]+'_yo.txt')
        #yo = np.array(file_yo.readline().split(),dtype=float)
        #file_yo.close()

        #print 'Reading Vx file...'
        #file_Vx = open(mnestfile[0:-4]+'_Vx.txt')
        #Vx = -np.array(file_Vx.readline().split(),dtype=float)
        #file_Vx.close()

        #print 'Reading Vy file...'
        #file_Vy = open(mnestfile[0:-4]+'_Vy.txt')
        #Vy = np.array(file_Vy.readline().split(),dtype=float)
        #file_Vy.close()

        #print 'Reading Vz file...'
        #file_Vz = open(mnestfile[0:-4]+'_Vz.txt')
        #Vz = np.array(file_Vz.readline().split(),dtype=float)
        #file_Vz.close()

        #print 'Reading D file...'
        #file_D = open(mnestfile[0:-4]+'_D.txt')
        #D = np.array(file_D.readline().split(),dtype=float)
        #file_D.close()

   

    axisLabel_x0 = 'Sgr A* $\Delta$RA Position (mas)'
    axisLabel_y0 = 'Sgr A* Y Position (mas)'
    axisLabel_D = 'Ro (kpc)'
    axisLabel_P = 'Period (yr)'
    axisLabel_a = 'Semi-Major Axis (mas)' 
    axisLabel_e = 'Eccentricity'
    axisLabel_To = 'Epoch of Periapse (yr)'
    axisLabel_o = 'Argument of Periapse (deg)'
    axisLabel_i = 'Inclination (deg)'
    axisLabel_O = 'Angle of the Ascending Node (deg)'
    axisLabel_Vx = 'Sgr A* $\Delta$RA Velocity (mas/yr)'
    axisLabel_Vy = 'Sgr A* Y Velocity (mas/yr)'
    axisLabel_Vz = 'Sgr A* Z Velocity (km/s)'
    axisLabel_mass = 'Mass (million solar masses)' 

    # file to write out best fit and other stuff
    _out = open('plotParamsProb_limits'+star+'.txt', 'w')
    _out.write('%6s  %10s  %10s    ' % ('#Param', 'BestFit', 'PeakProb'))
    _out.write('%10s %10s   %10s %10s   ' % \
               ('1sigma_lo', '1sigma_hi', '1sig(-)', '1sig(+)'))
    _out.write('%10s %10s   %10s %10s   ' % \
               ('2sigma_lo', '2sigma_hi', '2sig(-)', '2sig(+)'))
    _out.write('%10s %10s   %10s %10s\n' % \
               ('3sigma_lo', '3sigma_hi', '3sig(-)', '3sig(+)'))

    pylab.rc('axes', titlesize=10, labelsize=10)
    pylab.rc('xtick', labelsize=8)
    pylab.rc('ytick', labelsize=8)



    # open the param file to get the best fit values from the multinest output
    if paramfile:
        table = asciidata.open(paramfile)

        # Make things into arrays of floats, etc.
        D_fit = float( table[0][0] )  # in pc
        x0_fit = -float( table[1][0] )  # in pix (working on abs, scale = 1)
        y0_fit = float( table[2][0] )  # in pix (working on abs, scale = 1)
        a_fit = float( table[3][0] )   # in mas
        P_fit = float( table[4][0] )   # in yrs
        e_fit = float( table[5][0] )
        To_fit = float( table[6][0] )
        Omega_fit = float( table[7][0] )
        i_fit = float( table[8][0] )
        omega_fit = float( table[9][0] )
        Vx_fit = -float( table[10][0] )  # mas/yr
        Vy_fit = float( table[11][0] )  # mas/yr
        Vz_fit = float( table[12][0] )  # km/s
        
        # convert semi-major axis and period into mass
        mass_fit = (a_fit * D_fit / 1000.0)**3 / P_fit**2

        # calculate periapse distance (in mpc)
        pdist_fit = a_fit * D_fit * (1.0 - e_fit) / (206265.)

        # calculate density (in Mo/pc^3)
        density_fit = mass_fit / ((4.0/3.0) * math.pi * pdist_fit**3)

        # Set up axis labels and units
        x0_fit *= 1000.0  # in mas (assumed scale = 1)
        y0_fit *= 1000.0  # in mas (assumed scale = 1)
        D_fit /= 1000.0  # in kpc
        mass_fit /= 1e6      # in millions of solar masses
        density_fit /= 1e6   # in 10^15 Mo/pc^3

    else:
        x0_fit,y0_fit,D_fit,Vx_fit,Vy_fit,Vz_fit,mass_fit,P_fit,a_fit,To_fit,e_fit,i_fit,omega_fit,Omega_fit = 0,0,0,0,0,0,0,0,0,0,0,0,0,0



    # if BH is not fixed, then plot the BH param PDFs
    # if BH_fixed=True, only plot the orbital parameter PDFs (that is done below)
    # all arrays were read in before!
    if not BH_fixed:
        print '****************************'
        print 'BH parameters were NOT fixed'
        print '****************************'
        pylab.figure(2, figsize=(11, 13))
        pylab.subplots_adjust(bottom=0.05, left=0.09, right=0.95, top=0.97)

        print '\nplotting xo'
        xo *= 1000
        pylab.subplot(5,3,1)
        print np.unique(xo)
        print np.unique(yo)
        print np.unique(D)
        print np.unique(Vx)
        print np.unique(Vy)
        print np.unique(Vz)
        print np.unique(mass)
        plotProb(xo, axisLabel_x0, 'x0', x0_fit, weights, BH_fixed, norange)

        print '\nplotting yo'
        yo *= 1000
        pylab.subplot(5,3,2)
        plotProb(yo, axisLabel_y0, 'y0', y0_fit, weights, BH_fixed, norange)

        print '\nplotting D'
        D /= (1.e3)
        pylab.subplot(5,3,3)
        plotProb(D, axisLabel_D, 'D', D_fit, weights, BH_fixed, norange)
    
        print '\nplotting Vx'
        Vx *= 1000
        pylab.subplot(5,3,4)
        plotProb(Vx, axisLabel_Vx, 'Vx', Vx_fit, weights, BH_fixed, norange)

        print '\nplotting Vy'
        Vy *= 1000
        pylab.subplot(5,3,5)
        plotProb(Vy, axisLabel_Vy, 'Vy', Vy_fit, weights, BH_fixed, norange)

        print '\nplotting Vz'
        pylab.subplot(5,3,6)
        plotProb(Vz, axisLabel_Vz, 'Vz', Vz_fit, weights, BH_fixed, norange)
    

        #a_AU = (mass*P**2.)**(1./3.)   # a in AU = (Mass*Period^2)^(1/3)
        #a = (a_AU/D)*(1.e3)  # a in mas  

        #pylab.subplot(5,3,7)
        #plotProb(a, axisLabel_a, 'a', a_fit, weights, BH_fixed)
        #if a_lim:
        #    pylab.xlim(a_lim)

        print '\nplotting mass'
        mass /= (1.e6)
        pylab.subplot(5,3,9)
        if not m_lim:
            m_lim = [mass.min(),mass.max()]
        lim_idx = np.where((mass>=m_lim[0]) & (mass<=m_lim[1]))
        plotProb(mass[lim_idx], axisLabel_mass, 'mass', mass_fit, weights[lim_idx], BH_fixed, norange)
        if m_lim:
            pylab.xlim(m_lim)


    # in any case, plot the orbital parameter PDFs
    else:
        pylab.figure(2, figsize=(11, 5))


    # Period
    if not BH_fixed:
        pylab.subplot(5,3,10)
    else:
        pylab.subplot(2,3,1)

        if not allparamsinput:
            print 'Reading P file...'
            file_P = open(mnestfile[0:-4]+'_P.txt')
            if 'new' in mnestfile or new:
                P = np.array(file_P.readlines(),dtype=float)
            else:
                P = np.array(file_P.readline().split(),dtype=float)
            file_P.close()

    print '\nplotting P'
    if not p_lim:
        p_lim = [P.min(),P.max()]
    lim_idx = np.where((P>=p_lim[0]) & (P<=p_lim[1]))
    print len(np.unique(P))
    plotProb(P[lim_idx], axisLabel_P, 'P', P_fit, weights[lim_idx], BH_fixed, norange)
    if p_lim:
        pylab.xlim(p_lim)
    del P


    # To
    if not BH_fixed:
        pylab.subplot(5,3,11)
    else:
        pylab.subplot(2,3,2)

        if not allparamsinput:
            print 'Reading To file...'
            file_To = open(mnestfile[0:-4]+'_To.txt')
            if 'new' in mnestfile or new:
                To = np.array(file_To.readlines(),dtype=float)
            else:
                To = np.array(file_To.readline().split(),dtype=float)
            file_To.close()

    print '\nplotting To'
    if not to_lim:
        to_lim = [To.min(),To.max()]
    lim_idx = np.where((To>=to_lim[0]) & (To<=to_lim[1]))
    plotProb(To[lim_idx], axisLabel_To, 't0', To_fit, weights[lim_idx], BH_fixed, norange)
    if to_lim:
        pylab.xlim(to_lim)
    to_axes = pylab.gca()
    xticks=to_axes.get_xticks()
    to_axes.set_xticklabels(np.array(xticks,dtype='string'))
    del To


    # Eccentricity
    if not BH_fixed:
        pylab.subplot(5,3,12)
    else:
        pylab.subplot(2,3,3)

        if not allparamsinput:
            print 'Reading e file...'
            file_e = open(mnestfile[0:-4]+'_e.txt')
            if 'new' in mnestfile or new:
                e = np.array(file_e.readlines(),dtype=float)
            else:
                e = np.array(file_e.readline().split(),dtype=float)
            file_e.close()

    print '\nplotting e'
    if not e_lim:
        e_lim = [e.min(),e.max()]
    lim_idx = np.where((e>=e_lim[0]) & (e<=e_lim[1]))
    plotProb(e[lim_idx], axisLabel_e, 'e', e_fit, weights[lim_idx], BH_fixed, norange)
    if e_lim:
        pylab.xlim(e_lim)
    del e


    # Inclination
    if not BH_fixed:
        pylab.subplot(5,3,13)
    else:
        pylab.subplot(2,3,4)    

        if not allparamsinput:
            print 'Reading i file...'
            file_i = open(mnestfile[0:-4]+'_i.txt')
            if 'new' in mnestfile or new:
                i = np.array(file_i.readlines(),dtype=float)
            else:
                i = np.array(file_i.readline().split(),dtype=float)
            file_i.close()

    print '\nplotting i'
    i *= (180./math.pi)
    if not i_lim:
        i_lim = [i.min(),i.max()]
    lim_idx = np.where((i>=i_lim[0]) & (i<=i_lim[1]))
    plotProb(i[lim_idx], axisLabel_i, 'i', i_fit, weights[lim_idx], BH_fixed, norange)
    if i_lim:
        pylab.xlim(i_lim)
    del i

    # omega
    if not BH_fixed:
        pylab.subplot(5,3,14)
    else:
        pylab.subplot(2,3,5)

        if not allparamsinput:
            print 'Reading omega file...'
            file_omega = open(mnestfile[0:-4]+'_omega.txt')
            if 'new' in mnestfile or new:
                omega = np.array(file_omega.readlines(),dtype=float)
            else:
                omega = np.array(file_omega.readline().split(),dtype=float)
            file_omega.close()

    print '\nplotting omega'
    omega *= (180./math.pi)
    if not o_lim:
        o_lim = [omega.min(),omega.max()]
    lim_idx = np.where((omega>=o_lim[0]) & (omega<=o_lim[1]))
    plotProb(omega[lim_idx], axisLabel_o, 'w', omega_fit, weights[lim_idx], BH_fixed, norange)
    if o_lim:
        pylab.xlim(o_lim)
    del omega


    # Omega
    if not BH_fixed:
        pylab.subplot(5,3,15)
    else:
        pylab.subplot(2,3,6)        

        if not allparamsinput:
            print 'Reading Omega file...'
            file_Omega = open(mnestfile[0:-4]+'_Omega.txt')
            if 'new' in mnestfile or new:
                Omega = np.array(file_Omega.readlines(),dtype=float)
            else:
                Omega = np.array(file_Omega.readline().split(),dtype=float)
            file_Omega.close()

    print '\nplotting Omega'
    Omega *= (180./math.pi)
    if not O_lim:
        O_lim = [Omega.min(),Omega.max()]
    lim_idx = np.where((Omega>=O_lim[0]) & (Omega<=O_lim[1]))
    plotProb(Omega[lim_idx], axisLabel_O, 'O', Omega_fit, weights[lim_idx], BH_fixed, norange)
    if O_lim:
        pylab.xlim(O_lim)
    del Omega

    pylab.rc('axes', titlesize=14, labelsize=14)
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)



def plot_SgrA_astrometric_contours(mnestfile='',color='black',linestyle='solid',star_num=0,bin=40,
                                   Vx_lim=None,Vy_lim=None,xo_lim=None,yo_lim=None,linewidth=None,
                                   inFile=None):

    if not inFile:
        inFile = asciidata.open(mnestfile)

    weights = inFile[0].tonumpy()
    logLike = inFile[1].tonumpy()

    xo = -inFile[3].tonumpy()
    xo *= 1000
    yo = inFile[4].tonumpy()
    yo *= 1000
    Vx = -inFile[5].tonumpy()
    Vx *= 1000
    Vy = inFile[6].tonumpy()
    Vy *= 1000

    axisLabel_x0 = 'Sgr A* $\Delta$RA Position (mas)'
    axisLabel_y0 = 'Sgr A* Y Position (mas)'
    axisLabel_Vx = 'Sgr A* $\Delta$RA Velocity (mas/yr)'
    axisLabel_Vy = 'Sgr A* Y Velocity (mas/yr)'

    pylab.figure(1,figsize=(10,5))
    pylab.subplot(121)
    (hist, yobins, xobins) = pylab.histogram2d(yo, xo, bins=(bin,bin),weights=weights)
    
    # Need to convert the 2d histogram into floats
    probDist = np.array(hist, dtype=float)
    levels = getContourLevels(probDist)

    pylab.contour(probDist, levels, origin=None, colors=color,
                  extent=[xobins[0], xobins[-1], yobins[0], yobins[-1]],
                  linestyles=linestyle,linewidths=linewidth)
    #pylab.axis('equal')
    pylab.xlabel(axisLabel_x0)
    pylab.ylabel(axisLabel_y0)
    if xo_lim:
        pylab.xlim(xo_lim)
    if yo_lim:
        pylab.ylim(yo_lim)


    pylab.subplot(122)
    (hist, Vybins, Vxbins) = pylab.histogram2d(Vy, Vx, bins=(bin,bin),weights=weights)
    
    # Need to convert the 2d histogram into floats
    probDist = np.array(hist, dtype=float)
    levels = getContourLevels(probDist)

    print linestyle
    pylab.contour(probDist, levels, origin=None, colors=color,
                  extent=[Vxbins[0], Vxbins[-1], Vybins[0], Vybins[-1]],
                  linestyles=linestyle,linewidths=linewidth)
    #pylab.axis('equal')
    pylab.xlabel(axisLabel_Vx)
    pylab.ylabel(axisLabel_Vy)
    if Vx_lim:
        pylab.xlim(Vx_lim)
    if Vy_lim:
        pylab.ylim(Vy_lim)


    pylab.subplots_adjust(bottom=0.15)#bottom=0.05, left=0.09, right=0.95, top=0.97)




def plot_SgrA_astrometric_env(summaryfile,color='black',star_num=0,row_num=0,label=''):

    inFile = asciidata.open(summaryfile)

    n_params = (inFile.ncols-2)/4

    x_idx,y_idx,vx_idx,vy_idx = 1,2,3,4
    xo, yo, vx, vy = -inFile[x_idx][row_num],inFile[y_idx][row_num],-inFile[vx_idx][row_num],inFile[vy_idx][row_num]
    dxo, dyo, dvx, dvy = inFile[x_idx+n_params][row_num],inFile[y_idx+n_params][row_num],inFile[vx_idx+n_params][row_num],inFile[vy_idx+n_params][row_num]

    to = 2000.0
    t1 = 1995.0
    t2 = 2015.0
    dt = 0.1

    x_max = []
    x_min = []
    y_max = []
    y_min = []
    times = np.arange(t1,t2,dt)
    for t in times:
        x_max.append(np.max([(xo+dxo) + (to-t)*(vx+dvx),(xo+dxo) + (to-t)*(vx-dvx),(xo-dxo) + (to-t)*(vx+dvx),(xo-dxo) + (to-t)*(vx-dvx)]))
        y_max.append(np.max([(yo+dyo) + (to-t)*(vy+dvy),(yo+dyo) + (to-t)*(vy-dvy),(yo-dyo) + (to-t)*(vy+dvy),(yo-dyo) + (to-t)*(vy-dvy)]))
        x_min.append(np.min([(xo+dxo) + (to-t)*(vx+dvx),(xo+dxo) + (to-t)*(vx-dvx),(xo-dxo) + (to-t)*(vx+dvx),(xo-dxo) + (to-t)*(vx-dvx)]))
        y_min.append(np.min([(yo+dyo) + (to-t)*(vy+dvy),(yo+dyo) + (to-t)*(vy-dvy),(yo-dyo) + (to-t)*(vy+dvy),(yo-dyo) + (to-t)*(vy-dvy)]))
        

    axisLabel_x0 = 'Sgr A* $\Delta$RA Position (mas)'
    axisLabel_y0 = 'Sgr A* Y Position (mas)'


    pylab.subplot(211)
    pylab.plot([t1,t2],[xo+(to-t1)*vx,xo+(to-t2)*vx],color=color,linestyle='solid',label=label)
    pylab.plot(times,x_max,linestyle='dotted',color=color)
    pylab.plot(times,x_min,linestyle='dotted',color=color)
    pylab.ylabel(axisLabel_x0)

    if label:
        pylab.legend()

    pylab.subplot(212)
    pylab.plot([t1,t2],[yo+(to-t1)*vy,yo+(to-t2)*vy],color=color,linestyle='solid')
    pylab.plot(times,y_max,linestyle='dotted',color=color)
    pylab.plot(times,y_min,linestyle='dotted',color=color)
    pylab.xlabel('Time (years)')
    pylab.ylabel(axisLabel_y0)


def get_xy_pred_orig(mnestfile='/u/syelda/research/gc/aligndir/test_schoedel/13_01_16/efit/chains/efit_S0-2_S0-38.txt',outDir='',star_num=1):
    '''
    star_num: 0 for the 1st star, 1 for the 2nd star, and so on
    '''

    # READ IN LINES
    f = open(mnestfile)
    lines = f.readlines()
    
    # INIT ARRAYS
    n_stars = 2
    n_pars = 7 + 6*n_stars

    arr = np.zeros((len(lines),n_pars),dtype='d')
    #arr = np.zeros((n_pars,len(lines)),dtype='d')
    weights = np.zeros(len(lines))
    loglikes = np.zeros(len(lines))

    # READ LINES INTO ARRAYS
    for i in range(len(lines)):
        line = np.array(lines[i].split(),dtype='d')
        weights[i] = line[0]
        loglikes[i] = line[1]
        for j in range(n_pars):
            arr[i][j] = line[j+2]
            #arr[j][i] = line[j+2]
        
    print arr.shape
    print weights.shape

    epochs_table = asciidata.open('epochs.txt')
    epochs = epochs_table[0].tonumpy()

    # Loop over epochs
    for e in range(len(epochs)):
        t = epochs[e]
        print 'Epoch '+str(e)+': '+str(t)
        xy_table = asciidata.create(6,len(lines))#np.zeros((len(lines),3))

        # Loop over lines in the multinest output
        for i in range(len(lines)):

            xy_table[0][i] = t
            xy_table[5][i] = weights[i]

            elem = np.zeros(8)    
            elem[0] = arr[i][6]  # Distance        
            elem[2] = arr[i][10+6*star_num] # Period
            elem[3] = arr[i][12+6*star_num] # Eccentricity
            elem[4] = arr[i][11+6*star_num] #t0
            elem[5] = arr[i][8+6*star_num]  #w
            elem[6] = arr[i][9+6*star_num]  #i
            elem[7] = arr[i][7+6*star_num]  #Omega

            # Get a from M and Period
            mass = arr[i][0]
            #a_AU = (mass*elem[2]**2.)**(1./3.)   # a in AU = (Mass*Period^2)^(1/3)
            #a = a_AU/elem[0]  # a in arcsec
            #elem[1] = a
        
            #xy_table[1][i],xy_table[2][i],xy_table[3][i] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

            # Add drift
            xo = arr[i][1]  # in arcsec
            yo = arr[i][2]
            Vxo = arr[i][3]  # in arcsec/yr
            Vyo = arr[i][4]
            Vzo = arr[i][5]

            drift_params = (xo, yo, Vxo, Vyo, Vzo)
        
            x, y, z, vx, vy, vz, v = get_orbit_prediction(elem, t, mass, drift_params)

            xy_table[1][i] = x
            xy_table[2][i] = y
            xy_table[3][i] = z
            xy_table[4][i] = vz

            #xy_table[1][i] = -xy_table[1][i] + xo+ Vxo * (t-2000.)	# arcsec in camera coordinates 
            #xy_table[2][i] += yo + Vyo  * (t-2000.)		# arcsec in camera coordinates 
            #xy_table[3][i] *= elem[0]			        # convert to AU

            #vz *=  elem[0] * 15./math.pi     # convert to km/s
            #vz += Vzo

            #xy_table[4][i] = vz

        print 'Writing out to: '+outDir+'model_epoch'+str(e)+'.out'
        xy_table.writeto(outDir+'model_epoch'+str(e)+'.out')
    


def read_next_value(file):
    '''
    Function to read the next parameter value in the one-line efit_[STAR]_combo_[PARAM].txt file outputted by master_wrapper.
    '''
    
    # Check if next char is '-'.  If so, the next value is negative and therefore takes up 25 bytes total instead of 24
    char1 = file.read(1)    
    if char1 == '-':
        nbytes_rest = 24
    else:
        nbytes_rest = 23

    rest = file.read(nbytes_rest)
    junk = file.read(1)
    
    return float(char1 + rest)
            


def get_xy_pred_BHfixed(mnestfile='/u/syelda/research/gc/aligndir/test_schoedel/13_01_16/efit/chains/efit_S0-2_S0-38.txt',
                        outDir='',epochs=None,epochs_idx=None):

    # READ IN FILES
    print 'Reading in files'
    file_P = open(mnestfile[0:-4]+'_P.txt')
    file_i = open(mnestfile[0:-4]+'_i.txt')
    file_e = open(mnestfile[0:-4]+'_e.txt')
    file_Omega = open(mnestfile[0:-4]+'_Omega.txt')
    file_omega = open(mnestfile[0:-4]+'_omega.txt')
    file_To = open(mnestfile[0:-4]+'_To.txt')

    file_xo = open(mnestfile[0:-4]+'_xo.txt')
    file_yo = open(mnestfile[0:-4]+'_yo.txt')
    file_Vx = open(mnestfile[0:-4]+'_Vx.txt')
    file_Vy = open(mnestfile[0:-4]+'_Vy.txt')
    file_Vz = open(mnestfile[0:-4]+'_Vz.txt')
    file_D = open(mnestfile[0:-4]+'_D.txt')
    file_mass = open(mnestfile[0:-4]+'_mass.txt')


    # GET ARRAY OF WEIGHTS
    print 'Getting array of weights'
    file_weights = open(mnestfile[0:-4]+'_weights.txt')
    weights = np.array(file_weights.readlines(),dtype=float)
    file_weights.close()

    #if not epochs:
    #    epochs_table = asciidata.open('epochs.txt')
    #    epochs = epochs_table[0].tonumpy()

    #print 'Looping over epochs'
    # Loop over epochs
    #if epochs_range:
    #    epochs_idx = range(epochs_range[0],epochs_range[1])
    #else:
    #    epochs_idx = range(len(epochs))


    for e in epochs_idx:
        t = epochs[e]
        print 'Epoch '+str(e)+': '+str(t)
        #xy_table = asciidata.create(6,len(weights))

        # go to beginning of the files
        file_P.seek(0)
        file_e.seek(0)
        file_To.seek(0)
        file_omega.seek(0)
        file_Omega.seek(0)
        file_i.seek(0)

        file_mass.seek(0)
        file_D.seek(0)
        file_xo.seek(0)
        file_yo.seek(0)
        file_Vx.seek(0)
        file_Vy.seek(0)
        file_Vz.seek(0)

        xy_outfile = open(outDir+'model_epoch'+str(e)+'.out','w')
        # Loop over lines in the multinest output
        nweights = len(weights)
        for i in range(nweights):

            #print i/(float(nweights))
        
            ##xy_table[0][i] = t
            ##xy_table[5][i] = weights[i]

            elem = np.zeros(8)    

        # Distance        
            elem[0] = float(file_D.readline())
        # Period
            elem[2] = float(file_P.readline())
        # Eccentricity
            elem[3] = float(file_e.readline())
        # t0
            elem[4] = float(file_To.readline())
        # w
            elem[5] = float(file_omega.readline())
        # i
            elem[6] = float(file_i.readline())
        # Omega
            elem[7] = float(file_Omega.readline())

        # Get a from M and Period
            mass = float(file_mass.readline())
        
            #xy_table[1][i],xy_table[2][i],xy_table[3][i] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

        # Add drift
            xo =  float(file_xo.readline()) # in arcsec
            yo = float(file_yo.readline())
            Vxo = float(file_Vx.readline())  # in arcsec/yr
            Vyo = float(file_Vy.readline())
            Vzo = float(file_Vz.readline())

            drift_params = (xo, yo, Vxo, Vyo, Vzo)
        
            x, y, z, vx, vy, vz, v = get_orbit_prediction(elem, t, mass, drift_params)

            line = '%8s\t%16s\t%16s\t%16s\t%16s\t%16s\n' % (t, x, y, z, vz, weights[i])
            xy_outfile.write(line)


        print 'Writing out to: '+outDir+'model_epoch'+str(e)+'.out'
        xy_outfile.close()
    


def get_xy_pred_BHfixed_old(mnestfile='/u/syelda/research/gc/aligndir/test_schoedel/13_01_16/efit/chains/efit_S0-2_S0-38.txt',
                        outDir='',epochs=None,epochs_idx=None,new=False):

    # READ IN FILES
    print 'Reading in files'
    file_P = open(mnestfile[0:-4]+'_P.txt')
    file_i = open(mnestfile[0:-4]+'_i.txt')
    file_e = open(mnestfile[0:-4]+'_e.txt')
    file_Omega = open(mnestfile[0:-4]+'_Omega.txt')
    file_omega = open(mnestfile[0:-4]+'_omega.txt')
    file_To = open(mnestfile[0:-4]+'_To.txt')

    file_xo = open(mnestfile[0:-4]+'_xo.txt')
    file_yo = open(mnestfile[0:-4]+'_yo.txt')
    file_Vx = open(mnestfile[0:-4]+'_Vx.txt')
    file_Vy = open(mnestfile[0:-4]+'_Vy.txt')
    file_Vz = open(mnestfile[0:-4]+'_Vz.txt')
    file_D = open(mnestfile[0:-4]+'_D.txt')
    file_mass = open(mnestfile[0:-4]+'_mass.txt')


    # GET ARRAY OF WEIGHTS
    print 'Getting array of weights'
    file_weights = open(mnestfile[0:-4]+'_weights.txt')
    weights = np.array(file_weights.readline().split(),dtype=float)
    file_weights.close()

    #if not epochs:
    #    epochs_table = asciidata.open('epochs.txt')
    #    epochs = epochs_table[0].tonumpy()

    #print 'Looping over epochs'
    # Loop over epochs
    #if epochs_range:
    #    epochs_idx = range(epochs_range[0],epochs_range[1])
    #else:
    #    epochs_idx = range(len(epochs))


    for e in epochs_idx:
        t = epochs[e]
        print 'Epoch '+str(e)+': '+str(t)
        #xy_table = asciidata.create(6,len(weights))

        # go to beginning of the files
        file_P.seek(0)
        file_e.seek(0)
        file_To.seek(0)
        file_omega.seek(0)
        file_Omega.seek(0)
        file_i.seek(0)

        file_mass.seek(0)
        file_D.seek(0)
        file_xo.seek(0)
        file_yo.seek(0)
        file_Vx.seek(0)
        file_Vy.seek(0)
        file_Vz.seek(0)

        xy_outfile = open(outDir+'model_epoch'+str(e)+'.out','w')
        # Loop over lines in the multinest output
        nweights = len(weights)
        for i in range(nweights):

            #print i/(float(nweights))
        
            ##xy_table[0][i] = t
            ##xy_table[5][i] = weights[i]

            elem = np.zeros(8)    

        # Distance        
            elem[0] = read_next_value(file_D)
        # Period
            elem[2] = read_next_value(file_P)
        # Eccentricity
            elem[3] = read_next_value(file_e)
        # t0
            elem[4] = read_next_value(file_To)
        # w
            elem[5] = read_next_value(file_omega)
        # i
            elem[6] = read_next_value(file_i)
        # Omega
            elem[7] = read_next_value(file_Omega)

        # Get a from M and Period
            mass = read_next_value(file_mass)
        
            #xy_table[1][i],xy_table[2][i],xy_table[3][i] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

        # Add drift
            xo =  read_next_value(file_xo) # in arcsec
            yo = read_next_value(file_yo)
            Vxo = read_next_value(file_Vx)  # in arcsec/yr
            Vyo = read_next_value(file_Vy)
            Vzo = read_next_value(file_Vz)

            drift_params = (xo, yo, Vxo, Vyo, Vzo)
        
            x, y, z, vx, vy, vz, v = get_orbit_prediction(elem, t, mass, drift_params)

            line = '%8s\t%16s\t%16s\t%16s\t%16s\t%16s\n' % (t, x, y, z, vz, weights[i])
            xy_outfile.write(line)


        print 'Writing out to: '+outDir+'model_epoch'+str(e)+'.out'
        xy_outfile.close()
    


def get_xy_pred_BHfixed_onecore(mnestfile,linerange,t,input_p,corenum):

    min,max = linerange
    nlines = max - min

    out_arr = np.zeros((nlines,2),dtype=float)

    # READ IN FILES
    print 'CORE '+str(corenum)+': Reading in files'
    file_P = open(mnestfile[0:-4]+'_P.txt')
    file_i = open(mnestfile[0:-4]+'_i.txt')
    file_e = open(mnestfile[0:-4]+'_e.txt')
    file_Omega = open(mnestfile[0:-4]+'_Omega.txt')
    file_omega = open(mnestfile[0:-4]+'_omega.txt')
    file_To = open(mnestfile[0:-4]+'_To.txt')

    file_xo = open(mnestfile[0:-4]+'_xo.txt')
    file_yo = open(mnestfile[0:-4]+'_yo.txt')
    file_Vx = open(mnestfile[0:-4]+'_Vx.txt')
    file_Vy = open(mnestfile[0:-4]+'_Vy.txt')
    file_Vz = open(mnestfile[0:-4]+'_Vz.txt')
    file_D = open(mnestfile[0:-4]+'_D.txt')
    file_mass = open(mnestfile[0:-4]+'_mass.txt')

    files = [file_P,file_i,file_e,file_Omega,file_omega,file_To,file_xo,file_yo,file_Vx,file_Vy,file_Vz,file_D,file_mass]

    print 'CORE '+str(corenum)+': Skipping lines...'
    for i in range(min):
        for file in files:
            file.next()

    print 'CORE '+str(corenum)+': Reading lines...'
    row = 0 # array row
    print min,max
    for i in range(min,max):
        print i

        elem = np.zeros(8)    

        # Distance        
        elem[0] = float(file_D.next())
        # Period
        elem[2] = float(file_P.next())
        # Eccentricity
        elem[3] = float(file_e.next())
        # t0
        elem[4] = float(file_To.next())
        # w
        elem[5] = float(file_omega.next())
        # i
        elem[6] = float(file_i.next())
        # Omega
        elem[7] = float(file_Omega.next())

        # Get a from M and Period
        mass = float(file_mass.next())
        
            #xy_table[1][i],xy_table[2][i],xy_table[3][i] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

        # Add drift
        xo =  float(file_xo.next()) # in arcsec
        yo = float(file_yo.next())
        Vxo = float(file_Vx.next())  # in arcsec/yr
        Vyo = float(file_Vy.next())
        Vzo = float(file_Vz.next())

        drift_params = (xo, yo, Vxo, Vyo, Vzo)
        
        x, y, z, vx, vy, vz, v = get_orbit_prediction(elem, t, mass, drift_params)

        out_arr[row][0] = x
        out_arr[row][1] = y

    print 'CORE '+str(corenum)+': Sending arr through pipe...'
    input_p.send(out_arr)


def get_xy_pred_BHfixed_multicore(mnestfile='/u/syelda/research/gc/aligndir/test_schoedel/13_01_16/efit/chains/efit_S0-2_S0-38.txt',outDir='',epochs=None,epochs_idx=None,ncores=1):

    start_time = time.strftime('%X %x %Z')
  
    # GET ARRAY OF WEIGHTS
    print 'Getting array of weights'
    file_weights = open(mnestfile[0:-4]+'_weights.txt')
    weights = np.array(file_weights.readlines(),dtype=float)
    file_weights.close()
    nweights = len(weights)

    #if not epochs:
    #    epochs_table = asciidata.open('epochs.txt')
    #   epochs = epochs_table[0].tonumpy()
    

    # get the line ranges for each core
    div = nweights/ncores
    mod = nweights % ncores

    lineranges = []
    max = 0
    for n in range(0,ncores):
        min = max
        if mod != 0:
            max = min + div + 1
            mod = mod - 1
        else:
            max = min + div

        lineranges.append((min,max))
        print min,max

    print lineranges

    print 'Looping over epochs'
    for e in epochs_idx:
        t = epochs[e]
        print 'Epoch '+str(e)+': '+str(t)

        processes = []
        outs = []
        ins = []

        div = nweights/ncores
        mod = nweights % ncores

        print div, mod

        max = 0     # this range calc doesn't really need to be done separately for each epoch...
        #output_dict = {}
        #for n in range(0,ncores):
        #    output_dict[str(n)] = 0.0
        for n in range(0,ncores):

            output_p, input_p = multiprocessing.Pipe()
            outs.append(output_p)
            ins.append(input_p)

            p = multiprocessing.Process(target=get_xy_pred_BHfixed_onecore,args=(mnestfile,lineranges[n],t,input_p,n))
            if n == 0:
                processes.append(p)
                p.start()   

        print 'Getting output'
        for n in range(ncores):
            output_dict[str(n)] = outs[n].recv()


        print 'Waiting for processes...'
        for p in processes:
            p.join()


        print 'Writing out file...'
        xy_outfile = open(outDir+'model_epoch'+str(e)+'.out','w')

        w = 0
        for k in output_dict.keys():
            arr = output_dict[k]
            siz = arr.shape
            for row in range(siz[0]):
                line = '%8s\t%16s\t%16s\t%16s\t%16s\t%16s\n' % (0.0, arr[row][0], arr[row][1], 
                                                                0.0, 0.0, weights[w])                
                xy_outfile.write(line)
                w += 0

        xy_outfile.close()

    end_time = time.strftime('%X %x %Z')

    print 'START -> STOP TIME:',start_time,end_time


def make_efit_outputfile(paramslimit_file,paramslimit_BHfile,star='S0-38',outputfile='orbit.S0-38.output'):

    table = asciidata.open(paramslimit_file)
    peakprob = table[2].tonumpy()
    
    P = peakprob[0]
    To = peakprob[1]
    ecc = peakprob[2]
    incl = peakprob[3]
    om = peakprob[4]
    OMEGA = peakprob[5]

    tableBH = asciidata.open(paramslimit_BHfile)
    #bestfitBH = tableBH[1].tonumpy()
    peakprobBH = tableBH[2].tonumpy()
    
    Xo = -peakprobBH[0]*1e-3
    Yo = peakprobBH[1]*1e-3
    Ro = peakprobBH[2]*1e3
    Vx = -peakprobBH[3]*1e-3
    Vy = peakprobBH[4]*1e-3
    Vz = peakprobBH[5]
    mass = peakprobBH[6]*1e6

    chi = 0.0

    efit_file = open(outputfile, "w")
    efit_file.write( "1 1 1 1 " + str(chi) + " 1 " + str(mass) + " " + str(Xo) + " " + str(Yo) + " " + str(Vx) + " " + 
                     str(Vy) + " " + str(Ro) + " " + str(Vz) + " "+star+" " + str(OMEGA) + " " + str(om) + " " + 
                     str(incl) + " " + str(ecc) + " " + str(P) + " " + str(To) + " 0 \n")    
    
    efit_file.close()

def make_model(summary_file,star,star_num=0,dt=0.050,t1=1995.,t2=2014.):
    ''' Make model file from the efit5 output '''

    table = asciidata.open(summary_file)
    nmodes = table.nrows
    nparams = table.ncols/4

    print 'Number of Modes Found:',nmodes
    print 'Number of Parameters in Fit:',nparams
    
    times = np.arange(t1,t2,dt)
    out_table = asciidata.create(7,len(times))

    for m in range(nmodes):
        for j in range(len(times)):
            t = times[j]
            out_table[0][j] = t
        
            # MAP parameters
            elem = np.zeros(8)
            elem[0] = table[6+nparams*3][m]  # Distance        
            elem[2] = table[10+6*star_num+nparams*3][m] # Period
            elem[3] = table[12+6*star_num+nparams*3][m] # Eccentricity
            elem[4] = table[11+6*star_num+nparams*3][m] #t0
            elem[5] = table[8+6*star_num+nparams*3][m]  #w
            elem[6] = table[9+6*star_num+nparams*3][m]  #i
            elem[7] = table[7+6*star_num+nparams*3][m]  #Omega

            if j == 0:
                print 'Best fit period:',elem[2]
            
            # Get a from M and Period
            mass = table[0+nparams*3][m]

            #a_AU = (mass*elem[2]**2.)**(1./3.)   # a in AU = (Mass*Period^2)^(1/3)
            #a = a_AU/elem[0]  # a in arcsec
            #elem[1] = a

            #out_table[1][j],out_table[2][j],out_table[3][j] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

            # Add drift
            xo = table[1+nparams*3][m]  # in arcsec
            yo = table[2+nparams*3][m]
            Vxo = table[3+nparams*3][m]  # in arcsec/yr
            Vyo = table[4+nparams*3][m]
            Vzo = table[5+nparams*3][m]
        
            drift_params = (xo, yo, Vxo, Vyo, Vzo)

            #out_table[1][j] = -out_table[1][j] + xo+ Vxo * (t-2000.);	# arcsec in camera coordinates 
            #out_table[2][j] += yo + Vyo  * (t-2000.);		# arcsec in camera coordinates 
            #out_table[3][j] *= elem[0];			        # convert to AU

            #vz *=  elem[0] * 15./math.pi;     # convert to km/s
            #vz += Vzo;

            #out_table[4][j],out_table[5][j],out_table[6][j] = vx, vy, vz

            x, y, z, vx, vy, vz, v = get_orbit_prediction(elem, t, mass, drift_params)

            out_table[1][j] = x
            out_table[2][j] = y
            out_table[3][j] = z
            out_table[4][j] = vx
            out_table[5][j] = vy
            out_table[6][j] = vz

        out_table[0].reformat('%4.3f')
        out_table.writeto('orbit.'+star+'.'+str(m+1)+'.model')
    

def make_model_lim(modelfile,rangefile,star,sigma=1,make_envelope=False):

    model = asciidata.open(modelfile)
    orbit_range = asciidata.open(rangefile)

    # flip the sign of Vz in the model file
    for i in range(model.nrows):
        model[6][i] = -model[6][i]

    model_times = model[0].tonumpy()
    range_times = orbit_range[0].tonumpy()

    model_lim_xy = asciidata.create(5,orbit_range.nrows)
    model_lim_vel = asciidata.create(7, orbit_range.nrows)
    if make_envelope:
        envelope_plus = asciidata.create(2, orbit_range.nrows)
        envelope_minus = asciidata.create(2, orbit_range.nrows)

    for range_idx in range(len(range_times)):
        time = range_times[range_idx]
        model_idx = np.where(model_times == time)[0][0]
        
        model_lim_xy[0][range_idx] = time
        model_lim_vel[0][range_idx] = time

        x = model[1][model_idx]
        y = model[2][model_idx]
        vz = model[6][model_idx]

        model_lim_xy[1][range_idx] = x - sigma*orbit_range[1][range_idx] #xmin
        model_lim_xy[2][range_idx] = x + sigma*orbit_range[1][range_idx] #xmax
        model_lim_xy[3][range_idx] = y - sigma*orbit_range[2][range_idx] #ymin
        model_lim_xy[4][range_idx] = y + sigma*orbit_range[2][range_idx] #ymax
    
        model_lim_vel[5][range_idx] = vz - sigma*orbit_range[3][range_idx] #vzmin
        model_lim_vel[6][range_idx] = vz + sigma*orbit_range[3][range_idx] #vzmax

        if orbit_range.ncols > 4 and make_envelope:
            r = math.sqrt(x**2 + y**2)
            phi = math.atan(y/x)  # radians

            r_plus = r + orbit_range[4][range_idx]
            r_minus = r - orbit_range[4][range_idx]
            
            x_plus = r_plus*math.cos(phi)
            y_plus = r_plus*math.sin(phi)

            x_minus = r_minus*math.cos(phi)
            y_minus = r_minus*math.sin(phi)
            
            if x < 0:
                x_minus = -x_minus
                y_minus = -y_minus
                x_plus = -x_plus
                y_plus = -y_plus

            envelope_plus[0][range_idx],envelope_plus[1][range_idx] = x_plus, y_plus
            envelope_minus[0][range_idx],envelope_minus[1][range_idx] = x_minus, y_minus

    model_lim_xy.writeto('orbit.'+star+'.model_lim_xy')
    model_lim_vel.writeto('orbit.'+star+'.model_lim_vel')

    if make_envelope:
        envelope_plus.writeto('envelope_plus_simple_'+star+'.txt')
        envelope_minus.writeto('envelope_minus_simple_'+star+'.txt')

    model.writeto(modelfile)


def weighted_stdev(values, weights):

    avg_weighted = np.average(values,weights=weights)
    var_weighted = np.sum(weights*(values - avg_weighted)**2)

    return var_weighted**0.5
    

def get_orbit_prediction(elem, t, mass, drift_params):
          
    # Get a from M and Period
    a_AU = (mass*elem[2]**2.)**(1./3.)   # a in AU = (Mass*Period^2)^(1/3)
    a = a_AU/elem[0]  # a in arcsec
    elem[1] = a

    x,y,z = orbit_position(elem, t)
    vx, vy, vz, v = orbit_vel(elem, t)

    xo,yo,Vxo,Vyo,Vzo = drift_params

    x = -x + xo + Vxo*(t-2000.)   # arcsec in camera coordinates 
    y += yo + Vyo*(t-2000.)       # arcsec in camera coordinates 
    z *= elem[0]                  # convert to AU
    
    vz *=  elem[0] * 15./math.pi  # convert to km/s
    vz += Vzo

    return x, y, z, vx, vy, vz, v


def make_model_orbitparams(orbit_params,star,tmin=1995.0,tmax=2013.0):
    ''' Make model file from list of orbital parameters '''

    #table = asciidata.open(summary_file)
    #nmodes = table.nrows
    #nparams = table.ncols/4

    #print 'Number of Modes Found:',nmodes
    #print 'Number of Parameters in Fit:',nparams
    
    times = np.arange(tmin,tmax,0.050)
    out_table = asciidata.create(7,len(times))

    #for m in range(nmodes):
    for j in range(len(times)):
        t = times[j]
        out_table[0][j] = t
        
        # MAP parameters
        elem = np.zeros(8)
        elem[0] = orbit_params[6]  # Distance        
        elem[2] = orbit_params[10]  # Period
        elem[3] = orbit_params[12] # Eccentricity
        elem[4] = orbit_params[11] #t0
        elem[5] = orbit_params[8]  #w
        elem[6] = orbit_params[9]  #i
        elem[7] = orbit_params[7]  #Omega

        # Get a from M and Period
        mass = orbit_params[0]

            #a_AU = (mass*elem[2]**2.)**(1./3.)   # a in AU = (Mass*Period^2)^(1/3)
            #a = a_AU/elem[0]  # a in arcsec
            #elem[1] = a

            #out_table[1][j],out_table[2][j],out_table[3][j] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

        # Add drift
        xo = orbit_params[1]
        yo = orbit_params[2]
        Vxo = orbit_params[3]
        Vyo = orbit_params[4]
        Vzo = orbit_params[5]
        
        drift_params = (xo, yo, Vxo, Vyo, Vzo)

            #out_table[1][j] = -out_table[1][j] + xo+ Vxo * (t-2000.);	# arcsec in camera coordinates 
            #out_table[2][j] += yo + Vyo  * (t-2000.);		# arcsec in camera coordinates 
            #out_table[3][j] *= elem[0];			        # convert to AU

            #vz *=  elem[0] * 15./math.pi;     # convert to km/s
            #vz += Vzo;

            #out_table[4][j],out_table[5][j],out_table[6][j] = vx, vy, vz

        x, y, z, vx, vy, vz, v = get_orbit_prediction(elem, t, mass, drift_params)

        out_table[1][j] = x
        out_table[2][j] = y
        out_table[3][j] = z
        out_table[4][j] = vx
        out_table[5][j] = vy
        out_table[6][j] = vz

    out_table.writeto('orbit.'+star+'.model')
    

def get_points(modelfile,epochs):

    model = asciidata.open(modelfile)
    times = model[0].tonumpy()

    for i in range(len(epochs)):
        epoch = epochs[i]
        idx = np.argmin(np.abs(epoch-times))
        
        
    return


#####################
# TEST/OLD FUNCTIONS
#####################

def test_r_lim():
    ''' test the code to get the envelope plot!  ultimately used in make_model_lim'''

    def get_envelope(x,y,delR):
        r = math.sqrt(x**2 + y**2)
        phi = math.atan(y/x)

        r_plus = r + delR
        r_minus = r - delR

        x_plus = r_plus*math.cos(phi)
        y_plus = r_plus*math.sin(phi)
       
        x_minus = r_minus*math.cos(phi)
        y_minus = r_minus*math.sin(phi)
            
        if x < 0:
            x_minus = -x_minus
            y_minus = -y_minus
            x_plus = -x_plus
            y_plus = -y_plus

        return ([x_minus, x, x_plus], [y_minus, y, y_plus])

    # pos x, pos y 
    x = 2.
    y = 2.
    delR = 1.

    x_vals,y_vals = get_envelope(x,y,delR)
    pylab.plot([0]+x_vals,[0]+y_vals,marker='.',color='black')

    # neg x, pos y
    x = -2.
    y = 2.
    delR = 1.

    x_vals,y_vals = get_envelope(x,y,delR)
    pylab.plot([0]+x_vals,[0]+y_vals,marker='.',color='red')

    # neg x, neg y
    x = -2.
    y = -2.
    delR = 1.

    x_vals,y_vals = get_envelope(x,y,delR)
    pylab.plot([0]+x_vals,[0]+y_vals,marker='.',color='blue')

    # pos x, neg y
    x = 2.
    y = -2.
    delR = 1.

    x_vals,y_vals = get_envelope(x,y,delR)
    pylab.plot([0]+x_vals,[0]+y_vals,marker='.',color='green')

    pylab.axis('equal')
    pylab.show()

def make_model_lim_old(mnestfile,summary_file,star,star_num=0,n_stars=1):

    outDir = ''

    # READ IN LINES
    f = open(mnestfile)
    lines = f.readlines()
    
    # INIT ARRAYS
    #n_stars = 2
    nparams = 7 + 6*n_stars

    summary = asciidata.open(summary_file)
    m = np.argmax(summary[summary.ncols-1].tonumpy())

    arr = np.zeros((len(lines),nparams),dtype='d')
    weights = np.zeros(len(lines))
    loglikes = np.zeros(len(lines))

    # READ LINES INTO ARRAYS
    for i in range(len(lines)):
        line = np.array(lines[i].split(),dtype='d')
        weights[i] = line[0]
        loglikes[i] = line[1]
        for j in range(nparams):
            arr[i][j] = line[j+2]

    #times = np.arange(1995.,2013.,0.050)
    times = np.arange(2013.5,2015.1,0.1)
    out_table = asciidata.create(7,len(times))
    out_table_vel = asciidata.create(7,len(times))
    #print out_table
    
    for j in range(len(times)):
        t = times[j]
        out_table[0][j] = t
        out_table_vel[0][j] = t
        print t

        # BEST FIT
        n_param = summary.ncols/4
        elem = np.zeros(8)
        elem[0] = summary[6+n_param*3][m]  # Distance        
        elem[2] = summary[10+6*star_num+n_param*3][m] # Period
        elem[3] = summary[12+6*star_num+n_param*3][m] # Eccentricity
        elem[4] = summary[11+6*star_num+n_param*3][m] #t0
        elem[5] = summary[8+6*star_num+n_param*3][m]  #w
        elem[6] = summary[9+6*star_num+n_param*3][m]  #i
        elem[7] = summary[7+6*star_num+n_param*3][m]  #Omega

        # Get a from M and Period
        mass = summary[0+n_param*3][m]

        # Add drift
        xo = summary[1+n_param*3][m]  # in arcsec
        yo = summary[2+n_param*3][m]
        Vxo = summary[3+n_param*3][m]  # in arcsec/yr
        Vyo = summary[4+n_param*3][m]
        Vzo = summary[5+n_param*3][m]

        drift_params = (xo, yo, Vxo, Vyo, Vzo)
        x_bestfit, y_bestfit, z_bestfit, vx_bestfit, vy_bestfit, vz_bestfit, v_bestfit = get_orbit_prediction(elem, t, mass, drift_params)


        x_arr = np.zeros(len(lines))
        y_arr = np.zeros(len(lines))
        vz_arr = np.zeros(len(lines))
    
        for i in range(len(lines)):
            # MAP parameters
            elem = np.zeros(8)    
            elem[0] = arr[i][6]  # Distance        
            elem[2] = arr[i][10+6*star_num] # Period
            elem[3] = arr[i][12+6*star_num] # Eccentricity
            elem[4] = arr[i][11+6*star_num] #t0
            elem[5] = arr[i][8+6*star_num]  #w
            elem[6] = arr[i][9+6*star_num]  #i
            elem[7] = arr[i][7+6*star_num]  #Omega

            # Get a from M and Period
            mass = arr[i][0]
            #a_AU = (mass*elem[2]**2.)**(1./3.)   # a in AU = (Mass*Period^2)^(1/3)
            #a = a_AU/elem[0]  # a in arcsec
            #elem[1] = a
        
            #xy_table[1][i],xy_table[2][i],xy_table[3][i] = orbit_position(elem, t)
            #vx, vy, vz, v = orbit_vel(elem, t)

            # Add drift
            xo = arr[i][1]  # in arcsec
            yo = arr[i][2]
            Vxo = arr[i][3]  # in arcsec/yr
            Vyo = arr[i][4]
            Vzo = arr[i][5]

            drift_params = (xo, yo, Vxo, Vyo, Vzo)

            x_arr[i], y_arr[i], z, vx, vy, vz_arr[i], v = get_orbit_prediction(elem, t, mass, drift_params)

        x_stdev = weighted_stdev(x_arr,weights)
        y_stdev = weighted_stdev(y_arr,weights)
        vz_stdev = weighted_stdev(vz_arr,weights)

        out_table[1][j] = x_bestfit - x_stdev
        out_table[2][j] = x_bestfit + x_stdev
        out_table[3][j] = y_bestfit - y_stdev
        out_table[4][j] = y_bestfit + y_stdev
        out_table[5][j] = 0
        out_table[6][j] = 0

        out_table_vel[1][j] = 0
        out_table_vel[2][j] = 0
        out_table_vel[3][j] = 0
        out_table_vel[4][j] = 0
        out_table_vel[5][j] = vz_bestfit - vz_stdev
        out_table_vel[6][j] = vz_bestfit + vz_stdev

        #print out_table
        out_table.writeto(outDir+'orbit.'+star+'.model_lim_xy')
        out_table_vel.writeto(outDir+'orbit.'+star+'.model_lim_vel')

    out_table.writeto(outDir+'orbit.'+star+'.model_lim_xy')
    out_table_vel.writeto(outDir+'orbit.'+star+'.model_lim_vel')


