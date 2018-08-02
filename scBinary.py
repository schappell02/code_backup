from gcwork import objects
from gcwork import starset
import sc_accel_class as acc
from gcwork import util
from gcwork import orbits
from gcwork import young
from pysqlite2 import dbapi2 as sqlite
import MySQLdb as mysqldb
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
from decimal import *
import numpy as np
import pylab as py
import math
import histNofill
#import matplotlib.axes3d as p3
import pdb
import scipy.optimize as opter
from scipy.optimize import fsolve
from matplotlib.ticker import ScalarFormatter
import polyfit2 as pf2

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
rsun = 6.955e10
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


def PeriodPlot(mtotSolar=[1.,101.,20.],sepAU=[0.01,10.,0.01]):
    mtot = np.arange(mtotSolar[0],mtotSolar[1],mtotSolar[2])
    sep = np.arange(sepAU[0],sepAU[1],sepAU[2])

    py.clf()
    num = len(mtot)
    for i in range(num):
        mtmp = mtot[i]
        color = np.array([(-0.51*i/float(num))+0.51,0.51*i/float(num),(-0.51*i/float(num))+0.51])
        #purple to green
        ptmp = np.sqrt(sep**3 / mtmp) * 365.0
        py.plot(sep,ptmp,color=color)

#    py.plot([0.0,sepAU[1]],[10.0,10.0],'k')
    py.xlabel('Seperation (AU)')
    py.ylabel('Period (days)')
    py.title(str(mtotSolar[0])+' to '+str(mtotSolar[1])+' solar masses')
    py.savefig('/u/schappell/plots/period_seperation.png')

def AmpPlot(mtot,numM,sepAU=[0.01,10.,0.01]):
    sepAU = np.arange(sepAU[0],sepAU[1],sepAU[2])
    periodYR = np.sqrt(sepAU**3 / mtot)

    periodSEC = periodYR * sec_in_yr
    mtotG = msun *mtot
    py.clf()
    fig1 = py.figure()
    ax1 = fig1.add_subplot(111)

    for i in range(numM):
        m2 = (mtotG * float(i+1) / numM)
        color = np.array([(-1.*i/float(numM))+1.,(-0.45*i/float(numM))+0.45,i/float(numM)])
        #orange to blue
        #print color
        v1 = m2 * (2.0*pi*G/(periodSEC * mtotG**2))**(1./3.)
        atmp = 2165.* v1/2.99e10

        ax1.plot(periodYR*365.0,atmp,color=color)

    ax1.yaxis.tick_left()
    py.xlabel('Period (days)')
    py.ylabel('Amplitude (nm)')
    py.xlim(0.0,max(periodYR*365.0))
    cond = np.where(periodYR*365.0 > 0.1)[0]
    py.ylim(0.0,max(atmp[cond]))
    ax2 = fig1.add_subplot(111, sharex=ax1,frameon=False)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    py.ylabel('Radial Velocity (km/s)')
    py.ylim(0.0,max(v1[cond]*1e-5))
    py.title('Mass of Binary: '+str(mtot))
    py.savefig('/u/schappell/plots/amp_period_'+str(mtot)+'.png')



class RVstar(object):
    """
    For given star:RV, error, dates, number of frames, field, mag, r2d, and prob of being young

    vz - np array of radial velocities in km/s, taken from database, only consider those with errors
    vzerr - np array of error in radial velocities in km/s
    date - np array of dates in years of radial velocity measurements
    frames - np array of number of frames in given observation
    name - name of star given
    fieldarray - np array of field star was observed in
    field - field where star was most recently observed
    mag - K' magnitude of star
    r2d - projected distance from Sgr A* in arcsec
    yng - star's probability of being early-type

    plotRV - plots w/ errorbars the RV measurements of given star, set plotFig to False
             to have figure open on screen as opposed to being saved

    """
    def __init__(self, sname):
        self.vz=np.array([])
        self.vzerr=np.array([])
        self.date=np.array([])
        self.fieldarray=np.array([])
        self.frames=np.array([])
        self.name = sname
        self.mag = -1.0

        #connection to online database
        database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
        cur = database.cursor()

        #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
        #connection = sqlite.connect(dbfile)
    # Create a cursor object
        #cur = connection.cursor()
        cur.execute("SELECT ddate, vlsr, vz_err,field,nframes FROM spectra WHERE name='%s'"%sname)
        for row in cur:
            try:
                if ((row[2] > 0.0) & (row[1] != None) & (row[2] != None)):
                    self.vzerr=np.append(self.vzerr,np.float(row[2]))
                    self.vz=np.append(self.vz,np.float(row[1]))
                    self.date=np.append(self.date,float(row[0]))
                    self.fieldarray=np.append(self.fieldarray,row[3])
                    self.field=row[3]
                    self.frames=np.append(self.frames,row[4])
            except:
                pass
            
        #connection = sqlite.connect(dbfile)
        #cur = connection.cursor()

        
        cur = database.cursor()
        try:
            cur.execute("SELECT probYngSimPrior FROM unknownSims WHERE name='%s'" %(sname))
            for row in cur:
                self.yng = np.float(row[0])

        except:
            pass

        #connection = sqlite.connect(dbfile)
        #cur = connection.cursor()
        cur = database.cursor()
        self.cnt = len(self.vzerr)

        cur.execute("SELECT young,old,kp,r2d,x,y FROM stars WHERE name='%s'" %(sname))
        for row in cur:
            self.mag=row[2]
            self.r2d=row[3]
            self.x=row[4]
            self.y=row[5]
            self.max_az = round((GM * 1.0e-5 * sec_in_yr / (row[3]*dist*cm_in_au)**2),2)
            self.max_vz = round(sqrt(GM / (row[3]*dist*cm_in_au))*1.0e-5,2)
            if ((row[0] == 'T') | (row[1] == 'F')):
                self.yng = 1.0
            elif ((row[0] == 'F') | (row[1] == 'T')):
                self.yng = 0.0

        if ((sname == 'S0-38') | (sname == 'S0-49')):
            self.yng = 0.0
        if (sname=='S0-32'):
            self.yng = 1.0

        try:
            cur = database.cursor()
            cur.execute("SELECT t0_spectra,vz,vz_err FROM bartko2009 WHERE ucla_name='%s'" %(sname))
            for row in cur:
                self.b09_vz = row[1]
                self.b09_vzerr = row[2]
                self.b09_date = row[0]
        except:
            print str(sname)+' does not have a RV in Bartko 2009'

        try:
            cur = database.cursor()
            cur.execute("SELECT ddate,vz,vz_err FROM gillessen2017 WHERE ucla='%s'" %(sname))
            tmpdate = []
            tmprv = []
            tmperr = []
            for row in cur:
                if ((row[2]==row[2]) & (row[2]!=0)):
                    try:
                        tmprv = np.append(tmprv,float(row[1]))
                        tmperr = np.append(tmperr,float(row[2]))
                        tmpdate = np.append(tmpdate,float(row[0]))
                    except:
                        continue
            if (len(tmpdate)>0):
                self.g17_date = tmpdate
                self.g17_vz = tmprv
                self.g17_vzerr = tmperr
        except:
            print str(sname)+' does not have any points in Gillessen 2017'
        #Read in disk probabilities
        #disk_names = np.loadtxt('/u/schappell/diskmembers.dat',usecols=(0,),dtype='|S15')
        #disk_prob = np.loadtxt('/u/schappell/diskmembers.dat',delimiter='&',usecols=(12,))
        #for i in range(len(disk_names)):
        #   if (str(sname)==str(disk_names[i])):
        #       self.diskProb = disk_prob[i]





    def appendRV(self,addRV,addRVerr,addDate):
        self.vz = np.append(self.vz,addRV)
        self.vzerr = np.append(self.vzerr,addRVerr)
        self.date = np.append(self.date,addDate)
    




    def futureRV(self,future_yr=2017.38,wBartko=True,wGill=True):
        if ((len(self.vz) >= 2) | ((wBartko==True) & (hasattr(self,'b09_vz'))) | ((wGill==True) & (hasattr(self,'g17_vz')))):
            py.clf()
            tmpvz = np.append(self.vz,[])
            tmpvzerr = np.append(self.vzerr,[])
            tmpdate = np.append(self.date,[])
            if ((wBartko == True) & (hasattr(self,'b09_vz'))):
                try:
                    tmpvz = np.append(tmpvz,self.b09_vz)
                    tmpvzerr = np.append(tmpvzerr,self.b09_vzerr)
                    tmpdate = np.append(tmpdate,self.b09_date)
                    py.errorbar(self.b09_date,self.b09_vz,yerr=self.b09_vzerr,fmt='.',label='Bartko 2009')
                except:
                    print self.name+' does not have a RV in Bartko 2009'
            if ((wGill==True) & (hasattr(self,'g17_vz'))):
                try:
                    tmpvz = np.append(tmpvz,self.g17_vz)
                    tmpvzerr = np.append(tmpvzerr,self.g17_vzerr)
                    tmpdate = np.append(tmpdate,self.g17_date)
                    py.errorbar(self.g17_date,self.g17_vz,yerr=self.g17_vzerr,fmt='.',label='Gillessen 2017')
                except:
                    print self.name+' does not have a RV in Gillessen 2017'
            py.errorbar(self.date,self.vz,yerr=self.vzerr,fmt='.',label='UCLA')
            tmpt0 = np.average(tmpdate,weights=[1.0/e/e for e in tmpvzerr])
            v_terms, v_conv = pf2.polyfit2(tmpdate,tmpvz,deg=1,errors=tmpvzerr,t0=tmpt0)
            plot_date = np.array([np.min(tmpdate) + i*((future_yr-np.min(tmpdate))/100.0) for i in range(101)])
            plot_rv = v_terms[0] + v_terms[1]*(plot_date - tmpt0)
            if (len(tmpdate)>= 3):
                j_terms, j_conv = pf2.polyfit2(tmpdate,tmpvz,deg=2,errors=tmpvzerr,t0=tmpt0)
                plot_jerk = j_terms[0] + j_terms[1]*(plot_date - tmpt0) + j_terms[2]*(plot_date - tmpt0)**2
                future_jerk = j_terms[0] + j_terms[1]*(future_yr - tmpt0) + j_terms[2]*(future_yr - tmpt0)**2
                jerk_error = np.abs(j_conv[0,0]) + (plot_date - tmpt0)**2 * np.abs(j_conv[1,1]) + np.abs(plot_date - tmpt0)*(np.abs(j_conv[0,1]) + np.abs(j_conv[1,0]))
                jerk_error += (plot_date - tmpt0)**4 * np.abs(j_conv[2,2]) + (plot_date - tmpt0)**2 * (np.abs(j_conv[0,2]) + np.abs(j_conv[2,0])) + np.abs(plot_date - tmpt0)**3 * (np.abs(j_conv[1,2]) + np.abs(j_conv[2,1]))
                jerk_error = np.sqrt(jerk_error)
                future_jverr = np.abs(j_conv[0,0]) + (future_yr - tmpt0)**2 * np.abs(j_conv[1,1]) + np.abs(future_yr - tmpt0)*(np.abs(j_conv[0,1]) + np.abs(j_conv[1,0]))
                future_jverr += (future_yr - tmpt0)**4 * np.abs(j_conv[2,2]) + (future_yr - tmpt0)**2 * (np.abs(j_conv[0,2]) + np.abs(j_conv[2,0])) + np.abs(future_yr - tmpt0)**3 * (np.abs(j_conv[1,2]) + np.abs(j_conv[2,1]))
                future_jverr = math.sqrt(future_jverr)
            else:
                future_jerk = 0.0
                future_jverr = 0.0
            future_vz = v_terms[0] + v_terms[1]*(future_yr - tmpt0)
            plot_error = np.sqrt(v_conv[0,0] + (plot_date - tmpt0)**2 * v_conv[1,1] + np.abs(plot_date - tmpt0)*(v_conv[0,1] + v_conv[1,0]))
            future_vzerr = math.sqrt(v_conv[0,0] + (future_yr - tmpt0)**2 * v_conv[1,1] + np.abs(future_yr - tmpt0)*(v_conv[0,1] + v_conv[1,0]))
            py.errorbar(future_yr,future_vz,yerr=future_vzerr,fmt='.',label='Predicted Accel')
            if (len(tmpdate)>= 3):
                py.errorbar(future_yr,future_jerk,yerr=future_jverr,fmt='.',label='Predicted Jerk')
                py.plot(plot_date,plot_jerk,'r',label='Best Fit Jerk')
                py.plot(plot_date,plot_jerk+jerk_error,'r:')
                py.plot(plot_date,plot_jerk-jerk_error,'r:',label='Jerk Error')
            py.plot(plot_date,plot_rv,'k',label='Best Fit Accel')
            py.plot(plot_date,plot_rv+plot_error,'k--')
            py.plot(plot_date,plot_rv-plot_error,'k--',label='Accel Error')
            py.legend()
            py.savefig('/u/schappell/plots/'+str(future_yr)+'_'+str(self.name)+'_RV.png')
            py.clf()
            #let's run the F test
            pass_ftest_aj = 0.0
            Fval_aj = 1e10
            if (len(tmpdate)> 3):
                pred_a = v_terms[0] + v_terms[1]*(tmpdate - tmpt0)
                pred_j = j_terms[0] + j_terms[1]*(tmpdate - tmpt0) + j_terms[2]*(tmpdate - tmpt0)**2
                chi2_a = np.sum((pred_a - tmpvz)**2/tmpvzerr/tmpvzerr)
                chi2_j = np.sum((pred_j - tmpvz)**2/tmpvzerr/tmpvzerr)
                dof_j = float(len(tmpvzerr)) - 3.0
                Fval_aj = ((chi2_a - chi2_j)/1.0)/(chi2_j/dof_j)
                Fprob_aj = stats.f.sf(Fval_aj,1,dof_j)
                signif = scipy.special.erfc(4.0/math.sqrt(2.0))
                if (Fprob_aj < signif):
                    pass_ftest_aj = 1.0
            
            
            rdex = np.argmax(tmpdate)
            recent_rv = tmpvz[rdex]
            recent_rverr = tmpvzerr[rdex]
            recent_date = tmpdate[rdex]

            return future_vz,future_vzerr,recent_date,recent_rv,recent_rverr,future_jerk,future_jverr, pass_ftest_aj, Fval_aj
        else:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0





    def accel(self,wBartko=True):
        if ((len(self.vz) >= 2) | ((wBartko==True) & (hasattr(self,'b09_vz')))):
            self.az = np.array([])
            self.azerr = np.array([])
            tmpvz = np.append(self.vz,[])
            tmpvzerr = np.append(self.vzerr,[])
            tmpdate = np.append(self.date,[])
            if (wBartko == True):
                try:
                    tmpvz = np.append(tmpvz,self.b09_vz)
                    tmpvzerr = np.append(tmpvzerr,self.b09_vzerr)
                    tmpdate = np.append(tmpdate,self.b09_date)
                except:
                    print self.name+' does not have a RV in Bartko 2009'
            adex = np.argsort(tmpdate)
            self.binary_cand = False
            for i in range(len(tmpdate)-1):
                for ii in np.arange(i+1,len(tmpdate),1):
                    j = adex[i]
                    jp = adex[ii]
                    self.az = np.append(self.az,(tmpvz[jp]-tmpvz[j])/(tmpdate[jp]-tmpdate[j]))
                    self.azerr = np.append(self.azerr,sqrt(tmpvzerr[jp]**2+tmpvzerr[j]**2)/(tmpdate[jp]-tmpdate[j]))
                    if ((self.az[-1] < -1.0*(self.max_az+self.azerr[-1])) | (self.az[-1] > (self.max_az+self.azerr[-1]))):
                        print self.name+' has accel in z between dates '+str(tmpdate[j])+' and '+str(tmpdate[jp])+' larger than az_max'
                        print str(self.az[-1])+' +/- '+str(self.azerr[-1])+' km/s/yr'
                        #print 'Sigma: '+str(round(abs(abs(self.az[-1])-self.max_az)/self.azerr[-1],2))+'     (|az - az_max|/az_err)'
                        self.binary_cand = True
            for i in range (len(self.az)-1):
                for j in np.arange(i+1,len(self.az),1):
                    if (abs(self.az[i] - self.az[j]) > (self.azerr[i]+self.azerr[j])):
                        print self.name+' has accelerations that are NOT consistent'
                        break
                    

            self.ext_mass = np.fabs(self.az)*1.0e5*(self.r2d*dist*cm_in_au)**2/(sec_in_yr*G*msun) #in solar masses
            self.ext_mass -= mass #subtract out Sgr A*
            self.EM_err = self.azerr*1.0e5*(self.r2d*dist*cm_in_au)**2/(sec_in_yr*G*msun)
            
            maxdex = np.argmax(np.fabs(self.az))
            print self.name
            print 'Max a_z: '+str(round(self.az[maxdex],2))+' +/- '+str(round(self.azerr[maxdex],2))+' km/s/yr'
            print 'Max a_z from Sgr A*: '+str(round(self.max_az,2))+' km/s/yr'
            maxdex = np.argmax(self.ext_mass)
            if ((self.ext_mass[maxdex] - self.EM_err[maxdex]) <= 0.0):
                print 'Change in RV can be explained by Sgr A* and z=0 within error'
                #self.binary_cand = True
            else:
               print 'Extended mass: %.2E +/- %.2E solar masses' %(Decimal(str(self.ext_mass[maxdex])),
                                                                   Decimal(str(self.EM_err[maxdex])))
               #self.binary_cand = False

        else:
            print 'Has less than 2 RV entries, extended mass cannot be calculated'
            self.binary_cand = False



    def plotRV(self,saveFig=True,wBartko=False,maxErr=1.0e6,spectra=''):
        if ((len(self.vz) > 1) | (hasattr(self,'b09_vz') & (wBartko==True))):
            py.clf()
            usedex=np.where(self.vzerr<=maxErr)[0]
            if (spectra!=''):
                py.figure(1)
                py.subplot(212)
            py.errorbar(self.date[usedex],self.vz[usedex],yerr=self.vzerr[usedex],fmt='b*',ms=15,label='Keck')
            if (self.name=='irs13E1'):
                py.errorbar([2009.3881],[134],yerr=10.0,fmt='r^',ms=10)
                #py.errorbar([2003.5],[71.0],yerr=20.0,xerr=1.5,fmt='r^',ms=10)
            if (wBartko==True):
                try:
                    if (self.b09_vzerr <= maxErr):
                        py.errorbar(self.b09_date,self.b09_vz,yerr=self.b09_vzerr,fmt='r^',ms=10,label='VLT')
                        py.legend(numpoints=1)
                except:
                    pass         
            py.xlabel('Time (years)')
            py.ylabel('RV (km/s)')
            if (spectra==''):
                py.xlim([2000,2017])
                py.title('Radius: '+str(self.r2d)+' as, Kp: '+str(self.mag)+', Prob Young: '+str(self.yng)+', max a$_Z$: '+
                         str(self.max_az)+' km/s/yr')
            else:
                #py.xlim([2007,2011])
                py.xticks([2007,2008,2009,2010,2011])
                locs,labels = xticks()
                xticks(locs,map(lambda x: "%g" % x, locs))
                #py.gca().ticklabel_format(style='plain',useOffset=False)
                py.yticks([-40,0,40,80,120,160])
                spec_array = np.loadtxt(spectra)
                py.subplot(211)
                py.plot(spec_array[:,0]*1000.0,spec_array[:,1])
                py.xlabel('Wavelength (nm)')
                py.ylabel('Flux')
                py.xlim([2120,2230])
            if (saveFig==False):
                py.show()
            else:
                py.savefig('/u/schappell/plots/RV_'+self.name+'.png')
        else:
            print 'Only 1 point, not enough to plot'


    def bestFitLine(self,plotFit=True,wBartko=True,maxErr=1.0e5):
        edex = np.where(self.vzerr<=maxErr)[0]
        use_vz = self.vz[edex]
        use_vzerr = self.vzerr[edex]
        use_date = self.date[edex]
        if (wBartko==True):
            try:
                use_vz = np.append(use_vz,self.b09_vz)
                use_vzerr = np.append(use_vzerr,self.b09_vzerr)
                use_date = np.append(use_date,self.b09_date)
            except:
                pass
        usedex = np.argsort(use_date)
        use_date = use_date[usedex]
        use_vz = use_vz[usedex]
        use_vzerr = use_vzerr[usedex]
        if (len(use_vz) > 2):
            lsq_res = opter.leastsq(chi2_line,[self.max_az,use_vz[0]],args=(use_date,use_vz,use_vzerr),full_output=1)
            coeff = lsq_res[0]
            best_fit = use_date * coeff[0] + coeff[1]
            residuals = use_vz - best_fit
            red_chi2 = np.sum(residuals**2 / use_vzerr**2 / (len(use_date) - 2.0))
            print 'Star: '+str(self.name)
            print 'Mag: '+str(self.mag)
            print 'R2d: '+str(self.r2d)+' arcsec'
            print 'Prob Young: '+str(self.yng)
            print 'Max a_z: '+str(self.max_az)+' km/s/yr'
            print ''
            print 'Best fit:'
            print 'SLOPE: '+str(coeff[0])+' km/s/yr'
            print 'Y intercept: '+str(coeff[1])+' km/s'
            print 'Reduced Chi^2: '+str(red_chi2)

            if (plotFit==True):
                py.clf()
                py.plot(use_date,best_fit,'k')
                py.errorbar(use_date,use_vz,yerr=use_vzerr,fmt='b.')
                if (wBartko==True):
                    try:
                        py.errorbar(self.b09_date,self.b09_vz,yerr=self.b09_vzerr,fmt='r.',label='Bartko 2009')
                        py.legend()
                    except:
                        pass
                py.xlabel('Time (years)')
                py.xlim([2000,2017])
                py.ylabel('RV (km/s)')
                py.title('Best fit a$_Z$: '+str(round(coeff[0],2))+' km/s/yr    Max a$_Z$: '+str(self.max_az)+' km/s/yr    Reduced $\chi^2$: '+str(round(red_chi2,2)))
                py.savefig('/u/schappell/plots/RV_'+str(self.name)+'_bestfitLINE.png')

                py.clf()
                py.plot([2000,2017],[0.0,0.0],'k')
                py.errorbar(use_date,residuals,yerr=use_vzerr,fmt='b.')
                if (wBartko==True):
                    try:
                        bdex = np.where(use_vzerr==self.b09_vzerr)[0]
                        py.errorbar(self.b09_date,residuals[bdex],yerr=self.b09_vzerr,fmt='r.',label='Bartko 2009')
                        py.legend()
                    except:
                        pass
                py.xlabel('Time (years)')
                py.xlim([2000,2017])
                py.ylabel('RV Residuals (km/s)')
                py.title('Best fit a$_Z$: '+str(round(coeff[0],2))+' km/s/yr    Max a$_Z$: '+str(self.max_az)+r' km/s/yr    Reduced $\chi^2$: '+str(round(red_chi2,2)))
                py.savefig('/u/schappell/plots/RV_'+str(self.name)+'_bfl_res.png')

            return coeff,use_vz,use_vzerr,use_date

        else:
            print 'Only 2 points or less, not enough to fit'


    def bestFitLine_sin(self,plotFit=True,guess=[30.0,10.0,0.0],maxErr=1.0e5,wBartko=True):
        coeff,use_vz,use_vzerr,use_date=self.bestFitLine(plotFit=plotFit,wBartko=wBartko,maxErr=maxErr)
        line_vz = use_date * coeff[0] + coeff[1]
        line_res = use_vz - line_vz

        if (len(use_vz) > 3):
            lsq_res = opter.leastsq(chi2_sin,guess,args=(use_date,line_res,use_vzerr),full_output=1)
            coeff_sin = lsq_res[0]
            best_fit = use_date * coeff[0] + coeff[1] + coeff_sin[0] * np.sin(2.0*pi*coeff_sin[1] * use_date - coeff_sin[2])
            residuals = use_vz - best_fit
            red_chi2 = np.sum(residuals**2 / use_vzerr**2 / (len(use_date) - 3.0))

            print ''
            print 'Best fit SIN:'
            print 'SIN AMP '+str(coeff_sin[0])+' km/s'
            print 'PERIOD '+str(1.0/coeff_sin[1])+' years'
            print 'Offset '+str(coeff_sin[2])+' radians'
            
            print 'Reduced Chi^2: '+str(red_chi2)

            if (plotFit==True):
                py.clf()
                py.plot(self.date,best_fit,'bs')
                py.errorbar(self.date,self.vz,yerr=self.vzerr,fmt='g.')
                py.xlabel('Time (years)')
                py.xlim([2000,2016])
                py.ylabel('RV (km/s)')
                py.title('Best fit a$_Z$: '+str(round(coeff[0],2))+' km/s/yr, Max a$_Z$: '+str(self.max_az)+' km/s/yr, Period: '
                         +str(round(1.0/coeff_sin[1],2))+' years')
                py.savefig('/u/schappell/plots/RV_'+str(self.name)+'_bestfitLINE_sin.png')

                py.clf()
                py.plot([2000,2016],[0.0,0.0])
                py.errorbar(self.date,residuals,yerr=self.vzerr,fmt='.')
                py.xlabel('Time (years)')
                py.xlim([2000,2016])
                py.ylabel('RV Residuals (km/s)')
                py.title('Best fit a$_Z$: '+str(round(coeff[0],2))+' km/s/yr, Max a$_Z$: '+str(self.max_az)+' km/s/yr, Period: '
                         +str(round(1.0/coeff_sin[1],2))+' years')
                py.savefig('/u/schappell/plots/RV_'+str(self.name)+'_bfl_sin_res.png')
        else:
            print 'Only 3 or less points, not enough to fit sin'


    
    def toDat(self,wBartko=True,label=''):
        use_vz = np.append(self.vz,[])
        use_vzerr = np.append(self.vzerr,[])
        use_date = np.append(self.date,[])
        if (wBartko==True):
            try:
                use_vz = np.append(use_vz,self.b09_vz)
                use_vzerr = np.append(use_vzerr,self.b09_vzerr)
                use_date = np.append(use_date,self.b09_date)
            except:
                pass

        np.savetxt('/u/schappell/'+label+'date_vz_err.dat',np.transpose([use_date,use_vz,use_vzerr]),delimiter=' ')




    def plot2017SinC(self,saveFig=True,wBartko=False,maxErr=1.0e6,period=2.0/365.25,phase0=0.0,RVmax=300.0,RV_os=0.0,eccen=0.0):
        if ((len(self.vz) > 1) | (hasattr(self,'b09_vz') & (wBartko==True))):
            py.clf()
            usedex=np.where(self.vzerr<=maxErr)[0]
            phase=self.date[usedex]*2.0*pi/period
            phase -= int(phase/2.0/pi)*2.0*pi
            py.errorbar(phase,self.vz[usedex],yerr=self.vzerr[usedex],fmt='b*',ms=15,label='Keck')
            #use_dates = np.append(self.date[usedex],[])
            #use_vz = np.append(self.vz[usedex],[])
            #err17 = np.zeros(len(dates17))+err17
            #use_error = np.append(self.vzerr[usedex],err17)
            if (self.name=='irs13E1'):
                phase=2009.3881*2.0*pi/period
                phase -= int(phase/2.0/pi)*2.0*pi
                py.errorbar(phase,[134],yerr=10.0,fmt='r^',ms=10)
                #py.errorbar([2003.5],[71.0],yerr=20.0,xerr=1.5,fmt='r^',ms=10)
                #use_dates=np.append(use_dates,2009.3881)
                #use_vz=np.append(use_vz,134.0)
                #use_error=np.append(use_error,[10.0,20.0])
            if (wBartko==True):
                try:
                    if (self.b09_vzerr <= maxErr):
                        phase=self.b09_date*2.0*pi/period
                        phase -= int(phase/2.0/pi)*2.0*pi
                        py.errorbar(phase,self.b09_vz,yerr=self.b09_vzerr,fmt='r^',ms=10,label='VLT')
                        py.legend(numpoints=1)
                        #use_dates=np.append(use_dates,self.b09_date)
                        #use_vz=np.append(use_vz,self.b09_vz)
                        #use_error=np.append(use_error,self.b09_vzerr)
                except:
                    pass
            #usedex = np.argsort(use_dates)
            #use_dates = use_dates[usedex]
            #use_vz = use_vz[usedex]

            #phase0 = (pi/2.0) - 2.0*pi*max_year/period
            #amp = (use_vz[-2] - use_vz[-1]) / 2.0 / math.sin(pi*(use_dates[-2]-use_dates[-1])/period) / math.cos(phase0 + pi*(use_dates[-2]+use_dates[-1])/period)
            #amp = abs(amp)
            #RV0 = use_vz[-1] - amp*math.sin(phase0 + 2.0*pi*use_dates[-1]/period)

            phases = np.linspace(-0.5/pi,2.0*pi,1000.0)
            RVs = RVmax * np.sin(phases+phase0) * np.sqrt(1.0+eccen*np.cos(phases+phase0)) + RV_os

            #use_dates = np.append(use_dates,dates17)
            #use_vz = amp*np.sin(phase0 + 2.0*pi*use_dates/period) + RV0
            #print use_dates
            #print use_vz
            #print ''
            #print 'Amplitude: '+str(amp)+' km/s'
            #print 'RV offset: '+str(RV0)+' km/s'
            #py.plot(use_dates,use_vz,'kD',ms=10)
            py.plot(phases,RVs,'k')
            #for i in range(len(dates17)):
            #    py.plot([dates17[i],dates17[i]],[-1000.0,1000.0],'k--')
            py.xlabel('Phase')
            #py.xlim([2007.5,2018])
            py.xlim([-0.1,2.0*pi])
            pymax = round(abs(RVmax) + RV_os + 20.0)
            pymin = round(-1.0*abs(RVmax) + RV_os - 20.0)
            py.ylim([pymin,pymax])
            py.ylabel('RV (km/s)')
            #py.title('Radius: '+str(self.r2d)+' as, Kp: '+str(self.mag)+', Prob Young: '+str(self.yng)+', max a$_Z$: '+
             #        str(self.max_az)+' km/s/yr')
            py.title(str(self.name))
            if (saveFig==False):
                py.show()
            else:
                py.savefig('/u/schappell/plots/RV_'+self.name+'_phase.png')
            #pdb.set_trace()
        else:
            print 'Only 1 point, not enough to plot'



    def plot2017SinB(self,saveFig=True,wBartko=False,maxErr=1.0e6,dates17=[2017.4,2017.4027,2017.4438,2017.4465],
                     period=2.0/365.25,phase0=0.0,RVmax=300.0,RV_os=0.0,eccen=0.0,):#max_year=2009.3881):
        if ((len(self.vz) > 1) | (hasattr(self,'b09_vz') & (wBartko==True))):
            py.clf()
            usedex=np.where(self.vzerr<=maxErr)[0]
            py.errorbar(self.date[usedex],self.vz[usedex],yerr=self.vzerr[usedex],fmt='b*',ms=15,label='Keck')
            #use_dates = np.append(self.date[usedex],[])
            #use_vz = np.append(self.vz[usedex],[])
            #err17 = np.zeros(len(dates17))+err17
            #use_error = np.append(self.vzerr[usedex],err17)
            if (self.name=='irs13E1'):
                py.errorbar([2009.3881],[134],yerr=10.0,fmt='r^',ms=10)
                #py.errorbar([2003.5],[71.0],yerr=20.0,xerr=1.5,fmt='r^',ms=10)
                #use_dates=np.append(use_dates,2009.3881)
                #use_vz=np.append(use_vz,134.0)
                #use_error=np.append(use_error,[10.0,20.0])
            if (wBartko==True):
                try:
                    if (self.b09_vzerr <= maxErr):
                        py.errorbar(self.b09_date,self.b09_vz,yerr=self.b09_vzerr,fmt='r^',ms=10,label='VLT')
                        py.legend()
                        #use_dates=np.append(use_dates,self.b09_date)
                        #use_vz=np.append(use_vz,self.b09_vz)
                        #use_error=np.append(use_error,self.b09_vzerr)
                except:
                    pass
            #usedex = np.argsort(use_dates)
            #use_dates = use_dates[usedex]
            #use_vz = use_vz[usedex]

            #phase0 = (pi/2.0) - 2.0*pi*max_year/period
            #amp = (use_vz[-2] - use_vz[-1]) / 2.0 / math.sin(pi*(use_dates[-2]-use_dates[-1])/period) / math.cos(phase0 + pi*(use_dates[-2]+use_dates[-1])/period)
            #amp = abs(amp)
            #RV0 = use_vz[-1] - amp*math.sin(phase0 + 2.0*pi*use_dates[-1]/period)

            dates = np.linspace(2000.0,2018.0,1000.0)
            phases = phase0 + 2.0*pi*dates/period
            RVs = RVmax * np.sin(phases) * np.sqrt(1.0+eccen*np.cos(phases)) + RV_os

            #use_dates = np.append(use_dates,dates17)
            #use_vz = amp*np.sin(phase0 + 2.0*pi*use_dates/period) + RV0
            #print use_dates
            #print use_vz
            #print ''
            #print 'Amplitude: '+str(amp)+' km/s'
            #print 'RV offset: '+str(RV0)+' km/s'
            #py.plot(use_dates,use_vz,'kD',ms=10)
            py.plot(dates,RVs,'k')
            for i in range(len(dates17)):
                py.plot([dates17[i],dates17[i]],[-1000.0,1000.0],'k--')
            py.xlabel('Time (years)')
            py.xlim([2007.5,2018])
            pymax = round(abs(RVmax) + RV_os + 20.0)
            pymin = round(-1.0*abs(RVmax) + RV_os - 20.0)
            py.ylim([pymin,pymax])
            py.ylabel('RV (km/s)')
            py.title('Radius: '+str(self.r2d)+' as, Kp: '+str(self.mag)+', Prob Young: '+str(self.yng)+', max a$_Z$: '+
                     str(self.max_az)+' km/s/yr')
            if (saveFig==False):
                py.show()
            else:
                py.savefig('/u/schappell/plots/RV_'+self.name+'.png')
            #pdb.set_trace()
        else:
            print 'Only 1 point, not enough to plot'


    def plot2017Sin(self,saveFig=True,wBartko=False,maxErr=1.0e6,dates17=[2017.4,2017.4027,2017.4438,2017.4465],
                    period=2.0/365.25,phase0=0.0):#max_year=2009.3881):
        if ((len(self.vz) > 1) | (hasattr(self,'b09_vz') & (wBartko==True))):
            py.clf()
            usedex=np.where(self.vzerr<=maxErr)[0]
            py.errorbar(self.date[usedex],self.vz[usedex],yerr=self.vzerr[usedex],fmt='b*',ms=15,label='Keck')
            use_dates = np.append(self.date[usedex],[])
            use_vz = np.append(self.vz[usedex],[])
            #err17 = np.zeros(len(dates17))+err17
            #use_error = np.append(self.vzerr[usedex],err17)
            if (self.name=='irs13E1'):
                py.errorbar([2009.3881],[134],yerr=10.0,fmt='r^',ms=10)
                #py.errorbar([2003.5],[71.0],yerr=20.0,xerr=1.5,fmt='r^',ms=10)
                use_dates=np.append(use_dates,2009.3881)
                use_vz=np.append(use_vz,134.0)
                #use_error=np.append(use_error,[10.0,20.0])
            if (wBartko==True):
                try:
                    if (self.b09_vzerr <= maxErr):
                        py.errorbar(self.b09_date,self.b09_vz,yerr=self.b09_vzerr,fmt='r^',ms=10,label='VLT')
                        py.legend(loc=2)
                        use_dates=np.append(use_dates,self.b09_date)
                        use_vz=np.append(use_vz,self.b09_vz)
                        #use_error=np.append(use_error,self.b09_vzerr)
                except:
                    pass
            usedex = np.argsort(use_dates)
            use_dates = use_dates[usedex]
            use_vz = use_vz[usedex]

            #phase0 = (pi/2.0) - 2.0*pi*max_year/period
            amp = (use_vz[-2] - use_vz[-1]) / 2.0 / math.sin(pi*(use_dates[-2]-use_dates[-1])/period) / math.cos(phase0 + pi*(use_dates[-2]+use_dates[-1])/period)
            amp = abs(amp)
            RV0 = use_vz[-1] - amp*math.sin(phase0 + 2.0*pi*use_dates[-1]/period)

            use_dates = np.append(use_dates,dates17)
            use_vz = amp*np.sin(phase0 + 2.0*pi*use_dates/period) + RV0
            print use_dates
            print use_vz
            print ''
            print 'Amplitude: '+str(amp)+' km/s'
            print 'RV offset: '+str(RV0)+' km/s'
            py.plot(use_dates,use_vz,'kD',ms=10)
            py.xlabel('Time (years)')
            py.xlim([2000,2018])
            #py.ylim([-100,200])
            py.ylabel('RV (km/s)')
            py.title('Radius: '+str(self.r2d)+' as, Kp: '+str(self.mag)+', Prob Young: '+str(self.yng)+', max a$_Z$: '+
                     str(self.max_az)+' km/s/yr')
            if (saveFig==False):
                py.show()
            else:
                py.savefig('/u/schappell/plots/RV_'+self.name+'.png')
        else:
            print 'Only 1 point, not enough to plot'


    def mockObs(self,wBartko=True,m1=40.0,q=0.5,period=6.0/365.0,eccen=0.0,incl=90.0,offset=0.0,plotFit=False,guess=[20.0,36.5,0.0],
                amp_bounds=[-1000.0,1000.0],f_bounds=[0.1,1.0e4],offset_bounds=[0.0,3.15]):
        #m1 is mass of larger one, in solar masses
        #max RV measured from 2nd (smaller or equal mass) star
        #offset and offset_sgr are offsets in sine motion of binary and in orbit around Sgr A*, in radians
        #q = m1/m2
        #period in years
        #inclination and phi in degress, inclination is that of binary system to line of sight (90 = edge on)
        #phi is inclination of binary's orbit around Sgr A*
        #take R2d of actual star as R2d of mock binary 

       # self.period = period ; self.m1 = m1 ; self.m2 = q * m1 ; self.q = q ; self.e = eccen ; self.i = incl ; self.offset = offset
        #self.phi = phi ; self.offset_sgr = offset_sgr

        period_s = period * sec_in_yr
        incl_rad = incl * pi / 180.0
        m1_g = m1 * msun
        #phi_rad = phi * pi / 180.0
        
        #r3d_cm = self.r2d * dist * cm_in_au / math.cos(phi_rad)
        #self.r3d_pc = r3d_cm / cm_in_pc
        #self.r3d_as = r3d_cm / dist / cm_in_au
        amp = (2.0*pi*G/(period_s*(m1_g+m1_g*q)**2))**(1.0/3.0) * m1_g * math.sin(incl_rad) / 1.0e5 #in km/s
        #amp_sgr = math.sqrt(GM/r3d_cm) * np.sin(phi_rad) / 1.0e5 #in km/s
        #period_sgr = 2.0 * pi * r3d_cm**(3.0/2.0) / math.sqrt(GM) / sec_in_yr

        #self.amp = amp ; self.amp_sgr = amp_sgr ; self.period_sgr = period_sgr

        use_vzerr = np.append(self.vzerr,[])
        use_date = np.append(self.date,[])
        if (wBartko==True):
            try:
                use_vzerr = np.append(use_vzerr,self.b09_vzerr)
                use_date = np.append(use_date,self.b09_date)
            except:
                pass

        print ''
        num_mock = raw_input("How many RV obs would you like to add? ") ; num_mock = int(num_mock)
        mock_date = np.zeros(num_mock)

        for i in range(num_mock):
            if (i==0):
                first_date = raw_input("When is the first new observation (enter year in decimal)? ")
                mock_date[0] = float(first_date)
            else:
                tmp_diff = raw_input("What is difference between last and next observation (enter years in decimal)? ")
                mock_date[i] = mock_date[i-1] + float(tmp_diff)

        mock_vzerr = raw_input("What is the error on the RV measurement(s) (in km/s)? ")
        mock_vzerr = np.zeros(len(mock_date)) + float(mock_vzerr)
        #self.mock_vzerr = mock_vzerr ; self.mock_date = mock_date

        use_date = np.append(use_date,mock_date) ; use_vzerr = np.append(use_vzerr,mock_vzerr)

        phase = (use_date * 2.0 * pi / period) + offset #in radians
        #phase_sgr = (use_date * 2.0 * pi / period_sgr) + offset #in radians
        vzmock=np.zeros(len(use_date))
        for i in range(len(use_date)):
            #true_val = amp * np.sin(phase[i])*np.sqrt(1.0+eccen*np.cos(phase[i])) + amp_sgr * np.sin(phase_sgr[i])
            true_val = amp * np.sin(phase[i])*np.sqrt(1.0+eccen*np.cos(phase[i]))
            vzmock[i] = np.random.normal(true_val,use_vzerr[i])
        #self.phase = phase ; self.phase_sgr = phase_sgr ; self.vzmock = vzmock


        #lsq_res = opter.leastsq(chi2_line,[self.max_az,vzmock[0]],args=(use_date,vzmock,use_vzerr),full_output=1)
        #coeff = lsq_res[0]
        #best_fit = use_date * coeff[0] + coeff[1]
        #residuals = vzmock - best_fit
        #red_chi2 = np.sum(residuals**2 / use_vzerr**2 / (len(use_date) - 2.0))
        #print ' '
        #print 'Star: '+str(self.name)
        #print 'Mag: '+str(self.mag)
        #print 'R2d: '+str(self.r2d)+' arcsec'
        #print 'Prob Young: '+str(self.yng)
        #print 'Max a_z: '+str(self.max_az)+' km/s/yr'
        #print ''
        #print 'Best fit:'
        #print 'SLOPE: '+str(coeff[0])+' km/s/yr'
        #print 'Y intercept: '+str(coeff[1])+' km/s'
        #print 'Reduced Chi^2: '+str(red_chi2)
        #print ' '
        #az0 = GM * sec_in_yr * math.sin(phi_rad) / r3d_cm**2 / 1.0e5
        #print 'Actual values:'
        #print 'SLOPE: '+str(az0)+' km/s/yr'
        #print 'Y intercept: '+str(amp_sgr)+' km/s'


        if (len(vzmock) > 3):
            #lsq_res = opter.leastsq(chi2_sin,guess,args=(use_date,residuals,use_vzerr),full_output=1)
            lsq_res = opter.leastsq(chi2_sin,guess,args=(use_date,vzmock,use_vzerr),full_output=1)
            coeff_sin = lsq_res[0]
            #best_fit_sin = use_date * coeff[0] + coeff[1] + (coeff_sin[0] * np.sin(2.0*pi*coeff_sin[1] * use_date - coeff_sin[2]))
            best_fit_sin = coeff_sin[0] * np.sin(2.0*pi*coeff_sin[1] * use_date - coeff_sin[2])
            residuals_sin = vzmock - best_fit_sin
            red_chi2_sin = np.sum(residuals_sin**2 / use_vzerr**2 / (len(use_date) - 3.0))

            print ''
            print 'Best fit SIN:'
            print 'SIN AMP '+str(coeff_sin[0])+' km/s'
            print 'PERIOD '+str(1.0/coeff_sin[1])+' years'
            print 'Offset '+str(coeff_sin[2])+' radians'
            print 'Reduced Chi^2: '+str(red_chi2_sin)
            print ' '
            print 'Actual values:'
            print 'SIN AMP '+str(amp)+' km/s'
            print 'PERIOD '+str(period)+' years'
            print 'Offset '+str(offset)+' radians'
            
            np.savetxt('/u/schappell/mockRV.dat',np.transpose([use_date,vzmock,use_vzerr]),delimiter=' ')

            if (plotFit==True):
                py.clf()
                pl_date=np.linspace(min(use_date),max(use_date),10000)
                pl_phase = (pl_date * 2.0 * pi / period) + offset #in radians
                #pl_phase_sgr = (pl_date * 2.0 * pi / period_sgr) + offset #in radians
                #pl_vz = amp * np.sin(pl_phase)*np.sqrt(1.0+eccen*np.cos(pl_phase)) + amp_sgr * np.sin(pl_phase_sgr)
                pl_vz = amp * np.sin(pl_phase)*np.sqrt(1.0+eccen*np.cos(pl_phase))
                #py.plot(use_date,best_fit,'r',label='Line fit')
                py.plot(use_date,best_fit_sin,'bs',label='Sine fit')
                py.plot(pl_date,pl_vz,'r')
                py.errorbar(use_date,vzmock,yerr=use_vzerr,fmt='g.',label='Mock data')
                py.legend(loc=2)
                py.xlabel('Time (years)')
                py.xlim([2000,2018])
                py.ylabel('RV (km/s)')

                #py.title('Best fit a$_Z$: '+str(round(coeff[0],2))+' km/s/yr, Actual a$_Z$: '+str(round(az0,2))+' km/s/yr, Best fit P: '
                #         +str(round(1.0/coeff_sin[1],2))+' years,  Actual P: '+str(round(period,2))+' years')
                py.savefig('/u/schappell/plots/RV_mock_bestfitLINE_sin.png')

                py.clf()
                py.plot([2000,2018],[0.0,0.0],'k')
                #py.errorbar(use_date,residuals,yerr=use_vzerr,fmt='r.',label='Line Fit')
                py.errorbar(use_date,residuals_sin,yerr=use_vzerr,fmt='b.',label='Sine fit')
                py.xlabel('Time (years)')
                py.xlim([2000,2018])
                py.ylabel('RV Residuals (km/s)')
                py.legend(loc=2)
                #py.title('Best fit a$_Z$: '+str(round(coeff[0],2))+' km/s/yr, Actual a$_Z$: '+str(round(az0,2))+' km/s/yr, Best fit P: '
                #         +str(round(1.0/coeff_sin[1],2))+' years,  Actual P: '+str(round(period,2))+' years')
                py.savefig('/u/schappell/plots/RV_mock_bfl_sin_res.png')
        else:
            print 'Have 3 or less observation points, not enough for sine fit'
        #pdb.set_trace()




def readLS(label='',savePlot=True):
    #use RVstar.toDat to save date, vz, and vzerr to .dat file
    #then run scPMN.outLombScargle, outputs .dat with period in years and power

    LSout = np.loadtxt('/u/schappell/'+label+'LS_period_power.dat')

    period = LSout[:,0]
    power = LSout[:,1]

    py.clf()
    py.plot(period,power)
    py.xlabel('Period (years)')
    py.ylabel('Power')
    py.xscale('log')
    py.ylim([0.0,1.0])
    #py.title('Radius: '+str(self.r2d)+' as, Kp: '+str(self.mag)+', Prob Young: '+str(self.yng)+', max a$_Z$: '+
    #         str(self.max_az)+' km/s/yr')
    if (savePlot==True):
        py.savefig('/u/schappell/plots/'+label+'LS_period.png')
    else:
        py.show()
        pdb.set_trace()

        

def chi2_line(constants, xvar, yvar, sigma):
    slope, intercept = constants
    chi2 = (yvar - (xvar*slope + intercept))**2 / sigma**2
    return chi2


def chi2_sin(constants, xvar, yvar, sigma,amp_bounds=[-1000.0,1000.0],f_bounds=[0.1,1.0e4],offset_bounds=[0.0,3.15]):
    amp, frequency, offset = constants
    try:
        if ((amp_bounds[0] < amp < amp_bounds[1]) & (f_bounds[0] < frequency < f_bounds[1]) & (offset_bounds[0] < offset < offset_bounds[1])):
            chi2 = (yvar - (amp*np.sin(2.0*pi*frequency*xvar + offset)))**2 / sigma**2
        else:
            chi2 = 1.0e100
    except:
        pdb.set_trace()
    return chi2


def chi2_line_sin(constants, xvar, yvar, sigma):
    slope, intercept, amp, frequency, offset = constants
    chi2 = (yvar - (xvar*slope + intercept + amp*np.sin(2.0 * pi * frequency * xvar + offset)))**2 / sigma**2
    return chi2


class RVsample(object):
    """
    """
    def __init__(self):
        rvnames = np.array([])

        database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
        cur = database.cursor()

        #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
        #connection = sqlite.connect(dbfile)
    # Create a cursor object
        #cur = connection.cursor()
        cur.execute('SELECT name,vz_err FROM spectra')
        for row in cur:
            try:
                tmpCheck = np.float(row[1])
                rvnames = np.append(rvnames,str(row[0]))
            except:
                pass

        rvnames = np.unique(rvnames)
        self.stars = np.array([])
        for tmpName in rvnames:
            tmpStar = RVstar(tmpName)
            if ((len(tmpStar.vzerr) > 0.0) & (tmpStar.mag >= 0.0)):
                self.stars = np.append(self.stars,tmpStar)

    def large_accel_sample(self,sigmaval=1.0,future_yr=2017.38,wBartko=True,Rcut=1.8):
        
        out = open('/u/schappell/tables/sig_accelZ_'+str(future_yr)+'.tex','w')
        out.write('\\usepackage{url,epstopdf,lscape,graphicx,subfig,caption,subcaption} \n')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{document} \n')
        out.write('\\clearpage \n')
        out.write('\\begin{landscape} \n')
        out.write('\\begin{deluxetable}{lccccccccc} \n')
        out.write('\\tabletypesize{\\tiny} \n')
        out.write('\\setlength{\\tabcolsep}{1.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        
        out.write('\\tablecaption{Stars with significant change in RV for '+str(future_yr)+'} \n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write("  \\colhead{K'} & \n")
        out.write('  \\colhead{R$_{2D}$} & \n')
        out.write('  \\colhead{Epochs} & \n')
        out.write('  \\colhead{Average} & \n')
        out.write('  \\multicolumn{2}{c}{Most Recent Measurement} & \n')
        out.write('  \\multicolumn{2}{c}{Predicted RV (km/s)} & \n')
        out.write('  \\colhead{F-test} & \\\\\n')
        
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{Error (km/s)} & \n')
        out.write('  \\colhead{Date} & \n')
        out.write('  \\colhead{RV (km/s)} & \n')
        out.write('  \\colhead{Accel} & \n')
        out.write('  \\colhead{Jerk} & \n')
        out.write('  \\colhead{Value} & \n')
        out.write('} \n')
        out.write('\\startdata \n')
        fmt = '%15s  %1s  %5.2f  %1s  %5.2f  %1s  %2d  %1s  %5.2f  %1s  %5.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %5s  %6.2f  %1s  %6.2f  %4s\n'

        for tmpStar in self.stars:
            future_vz,future_vzerr,recent_date,recent_vz,recent_vzerr,future_jerk,future_jverr,passFtest,ftestval=tmpStar.futureRV(future_yr=future_yr,wBartko=wBartko)
            if ((np.abs(future_vz-recent_vz)>(sigmaval*(future_vzerr+recent_vzerr))) & (tmpStar.r2d<Rcut) & (recent_date!=0.0) & (passFtest==0.0)):
                out.write(fmt % (tmpStar.name, '&', tmpStar.mag, '&', tmpStar.r2d, '&', len(tmpStar.vz), '&', np.average(tmpStar.vzerr), '&', recent_date, '&', recent_vz, '$\pm$', recent_vzerr, '&', future_vz, '$\pm$', future_vzerr, '&', future_jerk, '$\pm$', future_jverr, '&', ftestval, '\\\\'))
            elif((np.abs(future_jerk-recent_vz)>(sigmaval*(future_jverr+recent_vzerr))) & (tmpStar.r2d<Rcut) & (recent_date!=0.0) & (passFtest==1.0)):
                out.write(fmt % (tmpStar.name, '&', tmpStar.mag, '&', tmpStar.r2d, '&', len(tmpStar.vz), '&', np.average(tmpStar.vzerr), '&', recent_date, '&', recent_vz, '$\pm$', recent_vzerr, '&', future_vz, '$\pm$', future_vzerr, '&', future_jerk, '$\pm$', future_jverr, '&', ftestval, '\\\\'))
        out.write('\\hline \n')
        out.write('\\hline \n')
        out.write('\\multicolumn{10}{c}{Non-significant Change in RV} & \n')
        out.write('\\\\ \n')
        out.write('\\hline \n')
        for tmpStar in self.stars:
            future_vz,future_vzerr,recent_date,recent_vz,recent_vzerr,future_jerk,future_jverr,passFtest,ftestval=tmpStar.futureRV(future_yr=future_yr,wBartko=wBartko)
            if ((np.abs(future_vz-recent_vz)<=(sigmaval*(future_vzerr+recent_vzerr))) & (tmpStar.r2d<Rcut) & (recent_date!=0.0) & (passFtest==0.0)):
                out.write(fmt % (tmpStar.name, '&', tmpStar.mag, '&', tmpStar.r2d, '&', len(tmpStar.vz), '&', np.average(tmpStar.vzerr), '&', recent_date, '&', recent_vz, '$\pm$', recent_vzerr, '&', future_vz, '$\pm$', future_vzerr, '&', future_jerk, '$\pm$', future_jverr, '&', ftestval, '\\\\'))
            elif((np.abs(future_jerk-recent_vz)<=(sigmaval*(future_jverr+recent_vzerr))) & (tmpStar.r2d<Rcut) & (recent_date!=0.0) & (passFtest==1.0)):
                out.write(fmt % (tmpStar.name, '&', tmpStar.mag, '&', tmpStar.r2d, '&', len(tmpStar.vz), '&', np.average(tmpStar.vzerr), '&', recent_date, '&', recent_vz, '$\pm$', recent_vzerr, '&', future_vz, '$\pm$', future_vzerr, '&', future_jerk, '$\pm$', future_jverr, '&', ftestval, '\\\\'))

        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\clearpage \n')
        out.write('\\end{landscape} \n')
        out.write('\\end{document} \n')
        out.close()
        


                

    def plotBartko(self):
        py.clf()
        py.plot([0,0],[-1000,1000],'k')
        for tmpStar in self.stars:
            if (tmpStar.r2d > 1.0):
                try:
                    usedex = np.argmax(np.fabs(tmpStar.vz - tmpStar.b09_vz) / np.sqrt(tmpStar.vzerr**2 + tmpStar.b09_vzerr**2))
                    use_err = math.sqrt(tmpStar.vzerr[usedex]**2 + tmpStar.b09_vzerr**2)
                    if (str(tmpStar.name)=='S0-6'):
                        tmpStar.appendRV([102.11,95.11,97.37],[0.28,2.00,1.57],[2016.3675,2016.3702,2016.3730])
                    if (str(tmpStar.name)=='S2-17'):
                        tmpStar.appendRV(79.56,2.68,2014.5055)
                        usedex = np.argmax(np.fabs(tmpStar.vz - tmpStar.b09_vz) / np.sqrt(tmpStar.vzerr**2 + tmpStar.b09_vzerr**2))
                        use_err = math.sqrt(tmpStar.vzerr[usedex]**2 + tmpStar.b09_vzerr**2)
                        py.errorbar(tmpStar.vz[usedex]-tmpStar.b09_vz,tmpStar.vz[usedex],yerr=tmpStar.vzerr[usedex],xerr=use_err,fmt='r.')
                    else:
                        py.errorbar(tmpStar.vz[usedex]-tmpStar.b09_vz,tmpStar.vz[usedex],yerr=tmpStar.vzerr[usedex],xerr=use_err,fmt='b.')
                    if ((abs(tmpStar.vz[usedex] - tmpStar.b09_vz)) > (tmpStar.vzerr[usedex])):
                        print str(tmpStar.name)+' inconsistent between us and Bartko'
                except:
                    pass
        py.xlabel('RV - Bartko 2009 RV (km/s)')
        py.ylabel('RV (km/s)')
        py.ylim([-800,400])
        py.savefig('/u/schappell/plots/allRV_vsBartko.png')
        print ''
        py.clf()
        for tmpStar in self.stars:
            if (tmpStar.r2d > 1.0):
                #pdb.set_trace()
                use_vz = np.append(tmpStar.vz,[])
                use_vzerr = np.append(tmpStar.vzerr,[])
                try:
                    use_vz = np.append(tmpStar.vz,tmpStar.b09_vz)
                    use_vzerr = np.append(tmpStar.vzerr,tmpStar.b09_vzerr)

                except:
                    pass

                if (len(use_vz) > 1):
                    weights = 1./use_vzerr
                    weights /= np.sum(weights)
                    
                    weight_mean = np.average(use_vz,weights=weights)
                    red_chi2 = np.sum((weight_mean - use_vz)**2/use_vzerr**2)/(len(weights) - 1.0)
                    tmpStar.weight_mean = weight_mean ; tmpStar.constRV_rchi2 = red_chi2
                    delta_RV = 0.0
                    for i in range(len(use_vz)-1):
                        for j in np.arange(i+1,len(use_vz),1):
                            if (abs(use_vz[i]-use_vz[j])/math.sqrt(use_vzerr[i]**2+use_vzerr[j]**2) > delta_RV):
                                delta_RV = abs(use_vz[i]-use_vz[j])
                                delta_RVerr = math.sqrt(use_vzerr[i]**2+use_vzerr[j]**2)

                    tmpStar.max_deltaRV = delta_RV
                    tmpStar.max_deltaRVerr = delta_RVerr
                    
                    print str(tmpStar.name)+' mag: '+str(tmpStar.mag)+' num obs: '+str(len(use_vz))+' has weight mean: '+str(round(weight_mean,2))+ \
                        ' km/s reduced chi squared: '+str(round(red_chi2,2))+' and field: '+str(tmpStar.field)+' and yng?: '+str(tmpStar.yng)+ \
                        ' and R2d (as): '+str(tmpStar.r2d)+' and delta RV: '+str(delta_RV)+' km/s'
                    py.plot(weight_mean,red_chi2,'bo')
        py.xlabel('Weighted Mean RV (km/s)')
        py.ylabel(r'Reduced $\chi^2$')
        py.savefig('/u/schappell/plots/all_meanRV_chi2.png')

        py.clf()
        for tmpStar in self.stars:
            try:
                py.errorbar(tmpStar.constRV_rchi2,tmpStar.max_deltaRV,yerr=tmpStar.max_deltaRVerr,fmt='b.')
            except:
                pass
        py.xlabel(r'Reduced $\chi^2$')
        py.ylabel('Delta RV (km/s)')
        py.savefig('/u/schappell/plots/all_max_deltaRV_rchi2.png')


        py.clf()
        for tmpStar in self.stars:
            try:
                py.plot(tmpStar.constRV_rchi2,tmpStar.diskProb,'bo')
            except:
                pass
        py.xlabel(r'Reduced $\chi^2$')
        py.ylabel('Disk Probability')
        py.savefig('/u/schappell/plots/all_diskProb_rchi2.png')


    def plotAllRV(self,minPoints=2,wBartko=False):
        for tmpStar in self.stars:
            if ((len(tmpStar.vz) >= minPoints) | hasattr(tmpStar,'b09_vz')):
                tmpStar.plotRV(wBartko=wBartko)



    def binary_cat(self,wBartko=True):
        for tmpStar in self.stars:
            tmpStar.accel(wBartko=wBartko)
        print ''
        print 'Binary candidates:'
        if (wBartko==True):
            print 'Including points from Bartko 2009'
        for tmpStar in self.stars:
            if (tmpStar.binary_cand == True):
                print tmpStar.name


    def plotPositions(self):
        py.clf()
        py.plot([0.0],[0.0],'k+')
        for tmpStar in self.stars:
            if ((tmpStar.yng == 1.0) & (len(tmpStar.vz) >= 2)):
                py.plot(tmpStar.x,tmpStar.y,'bx')
            elif ((tmpStar.yng == 0.0) & (len(tmpStar.vz) >= 2)):
                py.plot(tmpStar.x,tmpStar.y,'rx')
            elif (len(tmpStar.vz) >= 2):
                py.plot(tmpStar.x,tmpStar.y,'gx')
            elif (tmpStar.yng == 1.0):
                py.plot(tmpStar.x,tmpStar.y,'bD')
            elif (tmpStar.yng == 0.0):
                py.plot(tmpStar.x,tmpStar.y,'rD')
            else:
                py.plot(tmpStar.x,tmpStar.y,'gD')
        py.xlabel('RA (arcsec)')
        py.ylabel('Dec (arcsec)')
        py.xlim([15.0,-5.0])
        py.savefig('/u/schappell/plots/RV_x_y.png')


    def makeTable(self):
        out = open('/u/schappell/tables/RVsample.tex','w')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{singlespace} \n')
        out.write('\\begin{deluxetable}{lccccc} \n')
    #out.write('\\rotate \n')
        out.write('\\tabletypesize{\\small} \n')
        out.write('\\setlength{\\tabcolsep}{3.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        out.write('\\begin{document} \n')
        out.write('\\tablecaption{}\n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write('  \\colhead{$Kp$} & \n')
        out.write('  \\colhead{R_{2D}} & \n')
        out.write('  \\colhead{Number} & \n')
        out.write('  \\colhead{Probability} & \n')
        out.write('  \\colhead{Field} & \n')
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(mag)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{Observations} & \n')
        out.write('  \\colhead{Young} & \n')
        out.write('  \\colhead{} & \n')
        out.write('} \n')
        out.write('\\startdata \n')
        fmt = '%15s  %1s  %5.2f  %1s  %6.3f  %1s  %6s  %1s  %6.2f  %1s  %7s  %4s\n'
        for tmpStar in self.stars:
            if (len(tmpStar.vzerr) >= 1):
                out.write(fmt % (tmpStar.name,'&',tmpStar.mag,'&',tmpStar.r2d,'&',len(tmpStar.vzerr),'&',tmpStar.yng,'&',
                                 tmpStar.field, '\\\\'))
        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\end{singlespace} \n')
        out.write('\\end{document} \n')
        out.close()

        out = open('/u/schappell/tables/RVsample_yng1.tex','w')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{singlespace} \n')
        out.write('\\begin{deluxetable}{lccccc} \n')
    #out.write('\\rotate \n')
        out.write('\\tabletypesize{\\small} \n')
        out.write('\\setlength{\\tabcolsep}{3.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        out.write('\\begin{document} \n')
        out.write('\\tablecaption{}\n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write('  \\colhead{$Kp$} & \n')
        out.write('  \\colhead{R_{2D}} & \n')
        out.write('  \\colhead{Number} & \n')
        out.write('  \\colhead{Probability} & \n')
        out.write('  \\colhead{Field} & \n')
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(mag)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{Observations} & \n')
        out.write('  \\colhead{Young} & \n')
        out.write('  \\colhead{} & \n')
        out.write('} \n')
        out.write('\\startdata \n')
        fmt = '%15s  %1s  %5.2f  %1s  %6.3f  %1s  %6s  %1s  %6.2f  %1s  %7s  %4s\n'
        for tmpStar in self.stars:
            if ((len(tmpStar.vzerr) == 1) & (tmpStar.yng == 1.0)):
                out.write(fmt % (tmpStar.name,'&',tmpStar.mag,'&',tmpStar.r2d,'&',len(tmpStar.vzerr),'&',tmpStar.yng,'&',
                                 tmpStar.field, '\\\\'))
        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\end{singlespace} \n')
        out.write('\\end{document} \n')
        out.close()

        out = open('/u/schappell/tables/RVsample_old1.tex','w')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{singlespace} \n')
        out.write('\\begin{deluxetable}{lccccc} \n')
    #out.write('\\rotate \n')
        out.write('\\tabletypesize{\\small} \n')
        out.write('\\setlength{\\tabcolsep}{3.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        out.write('\\begin{document} \n')
        out.write('\\tablecaption{}\n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write('  \\colhead{$Kp$} & \n')
        out.write('  \\colhead{R_{2D}} & \n')
        out.write('  \\colhead{Number} & \n')
        out.write('  \\colhead{Probability} & \n')
        out.write('  \\colhead{Field} & \n')
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(mag)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{Observations} & \n')
        out.write('  \\colhead{Young} & \n')
        out.write('  \\colhead{} & \n')
        out.write('} \n')
        out.write('\\startdata \n')
        fmt = '%15s  %1s  %5.2f  %1s  %6.3f  %1s  %6s  %1s  %6.2f  %1s  %7s  %4s\n'
        for tmpStar in self.stars:
            if ((len(tmpStar.vzerr) == 1) & (tmpStar.yng==0.0)):
                out.write(fmt % (tmpStar.name,'&',tmpStar.mag,'&',tmpStar.r2d,'&',len(tmpStar.vzerr),'&',tmpStar.yng,'&',
                                 tmpStar.field, '\\\\'))
        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\end{singlespace} \n')
        out.write('\\end{document} \n')
        out.close()

        out = open('/u/schappell/tables/RVsample_yng_multi.tex','w')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{singlespace} \n')
        out.write('\\begin{deluxetable}{lccccc} \n')
    #out.write('\\rotate \n')
        out.write('\\tabletypesize{\\small} \n')
        out.write('\\setlength{\\tabcolsep}{3.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        out.write('\\begin{document} \n')
        out.write('\\tablecaption{}\n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write('  \\colhead{$Kp$} & \n')
        out.write('  \\colhead{R_{2D}} & \n')
        out.write('  \\colhead{Number} & \n')
        out.write('  \\colhead{Probability} & \n')
        out.write('  \\colhead{Field} & \n')
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(mag)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{Observations} & \n')
        out.write('  \\colhead{Young} & \n')
        out.write('  \\colhead{} & \n')
        out.write('} \n')
        out.write('\\startdata \n')
        fmt = '%15s  %1s  %5.2f  %1s  %6.3f  %1s  %6s  %1s  %6.2f  %1s  %7s  %4s\n'
        for tmpStar in self.stars:
            if ((len(tmpStar.vzerr) >= 2) & (tmpStar.yng == 1.0)):
                out.write(fmt % (tmpStar.name,'&',tmpStar.mag,'&',tmpStar.r2d,'&',len(tmpStar.vzerr),'&',tmpStar.yng,'&',
                                 tmpStar.field, '\\\\'))
        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\end{singlespace} \n')
        out.write('\\end{document} \n')
        out.close()

        out = open('/u/schappell/tables/RVsample_old_mutli.tex','w')
        out.write('\\documentclass{aastex} \n')
        out.write('\\begin{singlespace} \n')
        out.write('\\begin{deluxetable}{lccccc} \n')
    #out.write('\\rotate \n')
        out.write('\\tabletypesize{\\small} \n')
        out.write('\\setlength{\\tabcolsep}{3.0mm} \n')
        out.write('\\tablewidth{0pt} \n')
        out.write('\\begin{document} \n')
        out.write('\\tablecaption{}\n')
        out.write('\\tablehead{ \n')
        out.write('  \\colhead{Star} & \n')
        out.write('  \\colhead{$Kp$} & \n')
        out.write('  \\colhead{R_{2D}} & \n')
        out.write('  \\colhead{Number} & \n')
        out.write('  \\colhead{Probability} & \n')
        out.write('  \\colhead{Field} & \n')
        out.write('%\n')
        out.write('  \\colhead{} & \n')
        out.write('  \\colhead{(mag)} & \n')
        out.write('  \\colhead{(arcsec)} & \n')
        out.write('  \\colhead{Observations} & \n')
        out.write('  \\colhead{Young} & \n')
        out.write('  \\colhead{} & \n')
        out.write('} \n')
        out.write('\\startdata \n')
        fmt = '%15s  %1s  %5.2f  %1s  %6.3f  %1s  %6s  %1s  %6.2f  %1s  %7s  %4s\n'
        for tmpStar in self.stars:
            if ((len(tmpStar.vzerr) >= 2) & (tmpStar.yng==0.0)):
                out.write(fmt % (tmpStar.name,'&',tmpStar.mag,'&',tmpStar.r2d,'&',len(tmpStar.vzerr),'&',tmpStar.yng,'&',
                                 tmpStar.field, '\\\\'))
        out.write('\\\\\n')
        out.write('\\enddata \n')
        out.write('\\end{deluxetable} \n')
        out.write('\\end{singlespace} \n')
        out.write('\\end{document} \n')
        out.close()


def PhaseCoverage(period,months=False,years=False,duration=1.0):
    #Give period in days or set months or years flag
    #duration given in days

    years_to_days = 365.0
    months_to_days = years_to_days / 12.0

    obs=np.array([])

    if (months==True):
        period *= months_to_days
    elif (years==True):
        period *= years_to_days 
    if ((months==True) & (years==True)):
        print "NOPE! Cannot do that, the period needs to be in either months OR years"
        stop

    database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
    cur = database.cursor()

    #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    #connection = sqlite.connect(dbfile)
    # Create a cursor object
    #cur = connection.cursor()
    cur.execute("SELECT ddate FROM observations WHERE field='GC Central'")
    for row in cur:
        obs=np.append(obs,np.float(row[0]))

    t0 = np.min(obs)
    time_days = (obs - t0) * years_to_days 
    phase = (time_days % period) * 360.0 / period
    #this assumes that velocity due to binary peaks at first date
    delta = duration * 360.0 / period
    degbins = np.arange(0.0,360.0+delta,delta)

    (num,binEnds,other) = py.hist(phase,degbins)

    covered = np.where(num > 0.0)[0]
    print ""
    print 'Total Converage:'
    print float(len(covered))/len(num)
    print ''
    print '1st quad:'
    binEnds = np.delete(binEnds,len(binEnds)-1)
    tmpDex = np.where((num > 0.0) & (binEnds < 90.0))[0]
    print float(len(tmpDex))*4.0 / len(num)
    print ''
    print '2nd quad:'
    tmpDex = np.where((num > 0.0) & (binEnds < 180.0) & (binEnds >= 90.0))[0]
    print float(len(tmpDex))*4.0 / len(num)
    print ''
    print '3rd quad:'
    tmpDex = np.where((num > 0.0) & (binEnds < 270.0) & (binEnds >= 180.0))[0]
    print float(len(tmpDex))*4.0 / len(num)
    print ''
    print '4th quad:'
    tmpDex = np.where((num > 0.0) & (binEnds < 360.0) & (binEnds >= 270.0))[0]
    print float(len(tmpDex))*4.0 / len(num)
    


def cadence(period,delta,numObs=10):
    obs = np.array([i*delta for i in range(numObs)])

    results = np.cos(2.0*pi*obs/period)
    py.clf()
    py.plot(obs,results,'o')
    py.show()



def mr_period_RV(m1=10.0,poom='days',eccent=0.0,incl=90.0,rverr=20.0):
    #m1 is the mass of the primary in solar masses
    #poom, period order of magnitude, days, years, tens, hundreds
    #mass ratio from zero to one
    #eccent is eccentricity
    #incl is inclination in degrees

    mratio = np.arange(0.0,1.0,0.01)

    if (poom=='days'):
        period = np.arange(0.0,10.0,0.1)*sec_in_yr/365.0
    elif (poom=='years'):
        period = np.arange(0.0,10.0,0.01)*sec_in_yr
    elif (poom=='tens'):
        period = np.arange(0.0,100.0,0.1)*sec_in_yr
    elif (poom=='hundreds'):
        period = np.arange(0.0,1000.0,1.0)*sec_in_yr
    else:
        print "NOPE: didn't give valid input for poom"
        sys.exit()

    RV = np.zeros([len(mratio),len(period)])
    for i in range(len(mratio)):
        for j in range(len(period)):
            RV[i,j] = 1e-5*(2.0*pi*G/period[j])**(1.0/3.0)*m1*math.sin(incl*pi/180.0)/(m1+m1*mratio[i])**(2.0/3.0)/math.sqrt(1.0-eccent**2.0)
            #RV[i,j] = 1e-5*(2.0*pi*G/(period[j]*(m1*msun)**2))**(1.0/3.0)*m1*msun*mratio[i]/np.sqrt(1.0-eccent**2)

    py.clf()

    period /= sec_in_yr
    if (poom=='days'):
        period *= 365.0
    levels=np.arange(0.0,1000.0,5.0)
    cont1 = py.contourf(period,mratio,RV,levels)
    py.contour(period,mratio,RV,levels)
    levels_sig=np.arange(rverr,rverr*5.0+1.0,rverr)
    cont2=py.contour(period,mratio,RV,levels_sig,colors='k')
    if (poom=='days'):
        py.xlabel('Period (days)')
    else:
        py.xlabel('Period (years)')
    py.ylabel('Mass Ratio')
    py.title('Mass Primary: '+str(m1))
    cbar1 = py.colorbar(cont1)
    cbar1.add_lines(cont2)
    cbar1.ax.set_ylabel('RV (km/s)')
    py.savefig('/u/schappell/plots/mratio_period_rv_m1'+str(m1)+'_e'+str(eccent)+'_i'+str(incl)+'.png')


def dist_RV(num=100000,m_min=10.0,m_max=40.0,m_slope=-1.7,p_min=0.0,p_max=3.0,
            p_slope=-0.55,period=False,inclination=True,a_min=0.01,a_max=0.2,
            a_slope=-0.54,logPeriod=True,
            binary_frac=0.8,flag='fromAS_everything_NEW',p_peak=5.5,p_meu=5.2,
            giants=False,fromAS=True,eccent=True):

    #period needs to be in log days if doing it in log space
    #period in years if doing lognormal
    #seperation needs to be in AU
    #mass needs to be in solar masses

#    if (q_slope==0.0):
#        q_rand = np.random.rand(num)*0.9 + 0.1
#    else:
#        q_rand = sample_dist(q_min,q_max,q_slope,num=num)
    if(fromAS==True):
        period=False

    if (fromAS==False):
        m_rand = sample_dist(m_min,m_max,m_slope,num=num) * msun #mass in grams
        m_2_rand = sample_dist(m_min,m_max,m_slope,num=num) * msun #mass in grams
        for i in range(len(m_rand)):
            if (m_rand[i] < m_2_rand[i]):
                m_high = m_2_rand[i] * 1.0
                m_low = m_rand[i] * 1.0
                m_rand[i] = m_high
                m_2_rand[i] = m_low
        mtot = m_rand + m_2_rand

    if (period==True):
        if (logPeriod==True):
            print 'Using period distribution: power law in log space'
            p_rand = 10.0**(sample_dist(p_min,p_max,p_slope,num=num)) * 8.64e4 #period in seconds
        else:
            print 'Using period distribution: log normal'
            p_rand = np.random.lognormal(mean=p_peak,sigma=p_meu,size=num) * sec_in_yr #period in seconds
        a_rand = (p_rand**2*G*mtot/(4.0*pi**2))**(1.0/3.0)

    elif (fromAS==True):
        print 'Using seperation distribution: from array in .txt file'
        a1_matrix = np.loadtxt('/u/schappell/Downloads/binary_parameters_final_all.txt')
        a1_array=a1_matrix[:,1]
        p_array = a1_matrix[:,0] * 86400.0 #days to seconds
        a1dex = np.random.randint(len(a1_array),size=num)
        p_rand = p_array[a1dex]
        a_rand = a1_array[a1dex] * cm_in_au
        #a_rand = a1_array[a1dex] * cm_in_au
        m1_array = a1_matrix[:,2]
        m2_array = a1_matrix[:,3]
        r1_array = a1_matrix[:,6] * cm_in_au
        r2_array = a1_matrix[:,7] * cm_in_au
        m1dex = np.random.randint(len(m1_array),size=num)
        m2dex = np.random.randint(len(m2_array),size=num)
        m_rand = m1_array[m1dex] * msun
        m_2_rand = m2_array[m2dex] * msun
        mtot = m_rand + m_2_rand
        r1_rand = r1_array[m1dex]
        r2_rand = r2_array[m2dex]
        #p_rand = np.sqrt(4.0*pi**2*a_rand**3 / (G*mtot))
    else: #if do not use period distribution, use separation distribution
        print 'Using separation distribution: given power law'
        a_rand = sample_dist(a_min,a_max,a_slope,num=num) * cm_in_au
        p_rand = np.sqrt(4.0*pi**2*a_rand**3 / (G*mtot))

    seperation = ((m_rand/(40.0*msun))**0.78 + (m_2_rand/(40.0*msun))**0.78)*18.0*rsun
    if (giants==True):
        seperation = rsun * 200.0

    if (eccent==True):
        e1_matrix = np.loadtxt('/u/schappell/Downloads/binary_parameters_final_all.txt')
        e1_array = e1_matrix[:,5] #eccentricity of the second star, the smaller star, thus highest RV
        e1dex = np.random.randint(len(e1_array),size=num)
        e_rand = e1_array[e1dex]
        
        not_touching = np.where((r1_rand+r2_rand) < a_rand)[0]
        #not_touching = np.where(seperation < (a_rand*(1.0 - e_rand)))[0]
        RVel = (2.0*pi*G/(p_rand*mtot**2))**(1.0/3.0)*m_rand/np.sqrt(1.0-e_rand**2)
        RVel *= 1e-5 #km/s
        e_rand = e_rand[not_touching]
    else:    
        not_touching = np.where(seperation < a_rand)[0]
        RVel = (m_rand**3 * 2.0 * G * pi / (p_rand * mtot**2))**(1.0/3.0) / 1e5 #km/s
    
    RVel = RVel[not_touching]
    p_rand = p_rand[not_touching]
    m_rand = m_rand[not_touching]
    mtot = mtot[not_touching]
    if (inclination==True):
        incl = np.arcsin(np.random.rand(len(not_touching))*2.0 - 1.0)
        RVel *= abs(np.sin(incl))

    py.clf()
#    bins = np.arange(0.0,np.max(RVel+1.0),10.0)
#    py.hist(RVel,bins,normed=True)
    hist,bins = np.histogram(RVel,bins=70)
    widths = np.diff(bins)
    hist = hist* binary_frac / np.sum(hist*widths)
    py.bar(bins[:-1],hist,widths)
    py.xlabel('Radial Velocity (km/s)')
    py.savefig('/u/schappell/plots/RV_dist_'+flag+'.png')



#    bsize = int(np.max(p_rand/8.64))/5.0e5
#    bins = np.arange(0.0,np.max(p_rand/8.64e4)+1.0,0.2)
#    py.hist(p_rand/8.64e4,bins,normed=True)
    py.clf()
    #pdb.set_trace()
    if (giants==True):
        hist,bins=np.histogram(p_rand/sec_in_yr,range=[0.0,4.0],bins=70)
    else:
        hist,bins=np.histogram(p_rand/8.64e4,range=[0.0,1000.0], bins=5000)
    widths = np.diff(bins)
    hist = hist * binary_frac / float(len(p_rand)) / widths
    py.bar(bins[:-1],hist,widths)
    if (giants==True):
        py.xlabel('Period (Years)')
    else:
        #py.xlim([0.0,100.0])
        py.xlabel('Period (Days)')
    py.xlim([0.0,12.0])
    py.savefig('/u/schappell/plots/period_dist_'+flag+'.png')

    #get distribution of expect changes in RV with two observations
    obs_dex = np.random.randint(len(RVel),size=num)
    phase_1 = np.random.rand(num)*2.0*pi
    phase_2 = np.random.rand(num)*2.0*pi
    if (eccent==True):
        e_obs = e_rand[obs_dex]
        delta_RV = RVel[obs_dex]*abs(np.sin(phase_1)*np.sqrt(1.0+e_obs*np.cos(phase_1))-np.sin(phase_2)*np.sqrt(1.0+e_obs*np.cos(phase_2)))
    else:
        delta_RV = RVel[obs_dex]*abs(np.sin(phase_1) - np.sin(phase_2))

    py.clf()
    hist,bins = np.histogram(delta_RV,bins=300)
    widths = np.diff(bins)
    hist = hist* binary_frac / np.sum(hist*widths)
    py.bar(bins[:-1],hist,widths)
    py.xlabel('Change in RV (km/s)')
    py.ylabel('Normalized Frequency')
    py.xlim([0,200])
    py.savefig('/u/schappell/plots/deltaRV_dist_'+flag+'.png')

    #pdb.set_trace()

    py.clf()
    pl_phase = np.arange(0.0,6.29,0.01)
    py.plot(pl_phase,100.0*np.sin(pl_phase)*np.sqrt(1.0+0.9*np.cos(pl_phase)))
    plph_rand = np.random.rand(2)*2.0*pi
    py.errorbar(plph_rand,100.0*np.sin(plph_rand)*np.sqrt(1.0+0.9*np.cos(plph_rand)),yerr=[10.0,10.0],fmt='o')
    py.xlabel('Phase (radian)')
    py.ylabel('Radial Velocity (km/s)')
    py.xlim([0.0,2.0*pi])
    py.savefig('/u/schappell/plots/phase_RV_one.png')

    name,kp_tab,xtab,ytab,r2d_tab,field,RVz,RVerr,ddate_tab,numobs = yng_RVerr()
    p_detect = np.zeros(len(RVerr))
    for i in range(len(RVerr)):
        tmp_err=RVerr[i]
        det_dex = np.where(delta_RV >= 2.0*tmp_err)[0]
        p_detect[i] = binary_frac * np.float(len(det_dex)) / len(delta_RV)

    fieldex = np.where((field=='S') | (field=='E') |(field=='SE') | (name=='S1-12') | (name=='S1-14') |
                       (name=='S1-2') | (name=='S1-21') | (name=='S1-33') | (name=='S1-8'))[0]

    #pdb.set_trace()

    #fieldex = np.where(r2d >= 1.0)[0]

    out = open('/u/schappell/tables/RV_prob.tex','w')
    out.write('\\documentclass{aastex} \n')
    out.write('\\begin{singlespace} \n')
    out.write('\\begin{deluxetable}{lcccccccccc} \n')
    #out.write('\\rotate \n')
    out.write('\\tabletypesize{\\small} \n')
    out.write('\\setlength{\\tabcolsep}{1.0mm} \n')
    out.write('\\tablewidth{0pt} \n')
    out.write('\\begin{document} \n')
    out.write('\\tablecaption{}\n')
    out.write('\\tablehead{ \n')
    out.write('  \\colhead{Star} & \n')
    out.write('  \\colhead{$Kp$} & \n')
    out.write('  \\colhead{$X$} & \n')
    out.write('  \\colhead{$Y$} & \n')
    out.write('  \\colhead{R_{2D}} & \n')
    out.write('  \\colhead{Field} & \n')
    out.write('  \\colhead{RV} & \n')
    out.write('  \\colhead{$\sigma$_{RV}} & \n')
    out.write('  \\colhead{Date} & \n')
    out.write('  \\colhead{Prob} & \n')
    out.write('  \\colhead{Repeat} & \n')
    out.write('%\n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(mag)} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{(arcsec)} & \n')
    out.write('  \\colhead{} & \n')
    out.write('  \\colhead{(km/s)} & \n')
    out.write('  \\colhead{(km/s)} & \n')
    out.write('  \\colhead{(years)} & \n')
    out.write('  \\colhead{Detect} & \n')
    out.write('  \\colhead{Obs} & \n')
    out.write('} \n')
    out.write('\\startdata \n')
    fmt = '%15s  %1s  %5.2f  %1s  %6.3f  %1s  %6.3f  %1s  %6.3f  %1s  %4s  %1s  %5.2f  %1s  %5.2f  %1s  %5.2f  %1s  %2.4f  %1s  %3s  %4s\n'
    for j in fieldex:
        if (name[j]=='S1-21'):
            field[j] = 'C'

        if (field[j] == 'C'):
            tmpYN = 'Yes'
        else:
            tmpYN = 'No'

        out.write(fmt % (name[j],'&',kp_tab[j],'&',xtab[j],'&',ytab[j],'&',r2d_tab[j],'&',field[j],'&',RVz[j],'&',
                         RVerr[j],'&',ddate_tab[j],'&',p_detect[j],'&',tmpYN, '\\\\'))
    out.write('\\\\\n')
    out.write('\\enddata \n')
    out.write('\\end{deluxetable} \n')
    out.write('\\end{singlespace} \n')
    out.write('\\end{document} \n')
    out.close()
    

#    pdb.set_trace()
    return p_rand/8.64e4


def sample_dist(d_min,d_max,d_slope,num=1000):
    a_value = np.array([[d_min**(d_slope+1.0)/(d_slope+1.0), 1.0],[d_max**(d_slope+1.0)/(d_slope+1.0), 1.0]])
    b_value = np.array([0.0,1.0])
    d_const = np.linalg.solve(a_value,b_value)

    d_rand = np.random.rand(num)
    sample = ((d_rand - d_const[1])*(d_slope+1.0)/d_const[0])**(1.0/(d_slope+1.0))

    return sample


def obs_period(minFrames=3):

    database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
    cur = database.cursor()

    #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    #connection = sqlite.connect(dbfile)
    # Create a cursor object
    #cur = connection.cursor()
    cur.execute("SELECT ddate FROM observations WHERE field='GC Central' AND nFrames>%d" %(minFrames-1))
    obs=np.array([])
    for row in cur:
        obs=np.append(obs,np.float(row[0]))

    delta_t=np.array([])
    for i in range(len(obs)):
        for j in range(len(obs)-1-i):
            tmp = abs(obs[i] - obs[j+1+i]) *365.0
            delta_t=np.append(delta_t,tmp)

    bins = np.arange(0.0,np.max(delta_t)+1.0,20.0)
    py.clf()
    py.hist(delta_t,bins,normed=False)
    py.xlabel('Observing Period (Days)')
    py.savefig('/u/schappell/plots/obs_cadence.png')

    return delta_t


def compare_period(minFrames=3,num=100000,m_min=10.0,m_max=40.0,m_slope=-1.7,
                   p_min=0.0,p_max=3.0,p_slope=-0.55,period=False,inclination=True,
                   a_min=0.01,a_max=0.2,a_slope=-0.54,q_min=0.1,q_max=1.0,
                   q_slope=0.0,logPeriod=True,binary_frac=0.8,p_peak=5.5,p_meu=5.2,
                   flag='highM_fromAS',giants=False,fromAS=True,eccent=True):

    obs = obs_period(minFrames=minFrames)
    pred = dist_RV(num=num,m_min=m_min,m_max=m_max,m_slope=m_slope,p_min=p_min,
                   p_max=p_max,p_slope=p_slope,period=period,inclination=inclination,
                   a_min=a_min,a_max=a_max,a_slope=a_slope,q_min=q_min,q_max=q_max,
                   logPeriod=logPeriod,binary_frac=binary_frac,flag=flag,q_slope=q_slope,
                   p_peak=p_peak,p_meu=p_meu,giants=giants,fromAS=fromAS,eccent=eccent)

    py.clf()
    #bins=np.arange(0.0,3.01,0.05)
    hist, bins = np.histogram(np.log10(pred),bins=70)
    widths = np.diff(bins)
    hist = hist * binary_frac / np.sum(hist*widths)
    py.bar(bins[:-1],hist,widths,label='Model')
    #py.hist(np.log10(pred),bins,label='Model',normed=True)
    py.hist(np.log10(obs),bins,label='Observations',normed=True,facecolor='g')
    py.xlabel('Log Period (Days)')
    py.ylim([0.0,(int(np.max(hist*10.0))+2.0)/10.0])
    py.legend()
    py.savefig('/u/schappell/plots/compare_period_'+flag+'.png')
    #pdb.set_trace()


def vel_dis():

    stars=np.loadtxt('table_info.tab',usecols=[0],dtype='string')
    disk=np.loadtxt('table_info.tab',usecols=[14])

    database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
    #cur = database.cursor()

    #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    vx=np.array([])
    vy=np.array([])
    vz=np.array([])
    err=np.array([])
    prob=np.array([])

    for i in range(len(stars)):
        tmpStar = str(stars[i])
        #connection = sqlite.connect(dbfile)
        # Create a cursor object
        cur = database.cursor()
        cur.execute("SELECT vx,vy,vz,vz_err FROM stars WHERE name='%s'" &(tmpStar))
        for row in cur:
            try:
                vz=np.append(vz,np.float(row[2]))
                err=np.append(err,np.float(row[3]))
                vx=np.append(vx,row[0]*asy_to_kms/1000.0)
                vy=np.append(vy,row[1]*asy_to_kms/1000.0)
                prob=np.append(prob,disk[i])
            except:
                continue

    py.clf()
    py.errorbar(np.log(prob+1e-4),vz,yerr=er,fmt='.')
    py.show()


def yng_RVerr():

    database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
    cur = database.cursor()

    #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
    #connection = sqlite.connect(dbfile)
    # Create a cursor object
    #cur = connection.cursor()
    cur.execute("SELECT vz,vz_err,name,x,y,kp,r2d,ddate FROM stars WHERE young='T'")
    RVz = np.array([])
    RV_err = np.array([])
    name = np.array([])
    xtab = np.array([])
    ytab = np.array([])
    kp_tab = np.array([])
    r2d_tab = np.array([])
    ddate_tab = np.array([])
    for row in cur:
        try:
            RVz = np.append(RVz,np.float(row[0]))
            RV_err = np.append(RV_err,np.float(row[1]))
            name = np.append(name,str(row[2]))
            xtab = np.append(xtab,row[3])
            ytab = np.append(ytab,row[4])
            kp_tab = np.append(kp_tab,row[5])
            r2d_tab = np.append(r2d_tab,row[6])
            ddate_tab = np.append(ddate_tab,row[7])
        except:
            continue
    py.clf()
    py.hist(RV_err,bins=40)
    py.xlabel('RV error (km/s)')
    py.savefig('/u/schappell/plots/yng_RVerr_dist.png')

    #pdb.set_trace()

    field = np.chararray(len(name),itemsize=6)
    numobs = np.zeros(len(name))
    for i in range(len(name)):
        tmpName = name[i]
        cur = connection.cursor()
        cur.execute("SELECT field FROM spectra WHERE name='%s'" %(tmpName))
        for row in cur:
            if ((row[0]=='C') | (row[0]=='N') | (row[0]=='S') | (row[0]=='E') | (row[0]=='W') | 
                (row[0]=='NW') | (row[0]=='NE') | (row[0]=='SW') | (row[0]=='SE')):
                field[i] = row[0]

            numobs[i] += 1

    return name,kp_tab,xtab,ytab,r2d_tab,field,RVz,RV_err,ddate_tab,numobs



def equ_width_plots(starName):

    EW = np.array([])
    EW_err = np.array([])
    RV = np.array([])
    RV_err = np.array([])
    date = np.array([])

    #connection to online database
    database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
    cur = database.cursor()

        #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
        #connection = sqlite.connect(dbfile)
    # Create a cursor object
        #cur = connection.cursor()
    cur.execute("SELECT ddate,vlsr,vz_err,eq_width,eq_width_err FROM spectra WHERE name='%s'"%starName)
    for row in cur:
        try:
            if ((row[1] != None) & (row[2] != None) & (row[3] != None) & (row[4] != None)):
                date = np.append(date,row[0])
                RV = np.append(RV,row[1])
                RV_err = np.append(RV_err,row[2]) 
                EW = np.append(EW,row[3])
                EW_err = np.append(EW_err,row[4])
        except:
            continue

    py.clf()
    py.errorbar(date,EW,yerr=EW_err,fmt='.')
    py.xlabel('Date (year)')
    py.ylabel('Equivalent Width')
    py.savefig('/u/schappell/plots/'+str(starName)+'_EW_date.png')
    py.clf()
    py.errorbar(RV,EW,xerr=RV_err,yerr=EW_err,fmt='.')
    py.xlabel('RV (km/s)')
    py.ylabel('Equivalent Width')
    py.savefig('/u/schappell/plots/'+str(starName)+'_EW_RV.png')
    py.clf()




def SNR_plot(starName):

    SNR = np.array([])
    date = np.array([])

    #connection to online database
    database = mysqldb.connect(host="galaxy1.astro.ucla.edu",user="dbread",passwd="t36fCEtw",db="gcg")
    cur = database.cursor()

        #dbfile = '/u/ghezgroup/data/gc/database/stars.sqlite'
    # Create a connection to the database file
        #connection = sqlite.connect(dbfile)
    # Create a cursor object
        #cur = connection.cursor()
    cur.execute("SELECT ddate,SNR FROM spectra WHERE name='%s'"%starName)
    for row in cur:
        try:
            if ((row[1] != None) & (row[1] != 0)):
                date = np.append(date,row[0])
                SNR = np.append(SNR,row[1])
        except:
            continue

    if (starName == 'S0-2'):
        date = np.append(date,2016.3675)
        SNR = np.append(SNR,46.220861)

    py.clf()
    py.plot(date,SNR,'o')
    py.xlabel('Date (year)')
    py.ylabel('SNR')
    py.savefig('/u/schappell/plots/'+str(starName)+'_spec_SNR.png')
    py.clf()



def plotSpectra(file,name,file_type='.fits',telluric=False,tell_file='/u/schappell/090505_tell.dat'):
    if (file_type=='.fits'):
        star = pyfits.open(file)
        flux = star[0].data
        hdr = star[0].header

        min_wavelength = hdr['CRVAL1'] / 10.0
        delta_wavelength = hdr['CDELT1'] / 10.0

        wavelength = np.array(range(len(flux)))*delta_wavelength + min_wavelength #in nm

    elif (file_type=='.dat'):
        star = np.loadtxt(file)
        flux = star[:,1]
        wavelength = star[:,0] * 1000.0
        min_wavelength = np.min(wavelength)

    else:
        print 'Wrong input for file_type'

    if (telluric==True):
        tell = np.loadtxt(tell_file)
        flux /= tell

    max_wavelength = np.max(wavelength)
    py.clf()
    py.plot(wavelength, flux)
    py.xlabel('Wavelength (nm)')
    py.ylabel('Flux')
    py.title(name)
    py.xlim([min_wavelength,max_wavelength])
    py.savefig('/u/schappell/plots/spectra_'+str(name)+'.png')
    pdb.set_trace()
