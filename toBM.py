from gcwork import objects
from gcwork import starset
from gcwork import util
from gcwork import orbits
from gcwork import young
from pysqlite2 import dbapi2 as sqlite
import scipy
import pyfits # can import whole libraries
from scipy import stats # can import just a sub-library
from gcwork import starTables
import pickle
import nmpfit_sy
import asciidata, os, sys, pickle
import nmpfit_sy2 as nmpfit_sy
from pylab import *
import numpy as np
import pylab as py #can give another name to this library
import math
import histNofill
#import matplotlib.axes3d as p3
import pdb #this is loading the debugger
import scipy.optimize as opter

# import libraries to use built in functions

#<-----Use this symbol to comment out a line, used for commenting code

"""
use above to comment out a whole section
of text and use below to close that section
"""

#Warning, this code is not meant to work, so don't compile it

#A .py file may have multiple functions in it

def test_fuction(input, option = True):
    does things
    return other_things

def another_function():

    py.clf() #cleans the figue from anything that mayhave been left from previous plots
    py.figure(figsize=(12,5))
    py.subplots_adjust(wspace=0.25, hspace=0.25, left=0.08, right=0.95, top=0.95)
    #this will be multiple plots to one image
    py.subplot(1,3,1)
    py.semilogy(np.hypot(x,y), np.abs(vx), 'r.') #semilogy is some special plotting function, don't remember what it does
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('|X Velocity| (mas/yr)') #label axis
    py.subplot(1,3,2)
    py.semilogy(np.hypot(x,y), np.abs(vy), 'b.')
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('|Y Velocity| (mas/yr)')
    py.subplot(1,3,3)
    py.semilogy(np.hypot(x,y), np.hypot(vx,vy), 'k.')
    py.xlabel('Projected Radius (arcsec)')
    py.ylabel('Proper Motion (mas/yr)')
    py.savefig(outdir + 'pm_vs_r2d.png') #this saves it to the directory you give it
    # outdir would be previously defined

    
#an example of defining a function and all the variables you can set
"""
In python window:
import toBM
output = toBM.loadPop()

if you wanted to change the alnDir
output = toBM.loadPop(alnDir='different/')
"""
def loadPop(alnDir='13_08_21/',
            align = 'align/align_d_rms_1000_abs_t', root_tmp='/g/ghez/align/',
            poly='polyfit_c/fit',points='points_c/', starlist='yngstars'):

    does things

    return s

def something():

    py.clf()
    py.subplots_adjust(hspace=0.2, left=0.15, right=0.85,
                       top=0.9, bottom=0.1)
    py.subplot(2, 1, 1)
    py.plot(x,y, 'k.') #plots a scatter plot k calls the color (black) and . tells the function to use dots to plot the scatter
    #using o in .'s place would make them larger dots
    py.axis([-20, 20, 0, 45])
    #adjust axis range (if you need to)
    py.xlabel('Radial Acc. Sig. (sigma)')
    py.ylabel('Number of Epochs Detected')

    py.subplot(2, 1, 2)
    py.plot(at[idx] / ate[idx], cnt[idx], 'k.')
    py.axis([-20, 20, 0, 45])
    py.xlabel('Tangent. Acc. Sig. (sigma)')
    py.ylabel('Number of Epochs Detected')

    py.savefig(outdir + 'polyfit_hist_accel_nepochs%s.png' % tag)

    py.figure(6) #can label an individual figure
    py.clf()
    py.plot(ar*1e3,at*1e3,'k.')
    leg1=py.errorbar(ar[idnp]*1e3,at[idnp]*1e3,xerr=are[idnp]*3e3,yerr=ate[idnp]*3e3,fmt='.',label='Non Phys. Accel.') #can plot w/ errorbars
    leg2=py.errorbar(ar[idsig]*1e3,at[idsig]*1e3,xerr=are[idsig]*3e3,yerr=ate[idsig]*3e3,fmt='.',label='Sig. Accel.')
    py.plot([-10,10],[0,0],'k')
    py.plot([0,0],[-10,10],'k')
    py.axis([-1.2,4,-1.2,1.2])
    py.xlabel('Radial Acceleration (mas/yr/yr)')
    py.ylabel('Tangent Acceleration (mas/yr/yr)')
    py.legend() #compiles legend
    py.savefig(outdir + 'tan_rad.png')

