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

home='/u/schappell/'

pi = math.pi

G = 6.6726e-8
msun = 1.99e33
rsun = 6.955e10
sec_in_yr = 3.1557e7
cm_in_au = 1.496e13
cm_in_pc = 3.086e18
km_in_pc = 3.086e13
au_in_pc = 206265.0


def FixHeader(fileToFix,fileWithHead,set,frame,mjd_obs=0,filter='Kbb'):
#  Edits incomplete header of given fits file, fileToFix, with the header of
#  fileWithHead, saves resulting fits file as fileToFix+'_FH.fits'. Only
#  set, frame, and filter keyword is changed between the header of fileWithHead
#  and that of resulting output fits file.
#
#  Inputs:
#         fileToFix - Name of file w/ incomplete header w/o .fits extension
#                     ex: 's160515_a001002'
#         fileWithHead - Name of file w/ complete header and as much in common
#                        with fileToFix as possible, only set number, frame number,
#                        filter, and mjd-obs changed between these two fits files
#         set - set number of fileToFix
#         frame - frame number of fileToFix
#         mjd_obs - MJD-OBS keyword for fileToFix, only set if it is different
#                   from that of fileWithHead
#         filter - filter of fileToFix

#  Example: 'FixHeader('s160515_a009005','s160515_a009007',9,5)'


    command = 'cp '+fileToFix+'.fits '+fileToFix+'_FH.fits'
    os.system(command)

    header = pyfits.getheader(fileWithHead+'.fits')
    data = pyfits.getdata(fileToFix+'_FH.fits')

    header.update('DATAFILE',fileToFix+"_FH.fits")
    header.update('FRAMENUM',str(frame))
    header.update('SETNUM',str(set))
    header.update('SFILTER',filter)
    if (mjd_obs != 0):
        header.update('MJD-OBS',mjd_obs)

    pyfits.update(fileToFix+'_FH.fits',data,header,ext=0)
