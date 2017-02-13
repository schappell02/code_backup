import scipy
import pyfits
import asciidata, os, sys
import numpy as np
import pylab as py

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
