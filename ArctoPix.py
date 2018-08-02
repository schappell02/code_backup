import numpy as np
import pdb
from scipy.optimize import fsolve
from scipy.optimize import minimize


def skyToRefPix(align='test_schoedel/14_01_15/',points='efit/S0-38_97may_AOadderr1.7_constholoerr7.orig.points',printout=True):

    pointsTab=np.loadtxt('/g/ghez/align/'+align+points)

    years = pointsTab[:,0]
    x = pointsTab[:,1]
    y = pointsTab[:,2]

    xpix = np.zeros(len(years))
    ypix = np.zeros(len(years))

    trans_date = np.loadtxt('/g/ghez/align/'+align+'align/align_d_rms_1000_abs_t.date')
    trans = np.loadtxt('/g/ghez/align/'+align+'align/align_d_rms_1000_abs_t.trans')

    for i in range(len(years)):
        j = np.where(years[i] == trans_date)[0]

        const_array = np.array([[float(trans[j,5]),float(trans[j,7])],[float(trans[j,13]),float(trans[j,11])]])
        zero_array = np.array([x[i]-float(trans[j,3]),float(y[i]-trans[j,9])])

        #pdb.set_trace()
        solu = np.linalg.solve(const_array,zero_array)

        if (printout==True):
            print years[i]
            print solu
            print ''

        xpix[i] = solu[0]
        ypix[i] = solu[1]

    return years, xpix, ypix



def RefToPix(align='test_schoedel/14_01_15/',points='efit/S0-38_97may_AOadderr1.7_constholoerr7.orig.points'):

    years, xref, yref = skyToRefPix(align=align,points=points,printout=False)
    date = np.loadtxt('/g/ghez/align/'+align+'align/align_d_rms_1000.date')
    trans = np.loadtxt('/g/ghez/align/'+align+'align/align_d_rms_1000.trans')

    #x = np.zeros(len(years))
    #y = np.zeros(len(years))

    for i in range(len(years)):
        j = np.where(years[i] == date)[0]
        #pixels = fsolve(pix2pix,[1000.0,1000.0],args=([xref[i],yref[i]],trans[j,:].T))
        #print years[i]
        #print pixels
        #print ''

        const_array = np.array([[float(trans[j,5]),float(trans[j,7])],[float(trans[j,19]),float(trans[j,17])]])
        zero_array = np.array([xref[i]-float(trans[j,3]),float(yref[i]-trans[j,15])])

        #pdb.set_trace()
        solu = np.linalg.solve(const_array,zero_array)

        print years[i]
        print solu
        print ''



        
def pix2pix(pixels,*needed):
    refpix, trans = needed
    return1 = trans[3]+trans[5]*pixels[0]+trans[7]*pixels[1]+trans[9]*pixels[0]**2
    return1 += trans[11]*pixels[0]*pixels[1]+trans[13]*pixels[1]**2
    return2 = trans[15]+trans[17]*pixels[1]+trans[19]*pixels[0]+trans[21]*pixels[1]**2
    return2 += trans[23]*pixels[1]*pixels[0]+trans[25]*pixels[0]**2

    output = np.zeros(2)
    output[0] = abs(refpix[0] - return1)
    output[1] = abs(refpix[1] - return2)
    return output
