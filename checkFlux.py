import pyfits
import numpy as np
import math
import matplotlib.path as mplPath

def compareFlux(dark_sub_file,reduced_file,edges_file,pixFromMax=6):
    darkSub_file = pyfits.open(dark_sub_file)
    reduced_file = pyfits.open(reduced_file)

    darkSub = darkSub_file[0].data
    reduced = reduced_file[0].data

    edges = np.loadtxt(edges_file)

    #x1_edges = edges[:,0]
    #y1_edges = edges[:,1]
    #dex2 = [1,2,3,0]
    #x2_edges = x1_edges[dex2]
    #y2_edges = y1_edges[dex2]
    darkTot = 0.0

    quadAp = mplPath.Path(edges)

    for i in range(2048):
        for j in range(2048):
            #x1_tmp = x1_edges - i
            #x2_tmp = x2_edges - i
            #y1_tmp = y1_edges - j
            #y2_tmp = y2_edges - j
            #dp = (x1_tmp*x2_tmp) + (y1_tmp*y2_tmp)
            #cp = (x1_tmp*y2_tmp) - (y1_tmp*x2_tmp)
            #theta = np.arctan(cp/dp)

            #if (abs(np.sum(theta)) > math.pi):
            if (darkSub[i,j] > 0):
                darkTot += darkSub[i,j] * quadAp.contains_point((i,j))

    sum_lambda = np.sum(reduced,axis=2)
    sum_y = np.sum(sum_lambda,axis=1)
    maxdex = np.argmax(sum_y)

    lowDex = maxdex - pixFromMax
    highDex = maxdex - pixFromMax
    if (lowDex < 0):
        lowDex = 0
    if (highDex > 18):
        highDex = 18

    redTot = np.sum(sum_y[lowDex:highDex])

    return darkTot, redTot
