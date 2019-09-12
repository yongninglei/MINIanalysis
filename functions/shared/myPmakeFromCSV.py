import os
import glob as glob
import numpy as np
import subprocess as sp
import pandas as pd
import re

join = os.path.join

# tempdir = "/Users/gari/gDrive/BCBL/PROYECTOS/MINI/ANALYSIS/temp"
tempdir = "/Users/gari/gDrive/BCBL/PROYECTOS/MINI/ANALYSIS/temp"
os.chdir(tempdir)

# Read the csv
A = pd.read_csv('con2dist.csv')
for ii in range(len(A)):
    print ii
    vtx1 = A.vtx1[ii]
    vtx2 = A.vtx2[ii]
    if np.isnan(vtx1) or np.isnan(vtx2):
        print 'uno de los dos vertex es nan'
    else:
        cmd2 = str('mris_pmake --subject fsaverage --hemi lh --surface0 smoothwm '+
                   '--curv0 sulc --curv1 sulc --mpmOverlay euclidean ' + 
                   '--mpmProg pathFind --mpmArgs startVertex:'+ 
                   str(int(vtx1))+',endVertex:'+str(int(vtx2)))
    try:
        output = sp.check_output(cmd2, shell =True)
        A.surfDist[ii] = float(re.split('[\\[\\]]', output)[3].split()[0])
    except sp.CalledProcessError as e:
        print e.output
