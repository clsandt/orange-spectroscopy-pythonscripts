# Authors: 
# Christophe Sandt - christophe.sandt@synchrotron-soleil.fr
# Marko Toplak - @markotoplak
# Ferenc Borondics - @borondics

import scipy.optimize
import numpy as np
from copy import deepcopy
from orangecontrib.spectroscopy.data import getx


def nearest_index(array, n): return np.argmin(np.abs(array - n))

 
out_data = deepcopy(in_datas[0])
  
# Borns for subtraction
wn1=1216
wn2=1318
wn3=1271
  
wns_data = getx(in_datas[0])
  
id1 = nearest_index(wns_data, wn1) #peak boundary 1
id2 = nearest_index(wns_data, wn2) #peak boundary 2
id3 = nearest_index(wns_data, wn3) # position of the peak center for initial guess
 
#Baseline correction
sub_peak=deepcopy(in_datas[0].X[:,id1:id2])
lastid=sub_peak.shape[1]-1
 
#reference baseline correction
x1=[wns_data[id1],wns_data[id2]]
y1=[in_datas[1].X[0,id1],in_datas[1].X[0,id2]]
p1=np.polyfit(x1,y1,1)
ldb1=np.polyval(p1,wns_data[id1:id2])
reference=in_datas[1].X[0,id1:id2]-np.transpose(ldb1)
 
#peak baseline correction
x0=[wns_data[id1],wns_data[id2]]
 
for i in range(in_datas[0].X.shape[0]):
    y0=[sub_peak[i,0],sub_peak[i,lastid]]
    p0=np.polyfit(x0,y0,1)
    ldb0=np.polyval(p0,wns_data[id1:id2])
    sub_peak[i,:]=sub_peak[i,:]-ldb0
 
for i in range(in_datas[0].X.shape[0]):
    # function to minimize
    def f(x): return np.sum((sub_peak[i,:]-x*reference)**2)

    # initialize subtraction factor
    ifx=in_datas[0].X[i,id3]/in_datas[1].X[0,id3]

    # find optimum
    sf=scipy.optimize.fmin(f,ifx)
    print("subtraction factor",sf)

    if sf>0:
        out_data.X[i,:] = in_datas[0].X[i,:]-sf*in_datas[1].X
