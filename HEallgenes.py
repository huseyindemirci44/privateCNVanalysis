# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:31:51 2020

@author: huseyin.demirci
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Pyfhel import Pyfhel 
import time

HE = Pyfhel();
#HE.contextGen(p = 65537, m = 8192, flagBatching = True);
HE.contextGen(p = 65537, m = 32768, flagBatching = True);


m0= time.perf_counter();
HE.keyGen();
m1= time.perf_counter();
print("Key generation Time: ", m1-m0);





healthygenes = pd.read_csv('allgenes_healthy_csv.csv', sep = ',', names=[ 'geneid', 'genename', 'coverage' ]);
healthygenes.head()


patientgenes = pd.read_csv('allgenes_patient_csv.csv', sep = ',', names=[ 'geneid', 'genename', 'coverage' ]);
patientgenes.head()




r1 =  (patientgenes['coverage']<5) &  (healthygenes['coverage']<5)
r2 = np.logical_not(r1);



h1= healthygenes[r2]['coverage'];
p1= patientgenes[r2]['coverage'];



healthycoverage = np.log2(h1+1) ;
hcov_int = healthycoverage.astype(int);


patientcoverage = np.log2(p1+1) ;
pcov_int = patientcoverage.astype(int);


diff =  pcov_int - hcov_int ;

div = (p1+0.1)/(h1+0.1);





def sma(arr, n):
   weights = np.ones(n) / n

   return np.convolve(weights, arr)[n-1:-n+1]

h2 = sma(h1,3);
p2 = sma(p1,3);

div2 = (p1+0.1)/(h1+0.1);

div3 = sma(div2,2);


#i1 = 3;
#i2 = 7;
#
#c1 = HE.encryptInt(i1);
#c2 = HE.encryptInt(i2);
#
#csum = c1+c2;
#csub = c1-c2;
#cmul = c1*c2;
#
#dsum = HE.decryptInt(csum);
#dsub = HE.decryptInt(csub);
#dmul = HE.decryptInt(cmul);


th1 = 1.60 ;
th2 = 1.5 ;
th3 = 0.8 ;
th4 = 0.1 ;
l = len(h1);


ff =  np.zeros(l);
ff = ff.astype(np.int64);


b5 = np.array(round(h1*th1));
b5 = b5.astype(np.int64);
b5 = np.array(b5)

t0= time.perf_counter();

ca5 = HE.encryptArray(np.array(p1));
cb5 = HE.encryptArray(b5);
c5 = ca5-cb5;

d5 = c5.decrypt();
d5 = d5[0:l];
d5 = np.array(d5);


ff[(d5>0)] = 2 ;
#print(ff)



b52 = np.array(round(h1*th2));
b52 = b52.astype(np.int64);
b52 = np.array(b52)
cb52 = HE.encryptArray(b52);
c52 = ca5-cb52;
d52 = c52.decrypt();
d52 = d52[0:l];
d52 = np.array(d52);


ff[(d5<=0) & (d52>0) ] = 1 ;
#print(ff)


b53 = np.array(round(h1*th3));
b53 = b53.astype(np.int64);
b53 = np.array(b53)
cb53 = HE.encryptArray(b53);
c53 = ca5-cb53;
d53 = c53.decrypt();
d53 = d53[0:l];
d53 = np.array(d53);


ff[(d53<0) ] = -1 ;
#print(ff)

b54 = np.array(round(h1*th4));
b54 = b54.astype(np.int64);
b54 = np.array(b54)
cb54 = HE.encryptArray(b54);
c54 = ca5-cb54;
d54 = c54.decrypt();
d54 = d54[0:l];
d54 = np.array(d54);


ff[(d54<0) ] = -2 ;
#print(ff)


t1= time.perf_counter();
print("Encrypion and Devryption  Time: ", t1-t0);

regionbegin = 0 ; 

regionend = 350;



#plt.plot(h1[regionbegin:regionend], 'o', color='red', label='Healthy')
#plt.plot(p1[regionbegin:regionend], 'x', color='blue', label='Melanoma')

#plt.plot(hcov_int, 'o', color='red', label='Healthy')
#plt.plot(pcov_int, 'x', color='blue', label='Melanoma')

#plt.plot(diff[regionbegin:regionend], 'o', color='red', label='Difference')

#y = div3[regionbegin:regionend] ;
#y5 = sma(y,3) ; 
#x = np.arange(len(y));
#
#plt.scatter(x, y, color='blue', label='Ratio')
#
#plt.legend()
#
#plt.title('Melanoma  Tumor vs Healthy Coverage')
#
#
#plt.show()
#
#
#plt.scatter( x[regionbegin:regionend] , (d5>0)[regionbegin:regionend], color='blue')
#
##plt.legend()
#
#plt.title('Melanoma  vs Healthy ')
#
#plt.show()


f5 = sma(ff,3) ; 

x = np.arange(len(ff));

plt.scatter( x[0:l] , ff[0:l], color='blue')
plt.title('Copy Number Estimations ')
plt.show()

plt.scatter( x[0:l-5] , f5[0:l-5], color='blue')
plt.title('Smoothed Copy Number Estimations ')
plt.show()






