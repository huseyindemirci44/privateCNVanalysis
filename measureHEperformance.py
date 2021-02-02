# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:31:51 2020

@author: huseyin.demirci
to measure the time performance of HE experiment on genomics data.
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Pyfhel import Pyfhel 

import time
from scipy import stats

import seaborn as sns


HE = Pyfhel();
#HE.contextGen(p = 65537, m = 4096*8, flagBatching = True);
#HE.contextGen(p = 65537, m = 2048, flagBatching = True);
HE.contextGen(p = 65537, m = 1024, flagBatching = True);

m0= time.perf_counter();
HE.keyGen();
m1= time.perf_counter();
print("Key generation Time: ", m1-m0);





healthygenes = pd.read_csv('healthygenes_chomosome9.csv', sep = ',', names=['linenumber',  'geneid', 'genename', 'coverage' ]);
healthygenes.head()


patientgenes = pd.read_csv('patientgenes_chomosome9.csv', sep = ',', names=['linenumber',  'geneid', 'genename', 'coverage' ]);
patientgenes.head()




r1 =  (patientgenes['coverage']<5) &  (healthygenes['coverage']<5)
r2 = np.logical_not(r1);



h1= healthygenes[r2]['coverage'];
p1= patientgenes[r2]['coverage'];



th1 = 1.70 ;
th2 = 1.50 ;
th3 = 0.45 ;
th4 = 0.1 ;
l = len(h1);
l = 1000; ## comment this line if you want to test over all genes
# size of the array to be tested



ff =  np.zeros(l);
ff = ff.astype(np.int64);


b5 = np.array(round(h1*th1));
b5 = b5.astype(np.int64);
b5 = np.array(b5)



## to test 100 genes

#l = 100;

plain1 = p1[0:l];
b51 = b5[0:l];
t0= time.perf_counter();

ca5 = HE.encryptArray(np.array(plain1));
cb5 = HE.encryptArray(b51);
c5 = ca5-cb5;




t1 = time.perf_counter();

 
#print("Encrypion Time: ", t1) # CPU seconds elapsed (floating point)

d5 = c5.decrypt();
d5 = d5[0:l];
d5 = np.array(d5);
t2= time.perf_counter();

#print("Encrypion Time: ", t1-t0, "Decrypion Time: ", t1-t2) # CPU seconds elapsed (floating point)

ff[(d5>0)] = 2 ;
#print(ff)



b52 = np.array(round(h1*th2));
b52 = b52.astype(np.int64);
b52 = np.array(b52)
b52 = b52[0:l];


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
b53 = b53[0:l];

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
b54 = b54[0:l];

cb54 = HE.encryptArray(b54);
c54 = ca5-cb54;
d54 = c54.decrypt();
d54 = d54[0:l];
d54 = np.array(d54);

t3= time.perf_counter();

print("Encrypion Time: ", round(t1-t0,5), "Decrypion Time: ", round(t2-t1,5) ,"All operations finished Time: ", round(t3-t0,5 )); # CPU seconds elapsed (floating point)




ff[(d54<0) ] = -2 ;
#print(ff)


df =pd.read_csv('HEexperimentresults.txt')
print(df.head())
#figure1 = 
#sns.relplot(x= 'm' ,y = 'total_time', data = df, markers = True, kind = 'line')
df1 = df.melt('m', var_name='cols',  value_name='Execution Time (in seconds)')

g = sns.catplot(x="m", y="Execution Time (in seconds)", hue='cols', data=df1, kind='point')
#g.savefig('timeperf.eps')





