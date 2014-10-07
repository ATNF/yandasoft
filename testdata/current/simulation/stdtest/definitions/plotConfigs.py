#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

fin=open('ConfigurationData/ASKAP-SEIC-0005_Antenna_Configuration.csv','r')
lines=fin.readlines()
x=[]
y=[]
for i in range(4,40):
    cols=lines[i].split(',')
    x.append(cols[6])
    y.append(cols[7])

x=np.array(x,dtype=float)
y=np.array(y,dtype=float)
pads=np.array(np.arange(36)+1,dtype='|S2')

xav=x.mean()
yav=y.mean()

fig = plt.figure(1, figsize=(11.7,11.7), dpi=72)

offset=20
textsize='xx-small'
#offset=0
#textsize='medium'

plt.plot(x-xav,y-yav,'.')
plt.xlim(-3500,3500)
plt.ylim(-3500,3500)
#plt.text(x-xav+20,y-yav+20,pads)
for i in range(1,37):
#    plt.text(x[i-1]-xav,y[i-1]-yav,'%s'%i)
    plt.text(x[i-1]-xav+offset,y[i-1]-yav+offset,'%s'%i,size=textsize)

# plot BETA antennas in red
for i in [0,2,5,7,8,14]:
    plt.plot(x[i]-xav,y[i]-yav,'r.')

# plt ASKAP-12 antennas in green
for i in [1,3,4,13,9,11,23,26,29,12,15,27]:
    plt.plot(x[i]-xav,y[i]-yav,'g.')

theta=np.arange(0,1000)*2*np.pi/1000.

for radius in [500.,750.,1500.,3000.]:
    plt.plot(radius*np.cos(theta),radius*np.sin(theta))

plt.savefig('config.png')

