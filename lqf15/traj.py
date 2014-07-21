##!/usr/bin/python

import numpy as np
import pylab as plt

#with open("traj.dat") as f:
#    data = f.read()
#
#    data = data.split('\n')
#
#    x = [row.split(' ')[0] for row in data]
#    y = [row.split(' ')[1] for row in data]
#
#    fig = plt.figure()
#
#    ax1 = fig.add_subplot(111)
#
#    ax1.set_title("Plot title...")    
#    ax1.set_xlabel('your x label..')
#    ax1.set_ylabel('your y label...')
#
#    ax1.plot(x,y, c='r', label='the data')
#
#    leg = ax1.legend()
#fig = plt.figure()
plt.subplot(311)
#plt.ylim(0,4)
data = np.genfromtxt(fname='q.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,data.shape[1]):
    plt.plot(data[:,0],data[:,x])

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
plt.xlabel('time')
plt.ylabel('$x_i$')
plt.title('traj')

plt.subplot(312)
data = np.genfromtxt(fname='c.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,10):
    plt.plot(data[:,0],data[:,x])
plt.xlabel('time')
plt.ylim(-0.2,2)

plt.subplot(313)
data = np.genfromtxt(fname='norm')
plt.plot(data[:,0],data[:,1],'r-',linewidth=2)
plt.ylim(0,2)

#plt.subplot(2,2,4)
#data = np.genfromtxt(fname='wf.dat')
##data1 = np.genfromtxt(fname='wf0.dat')
#data0 = np.genfromtxt('../spo_1d/t500')
#plt.plot(data[:,0],data[:,1],'r--',linewidth=2)
#plt.plot(data0[:,0],data0[:,1],'k-',linewidth=2)
##plt.plot(data1[:,0],data1[:,1],'k-.',linewidth=2)
##plt.title('t=100')
#
##plt.figure(1) 
##plt.plot(x,y1,'-')
##plt.plot(x,y2,'g-')
#plt.xlim(0.8,2.1)
#plt.xlabel('x')
#plt.ylabel('$\psi^*\psi$')

plt.savefig('wf.pdf')
plt.show() 



