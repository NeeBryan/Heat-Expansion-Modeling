# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 04:13:44 2020

@author: Bryan
"""



import numpy as np 
import matplotlib.pyplot as plt
from vpython import *

Dair = 1.225 #(kg/m^3)
Dwater = 1000 #(kg/m^3)

g =  9.8 #(m/s^2)
Wout =32 #(Hz)
Cd = 0.47 
t=0
v=0
a=0
y=0
H0=0.008
r=1.5*(10**(-3))
viscosity=0.001
n=1#mole
R=0.08205
T=298.15
h=-0.002
m=Dair*(4/3)*(np.pi)*(r**3)
k=np.abs(Dwater*g*h*(4/3)*(np.pi)*(0.0015**3))


ball=sphere(pos=vector(0,-0.7238,0),radius=3)


ballorigin=sphere(pos=vector(4,-0.7238,0),radius=3)


b=box(pos=vec(0,0,0),size=vec(50,0.1,50))

scene = canvas(y=-0.002*10000000,width=10,height=20)


tt=[]
yy=[]
rr=[]
vv=[]
ff=[]
hh=[]
aa=[]
dd=[]
vo=[]
kk=[]
fo=[]
ww=[]
bb=[]
pp=[]
ll=[]
l=-0.002
ll. append(l)
aa=[]
a=1
aa. append(a)
bb=[]
b=1
bb. append(b)
while t<=0.001:
    

    dt=0.00001
    
    V=np.abs((k/(Dwater*g*h)))
    r=((3*V)/(4*np.pi))**(1/3)
    W=(-m)*g
    B=Dwater*V*g
    if v==0:
        D=0
    elif v>0:
        '''
        for i in range(75):
            a=1+i*0.01
            aa. append(a)
            D=-(1/2)*Dwater*(v)*Cd* (np.pi) *((aa[i]*r)**2)
            '''
        D=-(1/2)*Dwater*(v)*Cd* (np.pi) *((1.74*r)**2)
    else:
        '''
        for i in range(11):
            b=1+i*0.01
            bb. append(b)
            D=-(1/2)*Dwater*(v)*Cd* (np.pi) *((bb[i]*r)**2)
        '''   
        D=-(1/2)*Dwater*(v)*Cd* (np.pi) *((1.1*r)**2)
    Fout=m*H0*((Wout)**2)*np.cos((Wout*t+(np.pi)))
    Ft=W+B+D+Fout
    a=Ft/m
    v+=a*dt
    h+=v*dt
    y=h
    
    tt. append(t)
    yy. append(y)
    rr. append(r)
    vv. append(v)
    ff. append(Ft)
    hh. append(h)
    aa. append(a)
    dd. append(D)
    vo. append(V)
    kk. append(k)
    fo. append(Fout)
    ww. append(W)
    bb. append(B)
    t+=dt
'''
    ball.pos=vec(0,g,0)
    rate(1)
    t+=dt   
  

for i in range(99):
    p=10000000*(yy[i+1]-yy[i])
    pp. append(p)
    l=(ll[i]+pp[i])
    ll. append(l)



    ball.pos=vec(0,l,0)
    rate(3)


ballorigin=sphere(pos=vector(5,-5,0),radius=1)
'''
   
    

''' 
print(k)
print(m)

plt.plot(tt,yy)

plt.plot(tt,fo)

plt.plot(tt,ff)
'''

plt.xlabel("time (s)")
plt.ylabel("position (m)")
plt.plot(tt,yy)


