# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 11:00:31 2020

@author: Bryan
"""

import numpy as np
import matplotlib.pyplot as plt
import math as mth
import seaborn as sns

#第一題 第一小題 計算馬路最大應力============================================

#==========參數設定===========
l=3.5 #馬路寬度 單位:m
w=15*10**(-2) #馬路厚度 單位:m
r=(1/np.sin(10*np.pi/180))*(3.5/2)
r_=(r+w)*(20*np.pi/180) 
lmax=r_ #馬路最大可膨脹到的長度 單位:m
G=10**(8) #揚氏係數 單位:N/(m^2)

#==========計算最大應力===========
Δl=lmax-l
strainmax=(lmax-l)/l #馬路最大的應變 
stressmax=G*strainmax #馬路所能承受的最大應力 單位:N/(m^2) 
print('馬路所能承受的最大應力',stressmax,'N/m^2')
print('馬路最大的應變',strainmax)
print('馬路最大的形變量',Δl,'m')

#第一題 第二小題 計算多少溫度差的情況下，馬路會產生破裂========================

#===========參數設定===========
α=6.0*10**(-4) #熱膨脹係數 單位:K^(-1)
l=3.5 #馬路寬度 單位:m
w=15*10**(-2) #馬路厚度 單位:m
r=(1/np.sin(10*np.pi/180))*(3.5/2)
r_=(r+w)*(20*np.pi/180) 
lmax=r_ #馬路最大可膨脹到的長度 單位:m
G=10**(8) #揚氏係數 單位:N/(m^2)

#==========計算多少溫度差的情況下，馬路會產生破裂==========
ΔT=stressmax/(G*α)
print('最大溫度差',ΔT,'K')

#===========參數設定===========
α=6.0*10**(-4) #熱膨脹係數 單位:K^(-1)
l=3.5 #馬路寬度 單位:m
w=15*10**(-2) #馬路厚度 單位:m
z=0.0001
lmax=r_ #馬路最大可膨脹到的長度 單位:m
G=10**(8) #揚氏係數 單位:N/(m^2)
ΔT=stressmax/(G*α)
n=100 #將設定好的θ分成100段
posx1=[]
posy1=[]
posx2=[]
posy2=[]
zz=[]
for g in range(20):
    θ=z*np.pi/180 
    dθ=θ/n #分成100段時 每一段θ的大小
    x=(np.pi-θ)/2 #設定新的能代表θ的符號x
    for i in range(n+1):
        r=(1/abs(np.sin((z/2)*np.pi/180)))*(3.5/2)
        r_=(r+w)*(z*np.pi/180) 
        rtotal=r+w #對最上層馬路來說彎曲時的半徑 單位:m
        posx1.append((rtotal)*np.cos(x))
        posy1.append((rtotal)*abs(np.sin(x))-(r))
        posx2.append(r*np.cos(x))
        posy2.append(-w)
        x=x+dθ    
    plt.title("simulate road's thermal expandsion")
    plt.xlabel('relative x position')
    plt.ylabel('relative y position')
    
    plt.ylim(min(posy2)*2,max(posy1)*2)
    
    plt.plot(posx2,posy2,linewidth=10,color='#054E9F',label='bottom layer of the road')
    plt.plot(posx1,posy1,linewidth=10,color='red',label='top layer of the road')
    plt.legend()
    plt.show()
    z+=1
    zz.append(z)
    posx1.clear()
    posy1.clear()
    posx2.clear()
    posy2.clear()




#==========模擬馬路最上層熱膨脹==========
ΔT=stressmax/(G*α)
rtotal=r+w #對最上層馬路來說彎曲時的半徑 單位:m
θ=20*np.pi/180 
n=100 #將設定好的θ分成100段
dθ=θ/n #分成100段時 每一段θ的大小
x=(np.pi-θ)/2 #設定新的能代表θ的符號x
posx1=[]
posy1=[]
posx2=[]
posy2=[]
for i in range(n+1):
    posx1.append((rtotal)*np.cos(x))
    posy1.append((rtotal)*np.sin(x)-(r))
    posx2.append(r*np.cos(x))
    posy2.append(-w)
    x=x+dθ
plt.title("simulate road's thermal expandsion")
plt.xlabel('relative x position')
plt.ylabel('relative y position')
plt.ylim(min(posy2)*2,max(posy1)*2)
plt.plot(posx2,posy2,linewidth=10,color='#054E9F',label='bottom layer of the road')
plt.plot(posx1[0:45],posy1[0:45],linewidth=10,color='red',label='top layer of the road')
plt.plot(posx1[56:100],posy1[56:100],linewidth=10,color='red')
plt.legend()
plt.show()

#第一題 第三小題 計算水泥中間需要留多少的縫隙，才能避免溫差很大時，柏油路面產生破裂==================

#===========參數設定===========
α=6.0*10**(-4) #熱膨脹係數 單位:K^(-1)
l=3.5 #馬路寬度 單位:m
w=15*10**(-2) #馬路厚度 單位:m
r=(1/np.sin(10*np.pi/180))*(3.5/2)
r_=(r+w)*(20*np.pi/180) 
lmax=r_ #馬路最大可膨脹到的長度 單位:m
G=10**(8) #揚氏係數 單位:N/(m^2)
k=0.17*1000 #柏油的導熱係數 單位:w/(m*K)
N=100 #表示將柏油分成多少層
dx=w/N #一層柏油的厚度 單位:m
dh=0.001 #馬路的長度 單位:m
dt=0.001 #前一層的柏油傳給下一層柏油所需的時間 單位:s
Ac=l*dh #馬路的截面積 單位:m^2
Th=40+273 #最上層柏油的溫度 單位:K
Tc=-1+273 #最下層柏油的溫度 單位:K
Tnorm=25+273 #室溫 單位:K
ρ=2243 #柏油密度 單位:kg/m^3
V=l*dh*dx #一層柏油的體積 單位:m^3
m=ρ*V #一層柏油的質量 單位:kg
s=1670 #柏油的比熱 單位:J/(kg*K)

#==========計算水泥中間需要留多少的縫隙，才能避免柏油路面產生破裂===========

N=100 #表示將柏油分成多少層
dx=w/N
T=[Th]+[Tnorm for i in range(N-2)]+[Tc] #存放各層柏油的溫度
newT=[T[i] for i in range(N)] #用以儲存新計算的溫度
newT1=[Tc for i in range(N)] #用以儲存新計算的溫度
L=[l for j in range(N)]
newL=[L[j] for j in range(N)]
B=k/(s*ρ)
for M in range(30001):
    if M==30000:
        plt.title("simulate each roadlayer's temperature")
        plt.xlabel('relative x position')
        plt.ylabel('temperature')
        plt.plot(T)
        plt.show()

    for i in range(N):
        if i!=0 and i!=N-1: #更新除了熱源與冷源
            ΔT=-dt*B*((T[i]-T[i+1])+(T[i]-T[i-1]))/(dx**2) #計算delta T
            newT[i]=T[i]+ΔT #計算每個點的新溫度          
    T=[newT[i] for i in range(N)] #新溫度為下次計算時的原溫度
    if M==30000:
        ΔT_=[T[i]-newT1[i] for i in range(N)]
        for j in range(N):
            ΔL=α*l*ΔT_[j]
            if ΔL>=Δl:
                print((100-j)*w/100,'m','馬路此深度位置會破裂')
            else:
                print((100-j)*w/100,'m','馬路此深度位置不會破裂')
            newL[j]=L[j]+ΔL
        
        posx3=[]
        posy3=[]
        for j in range(N):
            for g in range(10):
                if j%10==0:
                    x=(newL[j]-L[j])/10*g
                    y=(100-j)*dx
                    posx3.append(x)
                    posy3.append(y)
        
        if j==10:
            for g in range(10):
                plt.plot(posx3,posy3)
        
            
        a=np.array([posx3,posy3])
        b=[]
        c=[]
        
        
    
        for j in range(N):
            b.append(a[0][j])
            c.append(a[1][j])
        plt.title("simulate tile's thermal expandsion in each layer")
        plt.xlabel('the length of the expansion')
        plt.ylabel('each layer')
        plt.plot(b[0:9],c[0:9])
        plt.plot(b[10:19],c[10:19])
        plt.plot(b[20:29],c[20:29])
        plt.plot(b[30:39],c[30:39])
        plt.plot(b[40:49],c[40:49])
        plt.plot(b[50:59],c[50:59])
        plt.plot(b[60:69],c[60:69])
        plt.plot(b[70:79],c[70:79])
        plt.plot(b[80:89],c[80:89])
        plt.plot(b[90:99],c[90:99])
        plt.show()

#第二題 第一小題 計算磁磚間隙要預留多少==================

#===========參數設定===========
α=9.0*10**(-6) #線性熱膨脹係數 單位:K^(-1)
β=3*α #非線性熱膨脹係數 單位:K^(-1)
l=0.5 #磁磚寬度 單位:m
w=1*10**(-2) #磁磚厚度 單位:m
G=0.5*10**(6) #揚氏係數 單位:N/(m^2)
k=1.5 #磁磚的導熱係數 單位:w/(m*K)
N=100 #表示將磁磚分成多少層
dx=w/N #一層磁磚的厚度 單位:m
dt=0.001 #前一層的磁磚傳給下一層磁磚所需的時間 單位:s
Ac=l**2 #馬路的截面積 單位:m^2
Th=100+273 #最上層磁磚的溫度 單位:K
Tc=0+273 #最下層磁磚的溫度 單位:K
Tnorm=25+273 #室溫 單位:K
V=(l**2)*w #磁磚的體積 單位:m^3
V_=(l**2)*dx #一層磁磚的體積 單位:m^3
s=1.05*10**(3) #磁磚的比熱 單位:J/kg
ρ=2.2*10**(3) #磁磚的密度 單位:kg/m^3

#==========計算水泥中間需要留多少的縫隙，才能避免柏油路面產生破裂===========

N=100 #表示將磁磚分成多少層
dx=w/N
T=[Th]+[Tnorm for i in range(N-2)]+[Tc] #存放各層柏油的溫度
newT=[T[i] for i in range(N)] #用以儲存新計算的溫度
newT1=[Tc for i in range(N)] #用以儲存新計算的溫度
V=[l**2 for j in range(N)]
newV=[V[j] for j in range(N)]
B=k/(s*ρ)

Δv=V_*0.002 #設定體積膨脹率最大只能為0.2%
for M in range(28001):
    if M==28000:
        plt.title("simulate each tilelayer's temperature")
        plt.xlabel('relative x position')
        plt.ylabel('temperature')
        plt.plot(T)
        plt.show()

    for i in range(N):
        if i!=0 and i!=N-1: #更新除了熱源與冷源
            ΔT=-dt*B*((T[i]-T[i+1])+(T[i]-T[i-1]))/(dx**2) #計算delta T
            newT[i]=T[i]+ΔT #計算每個點的新溫度          
    T=[newT[i] for i in range(N)] #新溫度為下次計算時的原溫度

    if M==28000:
        ΔT_=[T[i]-newT1[i] for i in range(N)]
        for j in range(N):
            ΔV=β*(V_)*ΔT_[j]
            if ΔV>=Δv:
                print((100-j)*w/100,'此深度位置磁磚會破裂',ΔV,'此深度位置磁磚的體膨脹量')
            else:
                print((100-j)*w/100,'此深度位置磁磚不會破裂',ΔV,'此深度位置磁磚的體膨脹量')
            newV[j]=V[j]+ΔV
            
rbposx1=[]
lbposx1=[]
mbposx1=[]
mbposx1_=[]
bposx1=[]
bposx2=[]
bposx2_=[]
    
rbposy1=[]
lbposy1=[]
mbposy1=[]
mbposy1_=[]
bposy1=[]
bposy2_=[]
bposy2=[]

for j in range(N-1):
    length=mth.sqrt(abs(newV[j+1]-newV[j]))
    
    if j==1:
        bx1=0
        by1=0
        for g in range(1,4):
            bx1+=(length/2)
            by1+=(length/2)
            mbposx1.append(bx1)
            mbposy1.append(by1)
 
        xm=mbposx1[1]
        ym=mbposy1[1]
        mbposx1_.append(xm)
        mbposy1_.append(ym)
        
for j in range(N-1):

    length=(mth.sqrt(abs(newV[j+1]-newV[j])))
    
    xm_=xm
    ym_=ym
    bx1_=(length/2)
    by1_=(length/2)
    rbposx1.append(xm_+bx1_)
    lbposx1.append(xm_-bx1_)
    rbposy1.append(ym_+by1_)
    lbposy1.append(ym_-by1_)
    bposx1=lbposx1+mbposx1_+rbposx1
    bposy1=lbposy1+mbposy1_+rbposy1
    bposy2.append(bposy1[0])
    bposy2.append(bposy1[0])
    bposy2.append(bposy1[0])
    bposy2_.append(bposy1[2])
    bposy2_.append(bposy1[2])
    bposy2_.append(bposy1[2])
    bposx2.append(bposx1[0])
    bposx2.append(bposx1[0])
    bposx2.append(bposx1[0])
    bposx2_.append(bposx1[2])
    bposx2_.append(bposx1[2])
    bposx2_.append(bposx1[2])
    if j%20==0:
        #colormap=plt.cm.gist_ncar
        #colors=[colormap(i) for i in range(6)]
        plt.title("simulate tile's thermal expandsion in each layer")
        plt.xlabel('relative x position')
        plt.ylabel('relative y position')
        plt.axis('equal')
        plt.plot(bposx1[0:3],bposy2[0:3],color='blue')
        plt.plot(bposx1[0:3],bposy2_[0:3],color='blue')
        plt.plot(bposx2[0:3],bposy1[0:3],color='blue')
        plt.plot(bposx2_[0:3],bposy1[0:3],color='blue')
        
    
    bposx1.clear()
    bposx2.clear()
    bposx2_.clear()
    bposy1.clear()
    bposy2_.clear()
    bposy2.clear()
    rbposx1.clear()
    lbposx1.clear()
    rbposy1.clear()
    lbposy1.clear()


    

        
            
        
        
        
        



  
    

    
    

    






    
    
    
    
    
    
    
    

    
    
    
    







