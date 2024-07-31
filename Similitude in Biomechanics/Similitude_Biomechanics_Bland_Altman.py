# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:11:24 2020

@author: Raul Ochoa Cabrero

Bland-Altman Statistical Analysis for Scaled Biomechanical Experiments

Copyright: Copyright (c) 2020 Raul Ochoa Cabrero
"""

import numpy as np
import matplotlib.pyplot as plt

# Repeatability A (Sawbones)

A=np.array([[0.075,0.073,0.076,0.077,0.07,0.073,0.07,0.073,0.073,0.074],[0.081,0.082,0.083,0.079,0.081,0.08,0.079,0.083,0.079,0.08],[0.084,0.085,0.082,0.08,0.085,0.079,0.086,0.081,0.08,0.084],[0.089,0.089,0.09,0.085,0.09,0.082,0.088,0.09,0.084,0.086]])

n=A.shape[1]
m=A.shape[0]

X=np.zeros((int(n*(n-1)/2)*m,1))
Y=np.zeros((int(n*(n-1)/2)*m,1))
Z=np.zeros((int(n*(n-1)/2)*m,2))
Sizes=np.zeros((int(n*(n-1)/2)*m,1))

for k in range(0,m):
    for i in range(0,n-1):
        for j in range(0,n-1-i):
            X[int(n*(n-1)/2)*k+sum(range(n-i,n))+j]=A[k,i]-A[k,i+(j+1)]


for k in range(0,m):
    for i in range(0,n-1):
        for j in range(0,n-1-i):
            Y[int(n*(n-1)/2)*k+sum(range(n-i,n))+j]=(A[k,i]+A[k,i+(j+1)])/2 

for k in range(0,int(n*(n-1)/2)*m):
    Z[k,:]=[X[k],Y[k]]

for k in range(0,int(n*(n-1)/2)*m):
    for j in range(0,k):
        if np.all(Z[k,:]==Z[j,:]):
            Sizes[k]=Sizes[k]+1
            
Sizes=Sizes/max(Sizes)+1
area = (6 * Sizes)**2

meanA=np.mean(X)
stdA=np.std(X)
rcoeffA=1.96*np.std(X)

upperlimitA=np.mean(X)+1.96*np.std(X)
lowerlimitA=np.mean(X)-1.96*np.std(X)

# Repeatability B (3D printed)

B=np.array([[0.07,0.072,0.075,0.069,0.07,0.069,0.074,0.07,0.073,0.073],[0.078,0.079,0.08,0.077,0.08,0.082,0.08,0.08,0.082,0.08],[0.085,0.08,0.085,0.079,0.084,0.082,0.081,0.081,0.084,0.083],[0.086,0.088,0.084,0.092,0.088,0.083,0.085,0.083,0.089,0.088]])

B=np.array([[0.067,0.065,0.073,0.065,0.067,0.073,0.07,0.068,0.066,0.068],[0.08,0.078,0.077,0.077,0.077,0.082,0.078,0.083,0.077,0.076],[0.085,0.08,0.08,0.076,0.081,0.077,0.081,0.083,0.084,0.081],[0.087,0.087,0.083,0.086,0.081,0.082,0.086,0.089,0.086,0.09]])

B=np.array([[0.078,0.078,0.076,0.074,0.079,0.078,0.078,0.078,0.073,0.073],[0.081,0.079,0.077,0.075,0.08,0.078,0.078,0.075,0.08,0.075],[0.084,0.085,0.081,0.078,0.080,0.082,0.087,0.085,0.082,0.078],[0.085,0.090,0.086,0.083,0.092,0.084,0.086,0.088,0.089,0.085]])

n=B.shape[1]
m=B.shape[0]

W=np.zeros((int(n*(n-1)/2)*m,1))
V=np.zeros((int(n*(n-1)/2)*m,1))
Z=np.zeros((int(n*(n-1)/2)*m,2))
Sizes=np.zeros((int(n*(n-1)/2)*m,1))

for k in range(0,m):
    for i in range(0,n-1):
        for j in range(0,n-1-i):
            W[int(n*(n-1)/2)*k+sum(range(n-i,n))+j]=B[k,i]-B[k,i+(j+1)]


for k in range(0,m):
    for i in range(0,n-1):
        for j in range(0,n-1-i):
            V[int(n*(n-1)/2)*k+sum(range(n-i,n))+j]=(B[k,i]+B[k,i+(j+1)])/2 

for k in range(0,int(n*(n-1)/2)*m):
    Z[k,:]=[W[k],V[k]]

for k in range(0,int(n*(n-1)/2)*m):
    for j in range(0,k):
        if np.all(Z[k,:]==Z[j,:]):
            Sizes[k]=Sizes[k]+1
            
Sizes=Sizes/max(Sizes)+1
area = (6 * Sizes)**2

meanB=np.mean(W)
stdB=np.std(W)
rcoeffB=1.96*np.std(W)

upperlimitB=np.mean(W)+1.96*np.std(W)
lowerlimitB=np.mean(W)-1.96*np.std(W)

# Plots

plt.scatter(Y,X,s=area,alpha=0.5,color='orangered')
plt.xlim(0.067, 0.092)
plt.ylim(-0.011, 0.011)
plt.axhline(linewidth=1, color='black')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("Repeatability Composite Bone", size=20)
plt.savefig('RepeatabilityAssesmentSawbones.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()

plt.scatter(V,W,s=area,alpha=0.5,color='orangered')
plt.xlim(0.07, 0.095)
plt.ylim(-0.011, 0.011)
plt.axhline(linewidth=1, color='black')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("Repeatability PLA/TPU 90A", size=20)
plt.savefig('RepeatabilityAssesment90A.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()

# Bland Altman FE models

C=np.array([[0.0741],[0.0776],[0.081],[0.0845]])
D=np.array([[0.0724],[0.0749],[0.0793],[0.0841]])
E=np.array([[0.0747],[0.078],[0.0812],[0.0845]])

n=A.shape[1]
m=A.shape[0]

X=np.zeros((n*m,1))
Y=np.zeros((n*m,1))
Z=np.zeros((n*m,2))
Sizes=np.zeros((n*m,1))

for k in range(0,m):
    for i in range(0,n):
        X[k*n+i]=A[k,i]-C[k,0]                                                                                                                                 

for k in range(0,m):
    for i in range(0,n):
        Y[k*n+i]=(A[k,i]+C[k,0])/2

for k in range(0,m*n):
    Z[k,:]=[X[k],Y[k]]

for k in range(0,m*n):
    for j in range(0,k):
        if np.all(Z[k,:]==Z[j,:]):
            Sizes[k]=Sizes[k]+1
            
Sizes=Sizes/max(Sizes)+1
area = (6 * Sizes)**2

meanAC=np.mean(X)
stdAC=np.std(X)
upperlimitAC=np.mean(X)+1.96*np.std(X)
lowerlimitAC=np.mean(X)-1.96*np.std(X)

plt.scatter(Y,X,s=area,alpha=0.5,color='orangered')
plt.xlim(0.069, 0.091)
plt.ylim(-0.0065, 0.006)
plt.axhline(linewidth=1, color='black')
plt.axhline(y=upperlimitAC,linewidth=1, color='r',linestyle='--')
plt.axhline(y=lowerlimitAC,linewidth=1, color='b',linestyle='--')
plt.axhline(y=meanAC,linewidth=1, color='black',linestyle='--')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("Composite Bone vs. FE Model", size=20)
plt.savefig('BlandAltmanSawbonesvsFE.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()

X=np.zeros((n*m,1))
Y=np.zeros((n*m,1))
Z=np.zeros((n*m,2))
Sizes=np.zeros((n*m,1))

for k in range(0,m):
    for i in range(0,n):
        X[k*n+i]=B[k,i]-D[k,0]

for k in range(0,m):
    for i in range(0,n):
        Y[k*n+i]=(B[k,i]+D[k,0])/2

for k in range(0,m*n):
    Z[k,:]=[X[k],Y[k]]

for k in range(0,m*n):
    for j in range(0,k):
        if np.all(Z[k,:]==Z[j,:]):
            Sizes[k]=Sizes[k]+1
            
Sizes=Sizes/max(Sizes)+1
area = (6 * Sizes)**2

meanBD=np.mean(X)
stdBD=np.std(X)
upperlimitBD=np.mean(X)+1.96*np.std(X)
lowerlimitBD=np.mean(X)-1.96*np.std(X)

plt.scatter(Y,X,s=area,alpha=0.5,color='orangered')
plt.xlim(0.063, 0.09)
plt.ylim(-0.009, 0.011)
plt.axhline(linewidth=1, color='black')
plt.axhline(y=upperlimitBD,linewidth=1, color='r',linestyle='--')
plt.axhline(y=lowerlimitBD,linewidth=1, color='b',linestyle='--')
plt.axhline(y=meanBD,linewidth=1, color='black',linestyle='--')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("PLA/TPU 95A vs. FE Model", size=20)
plt.savefig('BlandAltman95AvsFE.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()

X=np.zeros((n*m,1))
Y=np.zeros((n*m,1))
Z=np.zeros((n*m,2))
Sizes=np.zeros((n*m,1))

for k in range(0,m):
    for i in range(0,n):
        X[k*n+i]=B[k,i]-E[k,0]

for k in range(0,m):
    for i in range(0,n):
        Y[k*n+i]=(B[k,i]+E[k,0])/2

for k in range(0,m*n):
    Z[k,:]=[X[k],Y[k]]

for k in range(0,m*n):
    for j in range(0,k):
        if np.all(Z[k,:]==Z[j,:]):
            Sizes[k]=Sizes[k]+1
            
Sizes=Sizes/max(Sizes)+1
area = (6 * Sizes)**2

meanBE=np.mean(X)
stdBE=np.std(X)
upperlimitBE=np.mean(X)+1.96*np.std(X)
lowerlimitBE=np.mean(X)-1.96*np.std(X)

plt.scatter(Y,X,s=area,alpha=0.5,color='orangered')
plt.xlim(0.071, 0.091)
plt.ylim(-0.0065, 0.0095)
plt.axhline(linewidth=1, color='black')
plt.axhline(y=upperlimitBE,linewidth=1, color='r',linestyle='--')
plt.axhline(y=lowerlimitBE,linewidth=1, color='b',linestyle='--')
plt.axhline(y=meanBE,linewidth=1, color='black',linestyle='--')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("PLA/TPU 90A vs. FE Model", size=20)
plt.savefig('BlandAltman90AvsFE.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()


X=np.zeros((n*m,1))
Y=np.zeros((n*m,1))
Z=np.zeros((n*m,2))
Sizes=np.zeros((n*m,1))

for k in range(0,m):
    for i in range(0,n):
        X[k*n+i]=B[k,i]-C[k,0]

for k in range(0,m):
    for i in range(0,n):
        Y[k*n+i]=(B[k,i]+C[k,0])/2

for k in range(0,m*n):
    Z[k,:]=[X[k],Y[k]]

for k in range(0,m*n):
    for j in range(0,k):
        if np.all(Z[k,:]==Z[j,:]):
            Sizes[k]=Sizes[k]+1
            
Sizes=Sizes/max(Sizes)+1
area = (6 * Sizes)**2

meanBD=np.mean(X)
stdBD=np.std(X)
upperlimitBD=np.mean(X)+1.96*np.std(X)
lowerlimitBD=np.mean(X)-1.96*np.std(X)

plt.scatter(Y,X,s=area,alpha=0.5,color='orangered')
plt.xlim(0.063, 0.09)
plt.ylim(-0.008, 0.008)
plt.axhline(linewidth=1, color='black')
plt.axhline(y=0.0056,linewidth=1, color='r',linestyle='--')
plt.axhline(y=-0.0048,linewidth=1, color='b',linestyle='--')
plt.axhline(y=-0.0004,linewidth=1, color='black',linestyle='--')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("PLA/TPU 93A vs. FE Model", size=20)
plt.savefig('BlandAltman93AvsFE.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()

# Bland Altman

n=A.shape[1]
m=A.shape[0]

X=np.zeros((n*n*m,1))
Y=np.zeros((n*n*m,1))

for k in range(0,m):
    for i in range(0,n):
        for j in range(0,n):
            X[n*n*k+n*i+j]=A[k,i]-B[k,j]


for k in range(0,m):
    for i in range(0,n):
        for j in range(0,n):
            Y[n*n*k+n*i+j]=(A[k,i]+B[k,j])/2

meanAB=np.mean(X)
stdAB=np.std(X)
upperlimitAB=np.mean(X)+1.96*np.std(X)
lowerlimitAB=np.mean(X)-1.96*np.std(X)

plt.scatter(Y,X,s=area,alpha=0.5,color='orangered')
plt.xlim(0.071, 0.091)
plt.ylim(-0.0065, 0.0095)
plt.axhline(linewidth=1, color='black')
plt.axhline(y=upperlimitAB,linewidth=1, color='r',linestyle='--')
plt.axhline(y=lowerlimitAB,linewidth=1, color='b',linestyle='--')
plt.axhline(y=meanAB,linewidth=1, color='black',linestyle='--')
plt.xlabel("Average", size=14)
plt.ylabel("Differences", size=14)
plt.title("Composite Bone vs. PLA/TPU 93A", size=20)
plt.savefig('BlandAltman90AvsFE.jpg',dpi=300,quality=95,bbox_inches = "tight")
plt.show()



