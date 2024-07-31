# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 09:42:00 2022

@author: Raul Ochoa

Extended Finite Similitude in Finance: Size Effects

Copyright: Copyright (c) 2022 Raul Ochoa Cabrero
"""

import math
import sympy
from sympy import symbols
import numpy as np
import csv
import random

np.set_printoptions(precision=4)

#function to convert from cm to market proportions
def cmtoprop(x,n):
    r=pow(10,(1/3.9*x-n))
    return(r)

#function to convert from cm to returns
def cmtoret(x):
    r=.004/2.85*x
    return(r)

FigureData=np.zeros((25,2))
Means=np.zeros((4,2))

with open('FigureData.csv', 'r') as x:
    d = list(csv.reader(x,delimiter=","))
FigureData=np.array(d)
FigureData[0,0]=0.000041
FigureData=FigureData.astype(float)

for i in range(0,4):
    Means[i,0]=np.mean(FigureData[i*5:i*5+10,0])
    Means[i,1]=np.mean(FigureData[i*5:i*5+10,1])
    


Means=np.array([[0.00004,0.0034],[0.00012,0.0007],[0.00027,-0.0005],[0.00067,-0.001],[0.0052,0.00001]])
Betas=np.array([[1,2.7/6.7,1.2/6.7,0.4/6.7],[1,2.7/52,1.2/52,0.4/52],[1,6.7/52,1.2/52,0.4/52],[1,6.7/52,2.7/52,0.4/52],[1,6.7/52,2.7/52,1.2/52]])

epsilon=np.zeros((5,4))
for i in range(0,4):
    epsilon[0,i]=Means[3-i,1]*pow(Betas[0,i],-1)

epsilon[1,0]=Means[4,1]*pow(Betas[1,0],-1)
epsilon[1,1]=Means[2,1]*pow(Betas[1,1],-1)
epsilon[1,2]=Means[1,1]*pow(Betas[1,2],-1)
epsilon[1,3]=Means[0,1]*pow(Betas[1,3],-1)

epsilon[2,0]=Means[4,1]*pow(Betas[2,0],-1)
epsilon[2,1]=Means[3,1]*pow(Betas[2,1],-1)
epsilon[2,2]=Means[1,1]*pow(Betas[2,2],-1)
epsilon[2,3]=Means[0,1]*pow(Betas[2,3],-1)

epsilon[3,0]=Means[4,1]*pow(Betas[3,0],-1)
epsilon[3,1]=Means[3,1]*pow(Betas[3,1],-1)
epsilon[3,2]=Means[2,1]*pow(Betas[3,2],-1)
epsilon[3,3]=Means[0,1]*pow(Betas[3,3],-1)

for i in range(0,4):
    epsilon[4,i]=Means[3-i+1,1]*pow(Betas[4,i],-1)

R1, R2, R3 =symbols('R1 R2 R3')

eq1=sum((epsilon[:,1]-epsilon[:,0])*(epsilon[:,1]-epsilon[:,2]))+R1*sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,1]-epsilon[:,2]))+R2*sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,1]-epsilon[:,2]))-R2*R3*sum((epsilon[:,2]-epsilon[:,3])*(epsilon[:,1]-epsilon[:,2]))

eq2=sum((epsilon[:,1]-epsilon[:,0])*(epsilon[:,2]-epsilon[:,3]))+R1*sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,2]-epsilon[:,3]))+R2*sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,2]-epsilon[:,3]))-R2*R3*sum((epsilon[:,2]-epsilon[:,3])*(epsilon[:,2]-epsilon[:,3]))


sympy.solve(eq2.subs(R1,sympy.solve(eq1,R1)[0]),R2)[0]
#The above expression is the value of R2

eq3=sympy.solve(eq2.subs(R1,sympy.solve(eq1,R1)[0]),R2)[0]

eq4=sympy.solve(eq1,R1)[0]

#Hessian matrix
def Hessian(P3):
    H=np.zeros((3,3))
    H[0,1]=sum((epsilon[:,1]-epsilon[:,2])*((epsilon[:,1]-epsilon[:,2])-P3*(epsilon[:,2]-epsilon[:,3])))
    H[0,0]=sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,1]-epsilon[:,2]))
    H[0,2]=-eq3.subs(R3,P3)*sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,2]-epsilon[:,3]))
    H[2,0]=sum((epsilon[:,1]-epsilon[:,2])*(epsilon[:,2]-epsilon[:,3]))
    H[2,1]=sum((epsilon[:,2]-epsilon[:,3])*((epsilon[:,1]-epsilon[:,2])-P3*(epsilon[:,2]-epsilon[:,3])))
    H[2,2]=-eq3.subs(R3,P3)*sum((epsilon[:,2]-epsilon[:,3])*(epsilon[:,2]-epsilon[:,3]))
    H[1,0]=H[0,0]-H[2,0]
    H[1,1]=H[0,1]-H[2,1]
    H[1,2]=H[0,2]-H[2,2]
    det=np.linalg.det(H)
    return(H,det)

#substitute value of R3
P3=-21.49
#P2=0.0191662779286119/(-21.57)
P2=eq3.subs(R3,P3)
#P1=-0.0153214329375349
P1=eq4.subs([(R3,P3),(R2,P2)])

def ScaledResRet(e1,e2,e3,beta1,beta2,beta3):
    e0=e1*pow(beta1,-1)+P1*(e1*pow(beta1,-1)-e2*pow(beta2,-1))+P2*((e1*pow(beta1,-1)-e2*pow(beta2,-1))-P3*(e2*pow(beta2,-1)-e3*pow(beta3,-1)))
    return(e0)

#new approach starts here (4 points instead of 5)

Betas=np.array([1,0.0004501/0.002714,0.0001899/0.002714,0.00008107/0.002714])
Means=np.array([[0.002714,-0.00053],[0.0004501,-0.00078],[0.0001899,-0.00013],[0.00008107,0.00194]])
epsilon=np.zeros(4)
for i in range(0,4):
    epsilon[i]=Means[i,1]*pow(Betas[i],-1)
T1, T2, T3 =symbols('T1 T2 T3')
eq1=(epsilon[1]-epsilon[0])/(epsilon[1]-epsilon[2])+T1+T2-T3*(epsilon[2]-epsilon[3])/(epsilon[1]-epsilon[2])
eq2=-(epsilon[1]-epsilon[0])/(epsilon[2]-epsilon[3])-T1*(epsilon[1]-epsilon[2])/(epsilon[2]-epsilon[3])-T2*(epsilon[1]-epsilon[2])/(epsilon[2]-epsilon[3])+T3
eq3=eq2.subs(T1,sympy.solve(eq1,T1)[0])
t3=sympy.solve(eq3,T3)[0]
eq4=sympy.solve(eq1,T1)[0]
eq5=eq4.subs(T2,1)
t1=eq5.subs(T3,sympy.solve(eq3,T3)[0])
t2=1

H=np.zeros((3,3))
H[0,1]=(epsilon[1]-epsilon[2])*(epsilon[1]-epsilon[2])
H[0,0]=(epsilon[1]-epsilon[2])*(epsilon[1]-epsilon[2])
H[0,2]=-(epsilon[1]-epsilon[2])*(epsilon[2]-epsilon[3])
H[2,0]=-(epsilon[1]-epsilon[2])*(epsilon[2]-epsilon[3])
H[2,1]=-(epsilon[1]-epsilon[2])*(epsilon[2]-epsilon[3])
H[2,2]=(epsilon[2]-epsilon[3])*(epsilon[2]-epsilon[3])
H[1,0]=(epsilon[1]-epsilon[2])*(epsilon[1]-epsilon[2])
H[1,1]=(epsilon[1]-epsilon[2])*(epsilon[1]-epsilon[2])
H[1,2]=-(epsilon[1]-epsilon[2])*(epsilon[2]-epsilon[3])
det=np.linalg.det(H)

def ScaledResRet(e1,e2,e3,beta1,beta2,beta3):
    e0=e1*pow(beta1,-1)+t1*(e1*pow(beta1,-1)-e2*pow(beta2,-1))+t2*(e1*pow(beta1,-1)-e2*pow(beta2,-1))-t3*(e2*pow(beta2,-1)-e3*pow(beta3,-1))
    return(e0)
