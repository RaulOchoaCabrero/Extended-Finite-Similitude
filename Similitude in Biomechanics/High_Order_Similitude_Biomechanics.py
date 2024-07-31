# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:11:24 2020

@author: Raul Ochoa Cabrero

Zeroth Order & First Order Application Biomechanics

Copyright: Copyright (c) 2020 Raul Ochoa Cabrero

"""
# Definitions of scaling relationships and parameters.
def alpharho(beta,rho1,rhobeta):
    ar=rho1/(pow(beta,3)*rhobeta)
    return(ar)

def matpropt(rhobetac,rho1c,rho1t):
    rbt=rhobetac/rho1c*rho1t
    return(rbt)

def g(beta,rho1,rhobeta,E1,Ebeta):
    g=sqrt((E1*rhobeta*pow(beta,2))/(rho1*Ebeta))
    return(g)

def alphanu(beta,rho1,rhobeta,E1,Ebeta):
    av=alpharho(beta, rho1, rhobeta)*pow(beta,-1)*g(beta,rho1,rhobeta,E1,Ebeta)
    return(av)
    
def Forceb1(F1,beta1,beta2,rho1,rhobeta1,rhobeta2,E1,Ebeta1,Ebeta2,R1nu,C):
    Fb1=F1/(alphanu(beta1,rho1,rhobeta1,E1,Ebeta1)*g(beta1,rho1,rhobeta1,E1,Ebeta1)*(1+R1nu)-R1nu*alphanu(beta2,rho1,rhobeta2,E1,Ebeta2)*g(beta2,rho1,rhobeta2,E1,Ebeta2)*C)
    return(Fb1,C*Fb1)

def FirstOrderBio(F1,beta1,beta2,rho1,rhobeta1,rhobeta2,E1,Ebeta1,Ebeta2,R1nu,Fb1,rho1t,E1t):
    ar1=rho1/(pow(beta1,3)*rhobeta1)
    ar2=rho1/(pow(beta2,3)*rhobeta2)
    g1=sqrt((E1*rhobeta1*pow(beta1,2))/(rho1*Ebeta1))
    g2=sqrt((E1*rhobeta2*pow(beta2,2))/(rho1*Ebeta2))
    anu1=ar1*g1*pow(beta1,-1)
    anu2=ar2*g2*pow(beta2,-1)
    Fb2=(F1-(R1nu+1)*anu1*g1*Fb1)/(-R1nu*anu2*g2)
    rhob1t=rhobeta1/rho1*rho1t
    rhob2t=rhobeta2/rho1*rho1t
    Eb1t=Ebeta1/E1*E1t
    Eb2t=Ebeta2/E1*E1t
    Parameters=np.array([[ar1,ar2],[g1,g2],[anu1,anu2],[Fb1,Fb2],[rhob1t,rhob2t],[Eb1t,Eb2t]])
    return(Parameters.astype('float64'))

def ZerothOrderBio(F1,beta,rho1,rho1t,rhobeta,E1,E1t,Ebeta):
    ar=rho1/(pow(beta,3)*rhobeta)
    g=sqrt(E1/(beta*ar*Ebeta))
    anu=ar*g*pow(beta,-1)
    Fbeta=F1/(anu*g)
    rhobbetat=rhobeta/rho1*rho1t
    Ebetat=Ebeta/E1*E1t
    Parameters=np.array([[ar],[g],[anu],[Fbeta],[rhobbetat],[Ebetat]])
    return(Parameters.astype('float64'))

def ZerothOrderReverse(Fbeta,beta,rho1,rhobeta,rhobetat,E1,Ebeta,Ebetat):
    ar=rho1/(pow(beta,3)*rhobeta)
    g=sqrt((E1*rhobeta*pow(beta,2))/(rho1*Ebeta))
    anu=ar*g*pow(beta,-1)
    F1=Fbeta*anu*g
    E1t=Ebetat*E1/Ebeta
    rho1t=rhobetat*rho1/rhobeta
    Parameters=np.array([[pow(ar,-1)],[pow(g,-1)],[pow(anu,-1)],[F1],[rho1t],[E1t]])
    return(Parameters.astype('float64'))

# Application of similitude relationships to the experimental parameters
FirstOrderBio(2000,4/5,3/5,0.00000164,0.00000125,0.00000125,16700,3480,3286,.5,266.73,0.00000027,155)
ZerothOrderBio(2000,4/5,0.00000164,.00000027,0.00000125,16700,155,3286)
ZerothOrderReverse(160,3/5,.00000164,.00000125,0.00000121,16700,3286,20.84)

# Definition of Force scaling relationships
def FirstOrderForces(R1nu,a1,s,n,C,rhobeta2t,Ebeta2t,Ebeta1t):
    a=np.zeros(n)
    for i in range(0,n):
        a[i]=a1+(n-1-i)*s
    c=[FirstOrderBio(2000,4/5,3/5,0.00000164,0.00000125,0.00000125,16700,3480,3286,R1nu,C*a[0],0.00000027,155)[3,:]]
    for i in range(1,n):
        c=np.append(c,[FirstOrderBio(a[i],4/5,3/5,0.00000164,0.00000125,0.00000125,16700,3480,3286,R1nu,C*a[i],0.00000027,155)[3,:]],axis=0)
    d=[]
    for i in range(0,n):
        d=np.append(d,[ZerothOrderReverse(c[i,1],3/5,0.00000164,0.00000125,rhobeta2t,16700,3286,Ebeta2t)[3]])
    b=[]
    for i in range(0,n):
        b=np.append(b,[ZerothOrderReverse(c[i,0],4/5,0.00000164,0.00000125,rhobeta2t,16700,3480,Ebeta1t)[3]])
    return(a.reshape(-1,1),c,d.reshape(-1,1),b.reshape(-1,1))

# Scaling of forces with experimental parameters
FirstOrderForces(-.5,1200,100,9,.14,0.00000121,20.84,43.83)
x1=np.transpose(FirstOrderForces(-.5,1200,100,9,.14,0.00000121,20.84,43.83)[2])
x1=x1[0]
x2=np.transpose(FirstOrderForces(-.5,1200,100,9,.14,0.00000121,20.84,43.83)[3])
x2=x2[0]

"""
The next section requires importing the strain measurements from experimental and simulated data.

Read csv files and define strain tensor for each value of the force
"""
# Statistical optimisation of free parameter R1nu
s1=0
s2=0
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,6):
            s1=s1+(FEstrain[k,i,j]-FEstrain[k,i,j+6])*(FEstrain[k,i,j+6]-FEstrain[k,i,j+12])
            s2=s2+(-FEstrain[k,i,j+6]+FEstrain[k,i,j+12])*(-FEstrain[k,i,j+6]+FEstrain[k,i,j+12])
            
StatisticalR1nu=s1/s2

# Average Differences and Average Percentage Errors 
Difference=np.zeros((9,Length[0],6))
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,6):
            Difference[k,i,j]=abs(FEstrain[k,i,j]-((StatisticalR1nu+1)*FEstrain[k,i,j+6]-StatisticalR1nu*FEstrain[k,i,j+12]))
            
[Numbers,Bins]=np.histogram(Difference,bins='auto')
for i in range(0,len(Numbers)):
    if sum(Numbers[i:])/(sum(Numbers[0:i])+.001)<0.01:
        Number95=i
        break

Difference95Sum=0
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,6):
            if Difference[k,i,j]<Bins[Number95]:
                Difference95Sum=Difference95Sum+(Difference[k,i,j])
                
AverageDifference=Difference95Sum/(sum(Numbers[0:Number95])-sum(Numbers[Number95:]))

Difference95Sum=0
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,6):
            Difference95Sum=Difference95Sum+(Difference[k,i,j])
        
AverageDifference=Difference95Sum/(Length[0]*Length[1]*9)

PercentageErrorSum=np.zeros(9)
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,6):
            if Difference[k,i,j]<Bins[Number95]:
                PercentageErrorSum[k]=PercentageErrorSum[k]+(Difference[k,i,j]/max(abs(FEstrain[k,i,j]),.000001))

PercentageErrorSum=np.zeros(9)
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,6):
            PercentageErrorSum[k]=PercentageErrorSum[k]+(Difference[k,i,j]/max(abs(FEstrain[k,i,j]),.000001))

AveragePercentageError=PercentageErrorSum/(Length[0]*6)
                                                         
SumStrain=0
for k in range(0,9):
    for i in range(0,Length[0]):
        for j in range(0,Length[1]):
            SumStrain=SumStrain+abs(FEstrain[k,i,j])
            
AverageStrain=SumStrain/np.prod(FEstrain.shape)

sum(Difference[2,:,:])/Length[0]

Difference95Sum=np.zeros(7)
AverageDifference=np.zeros(10)

for k in range(0,7):
   for i in range(floor(Length[0]/7)*k,floor(Length[0]/7)*(k+1)):
       for j in range(0,6):
           if Difference[1,i,j]<Bins[Number95]:
               Difference95Sum[k]=Difference95Sum[k]+(Difference[1,i,j])

for i in range(0,7): 
    AverageDifference[i]=Difference95Sum[i]/((sum(Numbers[0:Number95])-sum(Numbers[Number95:]))/7)
    
