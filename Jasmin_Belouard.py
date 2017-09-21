# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:29:16 2017

@author: Pierre-Elie
"""

import numpy as np 
import numpy.random as nprd
import matplotlib.pyplot as plt 
import scipy.stats
import math
from operator import itemgetter

NOMBREBARRAGES = 1
ALPHA = 1E-6
T = 1. # horizon : un an
NOMBRE_DE_PRECIPITATIONS_PAR_AN = 50. 
DELTA2 = 0.7
RAPPORTDELTA =  10.
B = 0.2 
V0 = 200.
R = 1. 
R1 = 1. 
R2 = 1. 
CHANGEMENT_ECHELLE = 100
NOMBRE_SIMULATIONS = 1000
RATIO_CIBLE = 0.3 
N_RATIO_OPTIMAL = 100
NOMBRE_DIVISIONS = 20
EPSILON3=0.01
T1 = 40./365
T2 = 80./365 


# Variables globales non entrées

LAMBD = NOMBRE_DE_PRECIPITATIONS_PAR_AN
LAMBD2 = 2.*LAMBD
delta2 = DELTA2
delta1 = delta2*(1./RAPPORTDELTA) 
epsilon1 = ALPHA/10.
epsilon2 = RATIO_CIBLE/20. 
thetaMin = 0. 
thetaMax = delta1-EPSILON3

# Fonction privée 

def simulationProcessusPoisson(lambd) : # Simule un processus de Poisson de paramètre lambda sur la période T
    N = nprd.poisson(lam = lambd)
    dates = nprd.rand(N)
    dates = np.sort(dates) 
    return(dates,N) 
    
# Fonction privée 

def intensitePrecipitation(del1,del2,b) : # Simule une précipitation 
    x = nprd.random()
    if(x<b) : 
        return(nprd.exponential(1./del1))
    else : 
        return(nprd.exponential(1./del2))
        
# Fonction semi-privée : le résultat importe peu, mais il faut appeler la fonction avant d'appeler simulationVolumeBarrage
        
def simulationPrecipitations(del1,del2,lamb,b) : # Simule l'ensemble des précipitations sur la période [0,T] et renvoie le tableau des dates des précipitations, les quantités de pluie à chaque épisode et le nombre d'épisodes de précipitations
    (dates,N) = simulationProcessusPoisson(lamb)
    intensites = np.zeros(N) 
    for i in range(N) : 
        intensites[i] = intensitePrecipitation(del1,del2,b) 
    return(dates,intensites,N) 
  
# Fonction à utiliser après avoir appelée celle ci-dessus. Elle permet de tracer des trajectoires de remplissage du barrage sur un an.   
        
def simulationVolumeBarrage(dates,intensites,N,doitImprimer) :   # Effectue une simulation de l'évolution du volume d'eau contenu dans le barrage sur la période [0,T]
    resMax = np.zeros(N+1) 
    resMin = np.zeros(N+1)
    resMax[0] = V0
    resMin[0] = V0*math.exp(-R*dates[0])    
    for i in range(1,N) : 
        resMax[i] = resMin[i-1]+intensites[i-1]
        resMin[i] = resMax[i]*math.exp(-R*(dates[i]-dates[i-1]))        
    resMax[N] = resMin[N-1]+intensites[N-1]
    resMin[N] = resMax[N]*math.exp(-R*(T-dates[N-1]))
    if(doitImprimer) : 
        t = np.zeros(2*N+1) 
        y = np.zeros(2*N+1)
        t[0] = 0.
        t[1] = 0.5*dates[0] 
        y[0] = resMax[0] 
        y[1] = resMin[0] 
        for i in range(1,N) : 
            t[2*i] = dates[i-1] 
            t[2*i+1] = (dates[i-1] + dates[i]) /2. 
            y[2*i] = resMax[i]
            y[2*i+1] = resMin[i]
        t[2*N] = T
        y[2*N] = resMax[N]
        plt.plot(t,y) 
        plt.show()         
    return(resMax,resMin,np.max(resMax), resMin[N])

#Probabilité d'atteinte d'un seuil avec son intervalle de confiance, méthode naïve 
 
def proba_naif(vSeuil,boolMaxi,n):
    valeurs=np.zeros(n)
    for i in range(n):
        (dates,intensites,N) = simulationPrecipitations(delta1,delta2,LAMBD,B)
        (_,_,vMax,vFin) = simulationVolumeBarrage(dates,intensites,N, False)
        v0 = vFin    
        if(boolMaxi) : 
            v0 = vMax
        valeurs[i]=(v0>vSeuil)
    proba=np.mean(valeurs)        
    sigma=np.std(valeurs)
    delta=1.96*sigma/np.sqrt(n) #intervalle de confiance 
    print(delta)       
    return(proba) 



   
# Fonction publique qui calcule naïvement le quantile de niveau alpha

def quantileNaif() : # Donne de manière naïve la valeur du quantile de niveau alpha. 
    simulations = np.zeros(NOMBRE_SIMULATIONS) 
    for i in range(NOMBRE_SIMULATIONS) : 
        (dates,intensites,N) = simulationPrecipitations(delta1,delta2,LAMBD, B)
        (resMax, resMin,vMaxi,vFin) = simulationVolumeBarrage(dates,intensites,N,False) 
        simulations[i] = vMaxi    
    p = int(NOMBRE_SIMULATIONS*ALPHA)
    simulations = np.sort(simulations) 
    return(simulations[NOMBRE_SIMULATIONS-1-p]) #??
   
# Fonction privée pour faire un changement de loi
   
   
 
   
   
def changementLoi(theta,delta1,delta2,lambd,b) :   # Un changement de loi de paramètre thêta. 
    newLambd = lambd * b * delta1 /(delta1-theta) + lambd *(1-b) * delta2/(delta2-theta) 
    newDelta1 = delta1-theta
    newDelta2 = delta2-theta
    newB = b*lambd*delta1 /(newLambd*newDelta1) 
    return(newLambd, newDelta1, newDelta2, newB)    
       
# Fonction publique qui calcule la valeur qu'un volume dépasse un volume seuil.      
    
def probaChangementLoi(theta,vSeuil,boolMaxi,n) : 
    (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B) 
    ratio = 0. 
    valeurs=np.zeros(n) #stockage des valeurs pour intervalle de confiance
    for i in range(n) : 
        (dates,intensites,N) = simulationPrecipitations(newDelta1,newDelta2,newLambd, newB)
        (_,_,vMax,vFin) = simulationVolumeBarrage(dates,intensites,N, False)
        v0 = vFin    
        if(boolMaxi) : 
            v0 = vMax
        if(v0>vSeuil) :  
            x = np.sum(intensites) 
            b = -theta*x+(newLambd-LAMBD)
            a = math.exp(b)
            valeurs[i]=a
            ratio+=1
        else :    
            valeurs[i]=0
    proba=np.mean(valeurs)
    sigma=np.std(valeurs)
    delta=1.96*sigma/np.sqrt(n) #intervalle de confiance  
    print(delta)      
    return(proba, ratio/n,sigma) 

#Recherche du théta optimal et tests
#recherche solution gamma'(theta)=Vseuil 
        
def dgamma(theta,Vseuil):
    return(LAMBD*T*(B*delta1/(delta1-theta)**2+(1-B)*delta2/(delta2-theta)**2)-Vseuil)
 
def ddgamma(theta):
    return(2*LAMBD*T*(B*delta1/(delta1-theta)**3+(1-B)*delta2/(delta2-theta)**3))
    
def newton(Vseuil):    
    theta=delta1/2
    N=0
    while(dgamma(theta,Vseuil)>0.001 and N<10):
        theta=theta-dgamma(theta,Vseuil)/ddgamma(theta)
        #print(dgamma(theta,Vseuil))
        N+=1
    #print(N)    
    return(theta)

def intervalle_de_confiance(vSeuil,boolMaxi,n):
    THETA=np.linspace(0.01,0.06,20)
    SIGMA=np.zeros(20)
    for i in range(20):
        SIGMA[i]=probaChangementLoi(THETA[i],vSeuil,boolMaxi,n)[2]
    plt.plot(THETA,SIGMA)
    print(newton(vSeuil))
    
#print(newton(400))   
intervalle_de_confiance(400,True,1000)



    
# Fonction privée qui n'est pas utilisée.     
    
def rechercheThetaOptimal(thetaMin,thetaMax,ratioOptimal,vSeuil,questionA) : 
    thetaC = (thetaMin+thetaMax)/2. 
    ratioCourant = 0. 
    while(np.abs(ratioCourant-RATIO_CIBLE)>epsilon2) : 
        (_,ratioCourant) = probaChangementLoi(thetaC,vSeuil,questionA,N_RATIO_OPTIMAL)  
        if(ratioCourant<RATIO_CIBLE) : 
            thetaMax = thetaC
        else : 
            thetaMin = thetaC
        thetaC = (thetaMin+thetaMax)/2. 
    return(thetaC)
    
# Fonction publique qui recherche le quantile de niveau alpha pour le remplissage du barrage.         
     
def rechercheQuantileAlpha(mini,maxi,theta,boolPartie1,boolMaxi) : 
    c = (mini+maxi)/2. 
    probaMini = 1.
    while(np.abs(probaMini-ALPHA)>epsilon1) : 
        if(boolPartie1) : 
            (probaMini,_) = probaChangementLoi(theta,c, boolMaxi, NOMBRE_SIMULATIONS)
        else : 
            (probaMini,_) = probaChangementLoiDeuxBarrages(theta, c, boolMaxi, NOMBRE_SIMULATIONS)
        if(probaMini>ALPHA) : 
            mini = c
        else : 
            maxi = c
        c = (mini+maxi)/2. 
    return(c)
    
# Fonction publique qui calcule la fonction de répartition inverse dans le cas où le volume dépasse le volume seuil.     
    
def repartitionUnBarrage(mini,maxi,theta, boolMaxi) : 
    c = rechercheQuantileAlpha(mini, maxi, theta, True,boolMaxi)
    x = np.zeros(NOMBRE_DIVISIONS+1)
    y = np.zeros(NOMBRE_DIVISIONS+1) 
    z = np.zeros(NOMBRE_DIVISIONS+1)
    pas = 0.5*c*1./(NOMBRE_DIVISIONS+1)
    seuilCourant = c
    x[0] =  c
    y[0] = 1. 
    for i in range(NOMBRE_DIVISIONS) : 
        seuilCourant += pas
        (probaI,_) = probaChangementLoi(theta,seuilCourant,boolMaxi,NOMBRE_SIMULATIONS)
        x[i+1] = seuilCourant 
        z[i+1] = probaI/ALPHA
    for i in range(NOMBRE_DIVISIONS-1) : 
        y[i+1] = z[i+1]
    y[NOMBRE_DIVISIONS-1] = z[NOMBRE_DIVISIONS-1]
    plt.plot(x,y) 
    plt.show() 
    return() 
    
# Fonction privée qui simule deux processus de Poisson indépendants. 
    
    
def simulationDeuxProcessusPoissonIndependants(lambd,del1,del2,b) : 
    (dates,intensites,N) = simulationPrecipitations(del1,del2,lambd,b)
    datesBleues = []
    intensitesBleues = []
    datesRouges = []
    intensitesRouges = [] 
    nBleue = 0
    nRouge = 0
    for i in range(N) : 
        x = np.random.rand() 
        if(x<0.5) : 
            intensitesBleues.append(intensites[i])
            datesBleues.append(dates[i])
            nBleue +=1
        else : 
            intensitesRouges.append(intensites[i])
            datesRouges.append(dates[i])
            nRouge +=1
    intensitesRouges = np.array(intensitesRouges)
    intensitesBleues = np.array(intensitesBleues) 
    datesRouges = np.array(datesRouges) 
    datesBleues = np.array(datesBleues) 
    return(datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge)
    
# Fonction publique qi calcule une trajectoire sur le débit de la vallée.     
    
def debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, imprime) : 
    (resMaxRouge,_,_,_) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge, False) 
    (resMaxBleue,_,_,_)= simulationVolumeBarrage(datesBleues,intensitesBleues,nBleue, False) 
    datesRouges = datesRouges
    datesBleues = datesBleues
    datesRougesDecalees = [] 
    datesBleuesDecalees = [] 
    for i in range(nBleue) : 
        datesBleuesDecalees.append((datesBleues[i]+T1,True,i))
    for i in range(nRouge) : 
        datesRougesDecalees.append((datesRouges[i]+T2,False,i))
    datesDebit = [(0.,True,0)] + datesBleuesDecalees + datesRougesDecalees 
    datesDebit.sort(key=itemgetter(0))
    nombreDates = nRouge+nBleue+1
    debit = np.zeros(nombreDates) 
    indiceRouge = 0
    indiceRouge2 = 0
    indiceBleu = 0
    indice1An = 0
    for i in range(nombreDates) : 
        (temps,estBleue,indice) = datesDebit[i]
        if(estBleue) : 
            debit[i] = R1*resMaxBleue[indice] + R2*resMaxRouge[indiceRouge]*math.exp(-R2*(datesDebit[i][0]-datesDebit[indiceRouge2][0] ))
            indiceBleu2 = i
            indiceBleu = indice
        else : 
            debit[i] = R1*resMaxRouge[indice] + R2*resMaxBleue[indiceBleu]*math.exp(-R2*(datesDebit[i][0]-datesDebit[indiceBleu2][0] ))
            indiceRouge2 = i
            indiceRouge = indice
    datesDebit2 = np.zeros(nombreDates) 
    indice1An = 0
    booleenTrouve = False 
    for i in range(nombreDates) : 
        datesDebit2[i] = datesDebit[i][0]
        if(datesDebit[i][0] >T) and (not(booleenTrouve)): 
            booleenTrouve = True 
            indice1An = i            
    if(imprime) : 
        plt.plot(datesDebit2,debit) 
        plt.show() 
    return(np.max(debit), debit[indice1An]) 
    
# Fonction publique qui calcule la proba que le débit dans la vallée dépasse un certain débit seuil
      
      
def probaChangementLoiDeuxBarrages(theta,debitSeuil,boolMax,n) : 
    (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD2,B) 
    res = 0.
    ratio = 0. 
    for i in range(n) : 
        (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(newLambd,newDelta1,newDelta2,newB)
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, False) 
        debit0 = debitFin       
        if(boolMax) : 
            debit0 = debitMaxi
        if(debit0>debitSeuil) :  
            x = np.sum(intensitesBleues) + np.sum(intensitesRouges)
            b = -theta*x+(newLambd-LAMBD2)
            a = math.exp(b)
            res = res + a
            ratio+=1.
    return(res/n, ratio/n) 
    
def probaConditionnelleDeuxBarrages(mini, maxi, theta, vSeuil, debitSeuil, boolMax, n) : 
    (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD2,B) 
    res = 0.
    ratio = 0. 
    for i in range(n) : 
        (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(newLambd,newDelta1,newDelta2,newB)
        (_,_,_,volumeBleu) = simulationVolumeBarrage(datesBleues,intensitesBleues,nBleue,False)
        (_,_,_,volumeRouge) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge,False)
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, False) 
        debit0 = debitFin        
        if(boolMax) : 
            debit0 = debitMaxi
        if(debit0>debitSeuil and (volumeBleu>vSeuil or volumeRouge>vSeuil)) :  
            x = np.sum(intensitesBleues) +np.sum(intensitesRouges)
            b = -theta*x+(newLambd-LAMBD2)
            a = math.exp(b)
            res = res + a
            ratio+=1.
    return(res/(ALPHA*n),ratio/n) 
    
    
def estimationVolumeBarrageDebit(debitMini,debitMaxi,theta, boolMaxi,n) : 
    debitSeuil = rechercheQuantileAlpha(debitMini,debitMaxi,theta,boolMaxi)
    (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD2,B) 
    res = 0.
    for i in range(n) : 
        (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(newLambd,newDelta1,newDelta2,newB)
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, False) 
        debit0 = debitFin       
        if(boolMaxi) : 
            debit0 = debitMaxi
        if(debit0>debitSeuil) :  
            x = np.sum(intensitesBleues) + np.sum(intensitesRouges)
            b = -theta*x+(newLambd-LAMBD2)
            a = math.exp(b)
            (_,_,vMaxi,_) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge,False)
            res = res + a*vMaxi
    return(res/(ALPHA*n)) 
    
    
def repartitionDeuxBarrages(debitMini,debitMaxi,vMini,vMaxi,theta, boolMaxi) : 
    debitSeuil = rechercheQuantileAlpha(debitMini, debitMaxi, theta, False, boolMaxi)
    x = np.zeros(NOMBRE_DIVISIONS+1)
    y = np.zeros(NOMBRE_DIVISIONS+1) 
    z = np.zeros(NOMBRE_DIVISIONS+1)
    pas = (vMaxi-vMini)/(NOMBRE_DIVISIONS+1)
    volumeCourant = vMini
    for i in range(NOMBRE_DIVISIONS+1) : 
        volumeCourant += pas
        (probaI,_) = probaConditionnelleDeuxBarrages(debitMini,debitMaxi,theta,volumeCourant,debitSeuil,boolMaxi,NOMBRE_SIMULATIONS)
        x[i+1] = volumeCourant 
        z[i+1] = probaI/ALPHA
    for i in range(NOMBRE_DIVISIONS-1) : 
        y[i+1] = z[i+1]
    y[NOMBRE_DIVISIONS-1] = z[NOMBRE_DIVISIONS-1]
    plt.plot(x,y) 
    plt.show() 
    return() 
   
#   
#(dates,intensites,N) = simulationPrecipitations(delta1,delta2,LAMBD, B)
#simulationVolumeBarrage(dates,intensites,N,True) 

#print("Quantile Naïf :")  
#print(quantileNaif()) 
#
#print("Recherche quantile alpha pour le volume : ")
#print(rechercheQuantileAlpha(400., 1200., 0.5*thetaMax, True, True))
#print("Vérification grossière du quantile alpha pour le volume") 
#print(probaChangementLoi(thetaMax*0.55,555.,True,1000)) 
#print("Recherche quantile alpha pour le débit en vallée :") 
#print(rechercheQuantileAlpha(400., 1200., 0.5*thetaMax, False, True))
#print("Vérification grosière pour le quantile alpha du débit en vallée : ")
#print(probaChangementLoiDeuxBarrages(thetaMax*0.55, 850., True, 1000)) 
#print()
#     
#(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge) = simulationDeuxProcessusPoissonIndependants(2*LAMBD,delta1,delta2,B)
#print(debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, True))         
#    
#
#repartitionUnBarrage(400.,1200.,thetaMax*0.55,True)        
#print(estimationVolumeBarrageDebit(500.,1500.,0.5*thetaMax,True,500)) 
        
        
# Un théta qui marche bien : 0.5*thetaMax. 
    
    