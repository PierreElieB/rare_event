# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:29:16 2017

@author: Pierre-Elie
"""

import numpy as np 
import numpy.random as nprd
import matplotlib.pyplot as plt 
import math
from operator import itemgetter

print("Souhaitez-vous utiliser les paramètres par défaut ? Si oui, tapez 1, sinon, tapez 0.")
parametresParDefaut = int(input())

if(parametresParDefaut) : 
    ALPHA = float(1E-6)
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
    T1 = 0./365
    T2 = 0./365 
    BOOL_CORR = 0
    BOOL_AJOUT_PRECIPITATIONS = False 
    FACTEUR_PRECIPITATIONS = 0.1 
    BOOL_CRITERE_DU_MAXIMUM = True
else :
    print("Veuillez saisir le niveau alpha (par exemple : 1E-6)")
    ALPHA = float(input())
    print("Veuillez saisir l'horizon temporel de l'étude, en années (par exemple : 1.)")
    T = float(input()) # horizon : un an
    print("Veuillez saisir le nombre moyen de précipitations par an (par exemple : 50.)")
    NOMBRE_DE_PRECIPITATIONS_PAR_AN = float(input()) 
    print("Veuillez saisir la valeur de delta2 (l'énoncé suggère 0.7")
    DELTA2 = (float(input()))
    print("Veuillez saisir le rapport des delta, reflétant l'intensité relative des précipitations importantes par rapport aux faibles précipitations. L'énoncé suggère 10.")
    RAPPORTDELTA = float(input())
    print("Veuillez saisir la part des grandes précipitations (par exemple, B = 0.2)")
    B = float(input())
    print("Veuillez saisir le volume initiale du barrage. Unité : 10e4 m^3. Par exemple, 200.")
    V0 = float(input())
    R = 1. 
    R1 = 1. 
    R2 = 1. 
    CHANGEMENT_ECHELLE = 100
    print("Choisissez le nombre de simulations pour le calcul de chaque probabilité. Par exemple, 1000.")
    NOMBRE_SIMULATIONS = int(input())
    RATIO_CIBLE = 0.3 
    N_RATIO_OPTIMAL = 100
    NOMBRE_DIVISIONS = 20
    EPSILON3=0.01
    print("Quel est le temps de trajet (en années) de l'eau provenant du premier barrage ?")
    T1 = float(input())
    print("Quel est le temps de trajet (en années) de l'eau provenant du second barrage ?")
    T2 = float(input())
    print("Veuillez saisir le degré de corrélation entre les précipitations en amont des deux barrages. Si 0, précipitations totalement indépendantes. Si 1, précipitations identiques. Si 2, dates des précipitations identiques mais intensités indépendantes.")
    BOOL_CORR = int(input())
    print("Souhaitez-vous prendre en compte les précipitations entre le barrage et la vallée ? Si oui, tapez 1, si non, tapez False.")
    foo = int(input())
    BOOL_AJOUT_PRECIPITATIONS = bool(foo)
    FACTEUR_PRECIPITATIONS = 0.1 
    print("Souhaitez-vous vous baser sur le critère du volume maximum sur la prériode T (pour cela, tapez 1) ou sur le volume à la date T (pour cela, tapez 0) ?")
    t = int(input())    
    BOOL_CRITERE_DU_MAXIMUM = bool(t)

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
    N = nprd.poisson(lam = lambd*T)
    dates = nprd.rand(N)
    dates = T*np.sort(dates) 
    return(dates,N) 
    
# Fonction privée 

def intensitePrecipitation(del1,del2,b) : # Simule une précipitation 
    x = nprd.random()
    if(x<b) : 
        return(nprd.exponential(1./del1))
    else : 
        return(nprd.exponential(1./del2))
        
def impressionPrecipitations() : 
    (dates,N) = simulationProcessusPoisson(LAMBD) 
    z = np.zeros(N)
    for i in range(N)  :
        z[i] = intensitePrecipitation(delta1,delta2,B)
    width = 1./700.
    plt.bar(dates, z, width, color=(0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0))
    plt.grid()
    plt.ylabel('Volume précipitation (*10e5 m^3)')  
    plt.xlabel('Instant des précipitations (année)')
    plt.savefig('Precipitations.png', dpi = 720)
    plt.show()
    return()         
        
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
        plt.plot(t, y,color = 'blue')
        plt.grid()
        plt.ylabel('Volume du lac (*10e5 m^3)')  
        plt.xlabel('Temps (année)')
        plt.savefig('VolumeBarrage.png', dpi = 720)
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
    return(proba,delta) 





   
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
    boolMaxi = BOOL_CRITERE_DU_MAXIMUM
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
    return(proba,ratio/n,sigma,delta) 







    
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
        (probaI,_,_,_) = probaChangementLoi(theta,seuilCourant,boolMaxi,NOMBRE_SIMULATIONS)
        x[i+1] = seuilCourant 
        z[i+1] = probaI/ALPHA
    for i in range(NOMBRE_DIVISIONS-1) : 
        y[i+1] = z[i+1]
    y[NOMBRE_DIVISIONS-1] = z[NOMBRE_DIVISIONS-1]
    plt.plot(x,y) 
    plt.show() 
    return() 
    
# Fonction privée qui simule deux processus de Poisson indépendants, en simulant un grand processus de Poisson puis en utilisant la méthode du coloriage aléatoire. 
    
    
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
    
# Fonction publique qui calcule une trajectoire sur le débit de la vallée.     
    
def debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, imprime) : 
    if(BOOL_CORR == 1) : 
        (resMaxRouge,_,_,_) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge,False)
        resMaxBleue = resMaxRouge
        datesBleues = datesRouges
        nBleue = nRouge
        intensitesBleues = intensitesRouges
    else : 
        (resMaxRouge,_,_,_) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge, False) 
        (resMaxBleue,_,_,_)= simulationVolumeBarrage(datesBleues,intensitesBleues,nBleue, False) 
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
            if(BOOL_AJOUT_PRECIPITATIONS) : 
                debit[i] += resMaxBleue[indice] * FACTEUR_PRECIPITATIONS
        else : 
            debit[i] = R1*resMaxRouge[indice] + R2*resMaxBleue[indiceBleu]*math.exp(-R2*(datesDebit[i][0]-datesDebit[indiceBleu2][0] ))
            indiceRouge2 = i
            indiceRouge = indice
            if(BOOL_AJOUT_PRECIPITATIONS) : 
                debit[i] += resMaxRouge[indice] * FACTEUR_PRECIPITATIONS
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
        plt.grid()
        plt.ylabel('Débit en aval (*10e5 m^3/an)')  
        plt.xlabel('Temps (année)')
        plt.savefig('Debit.png', dpi = 720)
        plt.show() 
    return(np.max(debit), debit[indice1An]) 
    
# Fonction publique qui calcule la proba que le débit dans la vallée dépasse un certain débit seuil
      
      
def probaChangementLoiDeuxBarrages(theta,debitSeuil,boolMax,n) : 
    boolMax = BOOL_CRITERE_DU_MAXIMUM
    if(BOOL_CORR == 0) : 
        lambd2 = LAMBD2
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD2,B)
    if(BOOL_CORR == 0.5) : 
        lambd2 = LAMBD
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B)
    if(BOOL_CORR == 1) : 
        lambd2 = LAMBD
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B) 
    ratio = 0. 
    valeur=np.zeros(n) 
    for i in range(n) :  
        if(BOOL_CORR == 1) :  
            (datesRouges,intensitesRouges, nRouge) = simulationPrecipitations(newDelta1,newDelta2,newLambd,newB)
            datesBleues = np.copy(datesRouges)
            nBleue = nRouge
            intensitesBleues = np.copy(intensitesRouges)
        if(BOOL_CORR == 0.5) : 
            (datesRouges,intensitesRouges,nRouge,datesBleues,intensitesBleues,nBleue,intensitesAutres,NAutre) = simulationDeuxProcessusPoissonCorreles(newLambd,newDelta1,newDelta2,newB)
        if(BOOL_CORR == 0) :             
            (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(newLambd,newDelta1,newDelta2,newB)          
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, False) 
        debit0 = debitFin       
        if(boolMax) : 
            debit0 = debitMaxi
        if(debit0>debitSeuil) :  
            if(BOOL_CORR == 1) : 
                x = np.sum(intensitesRouges) 
            if(BOOL_CORR == 0) : 
                x = np.sum(intensitesBleues) + np.sum(intensitesRouges)   
            if(BOOL_CORR == 0.5) : 
                x = np.sum(intensitesBleues) + np.sum(intensitesRouges) + np.sum(intensitesAutres)
            b = -theta*x+(newLambd-lambd2)
            if(BOOL_CORR == 0.5) : 
                b = -theta*x +3.*(newLambd-lambd2)
            a = math.exp(b)
            ratio+=1.
            valeur[i]=a
        else:
            valeur[i]=0
        p=np.mean(valeur)
        sig=np.std(valeur)
        delta=1.96*sig/np.sqrt(n)
    return(p, ratio/n,sig,delta) 

#print(probaChangementLoiDeuxBarrages(newton(600),600,True,NOMBRE_SIMULATIONS))    

#Recherche du théta optimal et tests
#recherche solution gamma'(theta)=Vseuil 

#gamma'        
def dgamma(theta,Vseuil):
    return(LAMBD*T*(B*delta1/(delta1-theta)**2+(1-B)*delta2/(delta2-theta)**2)-Vseuil)


#gamma'' (pour la méthode de Newton) 
def ddgamma(theta):
    return(2*LAMBD*T*(B*delta1/(delta1-theta)**3+(1-B)*delta2/(delta2-theta)**3))
    
#Méthode de  Newton, renvoie theta optimal théorique pour un Vseuil donné   
def newton(Vseuil):    
    theta=delta1/2
    N=0
    while(dgamma(theta,Vseuil)>0.001 and N<10):
        theta=theta-dgamma(theta,Vseuil)/ddgamma(theta)
        #print(dgamma(theta,Vseuil))
        N+=1   
    #print(N)    
    return(theta)

#Recherche expérimentale de théta optimal
def intervalle_de_confiance(vSeuil,boolMaxi,n):
    boolMaxi = BOOL_CRITERE_DU_MAXIMUM
    THETA=np.linspace(0.015,0.045,10)
    SIGMA=np.zeros(10)
    for i in range(10):
        SIGMA[i]=probaChangementLoi(THETA[i],vSeuil,boolMaxi,n)[2]
    plt.plot(THETA,SIGMA)
    plt.title("Ecart-type de l'echantillon en fonction de theta")
    plt.xlabel("theta")
    plt.ylabel("Ecart-type")
    print("théta optimal théorique :")
    thetaOpt = newton(vSeuil)
    print(thetaOpt) 
    largeurIntervalle = 2* probaChangementLoi(thetaOpt, vSeuil, boolMaxi,n)[2]
    print("Intervalle de confiance :")
    print(largeurIntervalle) 
    return() 
    
    
#print(newton(559))   
#intervalle_de_confiance(559,True,20000)

#print(probaChangementLoi(newton(425),425,True,10000))
#print(proba_naif(425,True,10000))


# Fonction publique qui recherche le quantile de niveau alpha pour le remplissage du barrage.         



     
def rechercheQuantileAlpha(mini,maxi,boolPartie1,boolMaxi) : 
    boolMaxi = BOOL_CRITERE_DU_MAXIMUM
    c = (mini+maxi)/2. 
    probaMini = 1.
    while(np.abs(probaMini-ALPHA)>epsilon1) : 
        if(boolPartie1) :            
            (probaMini,_,_,_) = probaChangementLoi(newton(c),c, boolMaxi, NOMBRE_SIMULATIONS)
        else : 
            (probaMini,_,_,_) = probaChangementLoiDeuxBarrages(newton(c), c, boolMaxi, NOMBRE_SIMULATIONS)   
        if(probaMini>ALPHA) : 
            mini = c
        else : 
            maxi = c
        c = (mini+maxi)/2. 
    return(c)




def rechercheQuantileAlphaGraph(mini,maxi,boolPartie1,boolMaxi) :
    c=rechercheQuantileAlpha(mini,maxi,boolPartie1,boolMaxi)
    Q=np.linspace(c-20,c+20,20)
    P=np.zeros(20)
    E=np.zeros(20)
    if (boolPartie1) :
        for i in range(20):
            (P[i],_,_,E[i])=probaChangementLoi(newton(559),Q[i], boolMaxi, NOMBRE_SIMULATIONS)
    else :
        for i in range(20):
            print(i)
            (P[i],_,_,E[i])=probaChangementLoiDeuxBarrages(newton(559),Q[i], boolMaxi, NOMBRE_SIMULATIONS)
    P=P/ALPHA  
    E=E/ALPHA
    plt.errorbar(Q,P,xerr=0,yerr=E)
    plt.title("Recherche du Quantile Q(10**-6)")
    plt.xlabel("Q(alpha)")
    plt.ylabel("alpha(*10**-6)")
    plt.axhline(1)
    return(c)
   
    
#print(rechercheQuantileAlpha(500,1500,False,True))    
    
    
def probaConditionnelleDeuxBarrages(theta,vSeuil, debitSeuil, boolMax, n) : 
    boolMax = BOOL_CRITERE_DU_MAXIMUM
    if(BOOL_CORR == 0) : 
        lambd2 = LAMBD2
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD2,B)
    if(BOOL_CORR == 0.5) : 
        lambd2 = LAMBD
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B)
    if(BOOL_CORR == 1) : 
        lambd2 = LAMBD
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B) 
    
    ratio = 0. 
    valeur=np.zeros(n) 
    for i in range(n) : 
        if(BOOL_CORR == 1) :  
            (datesRouges,intensitesRouges, nRouge) = simulationPrecipitations(newDelta1,newDelta2,newLambd,newB)
            datesBleues = datesRouges
            nBleue = nRouge
            intensitesBleues = intensitesRouges 
        if(BOOL_CORR == 0.5) : 
            (datesRouges,intensitesRouges,nRouge,datesBleues,intensitesBleues,nBleue,intensitesAutres,NAutre) = simulationDeuxProcessusPoissonCorreles(newLambd,newDelta1,newDelta2,newB)
        if(BOOL_CORR == 0) : 
            (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(newLambd,newDelta1,newDelta2,newB)
        (_,_,_,volumeBleu) = simulationVolumeBarrage(datesBleues,intensitesBleues,nBleue,False)
        (_,_,_,volumeRouge) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge,False)
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, False) 
        debit0 = debitFin        
        if(boolMax) : 
            debit0 = debitMaxi
        if(debit0>debitSeuil and (volumeBleu>vSeuil or volumeRouge>vSeuil)) :  
            if(BOOL_CORR == 1) : 
                x = 2*np.sum(intensitesRouges) 
#            if(BOOL_CORR == 0.5) : 
#                x = np.sum(intensitesRouges) + np.sum(intensitesBleues) + np.sum(intensitesAutres)
            if(BOOL_CORR == 0) :                
                x = np.sum(intensitesBleues) + np.sum(intensitesRouges)
            b = -theta*x+(newLambd-lambd2)
#            if(BOOL_CORR == 0.5) : 
#                b = -theta*x+3*(newLambd-lambd2)
            a = math.exp(b)
            ratio+=1.
            if(BOOL_CORR == 0) : 
                valeur[i]=a/(2*ALPHA-ALPHA**2)
            if(BOOL_CORR == 1) : 
                valeur[i]=a/ALPHA   
        else:
            valeur[i]=0            
    p=np.mean(valeur)
    sig=np.std(valeur)
    delta=1.96*sig/np.sqrt(n)
    return(p, ratio/n,sig,delta) 
    
#print(probaConditionnelleDeuxBarrages(newton(559),200,845,True,10*NOMBRE_SIMULATIONS))   
    
    
def estimationVolumeBarrageDebit(debitSeuil,theta, boolMaxi,n) : 
    boolMaxi = BOOL_CRITERE_DU_MAXIMUM
    if(BOOL_CORR == 0) : 
        lambd2 = LAMBD2
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD2,B)
    if(BOOL_CORR == 0.5) : 
        lambd2 = LAMBD
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B)
    if(BOOL_CORR == 1) : 
        lambd2 = LAMBD
        (newLambd, newDelta1, newDelta2, newB) = changementLoi(theta,delta1,delta2,LAMBD,B)
    valeur=np.zeros(n)
    for i in range(n) : 
        if(BOOL_CORR == 1) :  
            (datesRouges,intensitesRouges, nRouge) = simulationPrecipitations(newDelta1,newDelta2,newLambd,newB)
            datesBleues = datesRouges
            nBleue = nRouge
            intensitesBleues = intensitesRouges 
        if(BOOL_CORR == 0.5) : 
            (datesRouges,intensitesRouges,nRouge,datesBleues,intensitesBleues,nBleue,intensitesAutres,NAutre) = simulationDeuxProcessusPoissonCorreles(newLambd,newDelta1,newDelta2,newB)
        if(BOOL_CORR == 0) : 
            (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(newLambd,newDelta1,newDelta2,newB)
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, False) 
        debit0 = debitFin       
        if(boolMaxi) : 
            debit0 = debitMaxi
        if(debit0>debitSeuil) :  
            if(BOOL_CORR == 1) : 
                x = np.sum(intensitesBleues) 
            if(BOOL_CORR == 0.5) : 
                x = np.sum(intensitesBleues) + np.sum(intensitesRouges) + np.sum(intensitesAutres)
            if(BOOL_CORR == 0) :
                x = np.sum(intensitesBleues) + np.sum(intensitesRouges)
            b = -theta*x+(newLambd-lambd2)
            if(BOOL_CORR == 0.5) : 
                b = -theta*x+3*(newLambd-lambd2)
            a = math.exp(b)
            (_,_,vMaxi,_) = simulationVolumeBarrage(datesRouges,intensitesRouges,nRouge,False)
            valeur[i]= a*vMaxi/ALPHA
        else:
            valeur[i]=0
    p=np.mean(valeur)
    sig=np.std(valeur)
    delta=1.96*sig/np.sqrt(n)            
    return(p,sig,delta) 

 
   
def repartitionDeuxBarrages(debitMini,debitMaxi,vMini,vMaxi,theta, boolMaxi) : 
    boolMaxi = BOOL_CRITERE_DU_MAXIMUM
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
   
   
def simulationDeuxProcessusPoissonCorreles(lambd,delta1,delta2,b) :
    (dates,N1) = simulationProcessusPoisson(lambd) 
    (_,N2) = simulationProcessusPoisson(lambd) 
    (_,N3) = simulationProcessusPoisson(lambd)
    intensitesRouges = np.zeros(N1) 
    intensitesBleues = np.zeros(N1) 
    NAutres = N2+N3-N1
    intensitesAutres = np.zeros(NAutres) 
    for i in range(N1) : 
        intensitesRouges[i] = intensitePrecipitation(delta1,delta2,b) 
        intensitesBleues[i] = intensitePrecipitation(delta1,delta2,b) 
    for i in range(NAutres) : 
        intensitesAutres[i] = intensitePrecipitation(delta1,delta2,b) 
    return(dates,intensitesRouges,N1,dates,intensitesBleues,N1,intensitesAutres,NAutres)
#    
while(True) : 
    print("Que souhaitez-vous faire ?")
    print("Si vous avez terminé, tapez 0.")
    print("Pour réaliser une simulation de l'évolution du volume du barrage au cours du temps, tapez 1.")
    print("Pour calculer le quantile de niveau alpha pour un barrage, tapez 2.")
    print("Pour calculer le quantile de niveau alpha pour le débit en vallée avec deux barrages, tapez 3.")
    print("Pour réaliser une simulation de l'évolution du débit à la confluence au cours du temps, tapez 4.")

    aFaire = int(input())
    
    if(aFaire == 0) : 
        break

    if(aFaire == 1) :
        print("Simulation d'une trajectoire pour le volume du barrage.")
        (dates,intensites,N) = simulationPrecipitations(delta1,delta2,LAMBD, B)
        simulationVolumeBarrage(dates,intensites,N,True) 
    
    if(aFaire == 2) : 
        print("Quantile naïf :") 
        print(quantileNaif()) 
        print()
        print("Quantile en utilisant les méthodes de simulation des événements rares :")
        q = rechercheQuantileAlpha(400., 1200., True, True)
        print(q)
        print()
        print("Vérification de la probabilité :")
        (p,_,_,delta)=probaChangementLoi(newton(559),q, True, NOMBRE_SIMULATIONS)
        print(p)
        print()
        print("Largeur de l'intervalle de confiance pour la probabilité :")
        largeur = 2.*delta
        print(largeur)
        print()
        print()
        
    if(aFaire == 3) : 
        print("Recherche quantile alpha pour le débit en vallée :") 
        q = rechercheQuantileAlpha(400., 1200., False, True)
        print(q)
        print()
        print("Vérification de la probabilité :")
        (p,_,_,delta) = probaChangementLoiDeuxBarrages(newton(559),q,True,NOMBRE_SIMULATIONS)
        print(p)
        print("Largeur de l'intervalle de confiance pour la probabilité :")
        largeur = 2.*delta
        print(largeur)
        print()
        print()
        
    if(aFaire == 4) : 
        print("Simulation d'une trajectoire pour le débit à la confluence.")
        (datesBleues, intensitesBleues, nBleue, datesRouges,intensitesRouges,nRouge) = simulationDeuxProcessusPoissonIndependants(LAMBD2,delta1,delta2,B)          
        (debitMaxi,debitFin) = debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, True) 
    
        
print()
print()
print()      
print("Nous vous remercions d'avoir utilisé cette application console développée en Python par Thomas et Pierre-Elie dans le cadre du Modal sur les événements rares !") 
   
   
#print()
#     
#(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge) = simulationDeuxProcessusPoissonIndependants(2*LAMBD,delta1,delta2,B)
#print(debitVallee(datesBleues,intensitesBleues, nBleue, datesRouges, intensitesRouges, nRouge, True))         
#    
#
#repartitionUnBarrage(400.,1200.,thetaMax*0.55,True)        
#print(estimationVolumeBarrageDebit(500.,1500.,0.5*thetaMax,True,500)) 
        
        
    
    