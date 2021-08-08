# from mpl_toolkits.mplot3d import Axes3D 

import numpy as np
import matplotlib.pyplot as plt


#----- Definition des constantes -----# 
   
a=5 # coefficients a et b des ellipses
b=4.7  #5 4.7
n=12  # Nombre d'ellipses sur le cone
portionCyc=[0,2*np.pi]# Domaine des portions de cyclides
#portionCyc=[0,(3/2)*np.pi]  

ThetaY=-0.4  # Angle d'inclinaison, -0.18, 0.80
# ThetaYC=0.18
# ThetaY > -pi/2 -arctan(-b/c)


#----- Choix des plots a afficher -----#

#True pour afficher
#False pour ne pas afficher

nAff=n  #  Nombre d'ellipses a afficher
mList1=[3]   #mList1=[0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
mList2=[3,9,15,21]          #mList2=[0, 3, 6] mList2=[3,6,9]
                                       # mList2=[3,7,11,15,19,23,27]
                                       # mList2=[3,9,15,21,27] 
afficher_ellipse=True  

afficher_Cyclides1=False 
afficher_Bord1=False
afficher_Cones1=False
afficher_Hyperbole1=False

afficher_Cyclides2=True   
afficher_Bord2=False 
afficher_Cones2=False
afficher_Hyperbole2=True
#----- Choix des couleurs (sauf cyclides) -----#

# couleurEllipse=(0.2,0.2,1)
couleurEllipse='black'
# couleurHyperbole1=(1,0,0)
couleurHyperbole1='orange'
couleurBord1=(0,0,0)
couleurCone1=(0.8,0.8,0.8)
alphaCone1=0.15
alphaCyclides1=1

couleurHyperbole2=couleurHyperbole1
couleurBord2=couleurBord1
couleurCone2=couleurCone1
alphaCone2=alphaCone1
alphaCyclides2=alphaCyclides1

epaisseurEllipse=2
epaisseurHyperbole=2
epaisseurBord=2.5

 
#============== Definition des fonctions ==============#

# nR Valeur entre 0 et 2 :
#    0 : rotation selon Ox,
#    1 : rotation selon Oy ,
#    2:  rotation selon Oz. 

# Rotation d'une surface
def Rot3DSurf(surface,angle,nR):
    c=np.cos(angle)
    s=np.sin(angle)
    RX = np.array( ((1,0,0),(0,c,-s),(0,s, c) ))
    RY=np.array( ((c,0,s),(0,1,0),(-s,0, c) ))
    RZ=np.array( ((c,-s,0),(s,c,0),(0,0, 1) ))
    MatRot=[RX,RY,RZ]
    surface_ROT=np.zeros(np.shape(surface))
    for j in np.arange(0,np.shape(surface)[1]):
        surface_rot_fun=[surface[0][j],surface[1][j],surface[2][j]]
        surface_rot_fun2=np.zeros(np.shape(surface_rot_fun))
        surface_rot_fun3=np.zeros(3)
        for i in np.arange(0,np.shape(surface)[2]):
            surface_rot_fun3=[surface_rot_fun[0][i],surface_rot_fun[1][i],surface_rot_fun[2][i]]
            surface_rot_fun2[:,i]=np.dot(surface_rot_fun3,MatRot[nR])
        surface_ROT[0][j]=surface_rot_fun2[0]
        surface_ROT[1][j]=surface_rot_fun2[1]
        surface_ROT[2][j]=surface_rot_fun2[2]
    return surface_ROT

# Rotation d'une courbe
def Rot3D(Ellipse,angle,nR):
    c=np.cos(angle)
    s=np.sin(angle)
    RX = np.array( ((1,0,0),(0,c,-s),(0,s, c) ))
    RY=np.array( ((c,0,s),(0,1,0),(-s,0, c) ))
    RZ=np.array( ((c,-s,0),(s,c,0),(0,0, 1) ))
    MatRot=[RX,RY,RZ]
    Ellipse_ROT=np.zeros(np.shape(Ellipse))
    for i in np.arange(0,np.shape(Ellipse)[1]):
        Ellipse_fun=[Ellipse[0][i],Ellipse[1][i],Ellipse[2][i]]
        Ellipse_ROT[:,i]=np.dot(Ellipse_fun,MatRot[nR])
    return Ellipse_ROT

# def EqCyclide(U,V,m):
#     x = c*np.cosh(V)-( (a*np.cosh(V)-m)*(c*np.cosh(V)-a*np.cos(U)) )/( a*np.cosh(V)-c*np.cos(U) )
#     y = (( a*np.cosh(V)-m )*b* np.sin(U))/(a* np.cosh(V)-c*np.cos(U)) 
#     z = b*np.sinh(V)-( ( a*np.cosh(V)-m )*b*np.sinh(V) ) / ( a*np.cosh(V)-c*np.cos(U) )
#     return [x,y,z]

def EqCyclide(U,V,m):
    
    x=c*np.cosh(V)+( (-m+a*np.cosh(V))*(-a*np.cos(U)-c*np.cosh(V)) )/( c*np.cos(U) + a*np.cosh(V) )
    y= ( b*( -m + a*np.cosh(V) ) * np.sin(U) ) / ( c*np.cos(U) + a* np.cosh(V) )
    z=b*np.sinh(V) - ( b*( -m +a*np.cosh(V) )*np.sinh(V) ) / ( c * np.cos(U) + a*np.cosh(V) )
    return [-x,y,z]

#def EqCyclide(U,V,m):
 #   x=-c*np.cosh(V)+( (-m+a*np.cosh(V))*(-a*np.cos(U)+c*np.cosh(V)) )/( -c*np.cos(U) + a*np.cosh(V) )
  #  y= ( b*( -m + a*np.cosh(V) ) * np.sin(U) ) / ( -c*np.cos(U) + a* np.cosh(V) )
   # z=b*np.sinh(V) - ( b*( -m +a*np.cosh(V) )*np.sinh(V) ) / ( -c * np.cos(U) + a*np.cosh(V) )
    #return [x,y,z]

    
def NewCyclideBord(m,V):
    if m < c :
        uLim= np.arccos(-m/c)
        U=np.linspace(np.pi-uLim,np.pi+uLim,nbUMax)
    else :
        U=np.linspace(0, 2*np.pi, nbUMax)
    
    return EqCyclide(U,V,m)
    
def NewCyclide(m,vMax):
    if m < c :
        uLim= np.arccos(-m/c)
        if portionCyc[0]>uLim:
            uMin=np.pi-portionCyc[0]
        else :
            uMin=np.pi-uLim
        if portionCyc[1]<np.pi+uLim :
            uMax=portionCyc[1]
        else :
            uMax=np.pi+uLim
        uTemp=np.linspace(uMin, uMax, nbUMax)
    else :
        uTemp=np.linspace(portionCyc[0], portionCyc[1], nbUMax)
    if m > a :
        vMin=np.arccosh(m/a)*np.sign(vMax)
        vTemp=np.linspace(vMin,vMax,nbVMax)
    else :
        vTemp=np.linspace(0,vMax,nbVMax)
    U,V=np.meshgrid(uTemp,vTemp)
    
    return EqCyclide(U,V,m)
    
def NewConeE(sommet):
    theta,t = np.meshgrid(np.linspace(0, 2*np.pi, nbUMax),np.linspace(0, 1, nbVMax))
    
    x=sommet[0]*(1-t)+t*a*np.cos(theta)
    y=sommet[1]*(1-t)+t*b*np.sin(theta)
    z=sommet[2]*(1-t)
    
    return [x,y,z]

def NewConeH(sommet):
    theta,t = np.meshgrid(np.linspace(-np.arccosh(a/c),np.arccosh(a/c),nbUMax),np.linspace(0, 1, nbVMax))
    
    x=sommet[0]*(1-t)-t*c*np.cosh(theta)
    y=sommet[1]*(1-t)
    z=sommet[2]*(1-t)+t*b*np.sinh(theta)
    
    return [x,y,z]

# def NewCone(sommet):
#     theta,t = np.meshgrid(np.linspace(0, 2*np.pi, nbUMax),np.linspace(0, 1, nbVMax))
    
#     x=sommet[0]*(1-t)+t*a*np.cos(theta)
#     y=sommet[1]*(1-t)
#     z=sommet[2]*(1-t)+t*b*np.sin(theta)
    
#     return [x,y,z]


def MonoCouleur(mList,couleur):
    couleurs=[]
    for m in mList :
        couleurs.append((couleur))
    return couleurs
    
def GradCouleur(mList,grad) :
    couleurs=[]
    if len(mList) == 1 :
        return grad[0]
    for i in np.arange(len(mList)):
        mNorm=(mList[i]-mList[0])/(mList[-1]-mList[0])
        m1=mNorm*(len(grad)-1)
        t1=int(m1//1)
        t2=m1%1
        if mList[i]==mList[-1] :
            c=grad[t1]
            couleurs.append(c)
        else :
            c=((1-t2)*grad[t1][0]+t2*grad[t1+1][0],(1-t2)*grad[t1][1]+t2*grad[t1+1][1],(1-t2)*grad[t1][2]+t2*grad[t1+1][2])
            couleurs.append(c)
    return couleurs    

def Ftest(v):
    return(np.tan(ThetaY)*((b*np.sinh(v))-(np.sin(ThetaY)*xc))+(c*np.cosh(v))-xc*np.cos(ThetaY))
    
def dico(inter):
    mid=(inter[0]+inter[1])/2
    interF=inter.copy()
    if Ftest(inter[0])*Ftest(mid)<0 :
        interF[1]=mid
    else :
        interF[0]=mid
    return interF


#=============== Debut du programme ===============#
#----- Couleurs des Cyclides -----#
# couleurCyclides1=MonoCouleur(mList1,(0.2,0.8, 0.8))
couleurCyclides2=GradCouleur(mList2,[(1,0,0),(0,1,0),(0,0,1)])


#----- Calcul des variables intermediaires -----#

nbUMax=50 #nbUMax=50
nbVMax=100 #nbVMax=10
ThetaZ=2*np.pi/n
c = np.sqrt(a**2 - b**2)

aProj=a*np.cos(ThetaY)
xp2=((pow(aProj,4)/pow(b,2))*pow(np.tan(ThetaZ/2),2))/(1+(pow(aProj,2)/pow(b,2))*pow(np.tan(ThetaZ/2),2))
xc=np.sqrt(xp2)+(pow(aProj,2)-xp2)*np.sqrt(1/xp2)


W=[0,20] # 20 au max
for i in range(20): # 20 au max
    W=dico(W) 
vHyperboleMax=W[0] # vHyperboleMax=W[0]
# vHyperboleMax=20
W=[-20,0] # -20 au max
for i in range(20): # 20 au max
    W=dico(W)
vHyperboleMin=W[0]
# vHyperboleMin=-20


#----- Echantillonage des parametres -----# 
# Choix des domaines de definition et creation des listes representant les parametres des courbes et surfaces

u = np.linspace(0, 2*np.pi, 100) #parametre des ellipses
# v1 = np.linspace(0, vHyperboleMax, 5*nbVMax) #parametre des hyperboles 5*nbVMax
v1 = np.linspace(-vHyperboleMax, vHyperboleMax, 5*nbVMax) #parametre des hyperboles 5*nbVMax

# v2 = np.linspace(vHyperboleMin, 0, 5*nbVMax) #parametre des hyperboles 5*nbVMax 

v2 = np.linspace(-vHyperboleMin, vHyperboleMin, 5*nbVMax)
# v2 = np.linspace(-np.pi/2, np.pi/2, 5*nbVMax)
# v2 = np.linspace(-vHyperboleMin, 0, 5*nbVMax)


#----- Calcul des courbes et surfaces -----#

Ellipse=[a * np.cos(u),  b*np.sin(u),np.zeros(len(u))]

Hyperbole1=[c*np.cosh(v1),np.zeros(len(v1)),b*np.sinh(v1)] 
# Hyperbole1=[-c*np.cosh(v1),np.zeros(len(v1)),b*np.sinh(v1)]
Hyperbole2=[-c*np.cosh(v2),np.zeros(len(v2)),b*np.sinh(v2)]


# sommetCone1=[Hyperbole1[0][-1],Hyperbole1[1][-1],Hyperbole1[2][-1]]
# Cone1=NewCone(sommetCone1)

# sommetCone2=[-xc,0,0]

sommetCone1=[Hyperbole2[0][0],Hyperbole2[1][0],Hyperbole2[2][0]] 
Cone1=NewConeE(sommetCone1) 

sommetCone2=[Hyperbole2[0][-1],Hyperbole2[1][-1],Hyperbole2[2][-1]] 
Cone2=NewConeE(sommetCone2) 

# sommetCone3=[Hyperbole2[0][-1],Hyperbole2[1][-1],Hyperbole2[2][-1]]

sommetCone3=[xc,0,0] 
Cone3=NewConeH(sommetCone3)

sommetCone4=[-xc,0,0] 
Cone4=NewConeH(sommetCone4)
 
Cyc1=[]
CycB1=[]
for m in mList1 :
    Cyc1.append(NewCyclide(m,vHyperboleMax))
    CycB1.append(NewCyclideBord(m,vHyperboleMax))
    

# Cyc2=[]
# CycB2=[]
# for m in mList2 :
#     Cyc2.append(NewCyclide(m,vHyperboleMax))
#     CycB2.append(NewCyclideBord(m,vHyperboleMax))

Cyc2Min=[]
CycB2=[]
for m in mList2 :
    Cyc2Min.append(NewCyclide(m,vHyperboleMin))
    CycB2.append(NewCyclideBord(m,vHyperboleMin))
    
Cyc2Max=[]
CycB2=[]
for m in mList2 :
    Cyc2Max.append(NewCyclide(m,vHyperboleMax)) 
    CycB2.append(NewCyclideBord(m,vHyperboleMax))


#----- Rotation des courbes et surfaces -----#

ER=Rot3D(Ellipse,ThetaY,1)
HR1=Rot3D(Hyperbole1,ThetaY,1)
HR2=Rot3D(Hyperbole2,ThetaY,1)
CoR1=Rot3DSurf(Cone1,ThetaY,1)
CoR2=Rot3DSurf(Cone2,ThetaY,1)
CoR3=Rot3DSurf(Cone3,ThetaY,1)
CoR4=Rot3DSurf(Cone4,ThetaY,1)

CycTemp1=[]
for C in Cyc1 :
    CycTemp1.append(Rot3DSurf(C,ThetaY,1))
Cyc1=CycTemp1.copy()

CycTemp2=[]
for C in Cyc2Min :
    CycTemp2.append(Rot3DSurf(C,ThetaY,1))
Cyc2Min=CycTemp2.copy()

CycTemp2=[]
for C in Cyc2Max : 
    CycTemp2.append(Rot3DSurf(C,ThetaY,1))
Cyc2Max=CycTemp2.copy()

CycBTemp1=[]
for CB in CycB1 :
    CycBTemp1.append(Rot3D(CB,ThetaY,1))
CycB1=CycBTemp1.copy()

CycBTemp2=[]
for CB in CycB2 :
    CycBTemp2.append(Rot3D(CB,ThetaY,1))
CycB2=CycBTemp2.copy()


#----- Translation des courbes et surfaces -----#

ER[0]=xc+ER[0]
HR1[0]=xc+HR1[0]
HR2[0]=xc+HR2[0]
CoR1[0]=xc+CoR1[0]
CoR2[0]=xc+CoR2[0]
CoR3[0]=xc+CoR3[0]
CoR4[0]=xc+CoR4[0] 

CycTemp1=[]
for C in Cyc1 :
    CycTemp1.append([C[0]+xc,C[1],C[2]])
Cyc1=CycTemp1.copy()
    
CycBTemp1=[]
for CB in CycB1 :
    CycBTemp1.append([CB[0]+xc,CB[1],CB[2]])
CycB1=CycBTemp1.copy()

CycTemp2=[]
for C in Cyc2Min :
    CycTemp2.append([C[0]+xc,C[1],C[2]])
Cyc2Min=CycTemp2.copy()

CycTemp2=[]
for C in Cyc2Max :
    CycTemp2.append([C[0]+xc,C[1],C[2]])
Cyc2Max=CycTemp2.copy()
    
CycBTemp2=[]    
for CB in CycB2 :
    CycBTemp2.append([CB[0]+xc,CB[1],CB[2]])
CycB2=CycBTemp2.copy()

#----- Duplication du modele et positionnement des duplicatas + configuration du plot  -----#

fig = plt.figure(frameon=False) 
# ax = fig.gca(projection='3d')
ax = fig.add_subplot(projection='3d')

# if nAff==1:
#     cst=Hyperbole1[2][-1]
#     ax.set_xlim(-cst+xc,cst+xc)
#     ax.set_ylim(-cst,cst)
#     ax.set_zlim(-cst, cst)
# else :
#     ax.set_box_aspect((1,1,1)) 
    
    
    
for i in np.arange(0,nAff):
    
    if afficher_ellipse :
        ERZ=Rot3D(ER,i*ThetaZ,2)
        ax.plot(ERZ[0], ERZ[1], ERZ[2],color=couleurEllipse,linewidth=epaisseurEllipse)
    
    if afficher_Hyperbole1 :
        HRZ1=Rot3D(HR1,i*ThetaZ,2)  
        ax.plot(HRZ1[0], HRZ1[1], HRZ1[2],color=couleurHyperbole1,linewidth=epaisseurHyperbole)
    
    if afficher_Hyperbole2 :
        HRZ2=Rot3D(HR2,i*ThetaZ,2)  
        ax.plot(HRZ2[0], HRZ2[1], HRZ2[2],color=couleurHyperbole2,linewidth=epaisseurHyperbole)
        
    # if afficher_Cones1 :
    #     CoRZ1=Rot3DSurf(CoR1,i*ThetaZ,2)
    #     ax.plot_surface(CoRZ1[0], CoRZ1[1], CoRZ1[2],color=couleurCone1,alpha=alphaCone1)
       
    if afficher_Cones2 :
        CoRZ1=Rot3DSurf(CoR1,i*ThetaZ,2)
        CoRZ2=Rot3DSurf(CoR2,i*ThetaZ,2)
        CoRZ3=Rot3DSurf(CoR3,i*ThetaZ,2)
        CoRZ4=Rot3DSurf(CoR4,i*ThetaZ,2)
        # ax.plot_surface(CoRZ1[0], CoRZ1[1], CoRZ1[2],color=couleurCone1,alpha=alphaCone1)
        # ax.plot_surface(CoRZ2[0], CoRZ2[1], CoRZ2[2],color=couleurCone2,alpha=0.5)
        # ax.plot_surface(CoRZ2[0], CoRZ2[1], -CoRZ2[2],color=couleurCone2,alpha=0.5)
        # ax.plot_surface(-CoRZ2[0], -CoRZ2[1], -CoRZ2[2],color=couleurCone2,alpha=alphaCone2)   
        # ax.plot_surface(CoRZ3[0], CoRZ3[1], CoRZ3[2],color=couleurCone2,alpha=0.5) 
        
        ax.plot_surface(CoRZ1[0], CoRZ1[1], CoRZ1[2],color='red',alpha=0.3)
        ax.plot_surface(CoRZ2[0], CoRZ2[1], CoRZ2[2],color='green',alpha=0.3)
        ax.plot_surface(CoRZ3[0], CoRZ3[1], CoRZ3[2],color='blue',alpha=0.4)   
        ax.plot_surface(CoRZ4[0], CoRZ4[1], CoRZ4[2],color='magenta',alpha=0.4)
   
    
   # if afficher_Cones2 :
    #     CoRZ2=Rot3DSurf(CoR2,i*ThetaZ,2)
    #     ax.plot_surface(CoRZ2[0], CoRZ2[1], CoRZ2[2],color=couleurCone2,alpha=alphaCone2)
        
    # if afficher_Cyclides1 :
    #     for j in np.arange(len(Cyc1)):
    #         CZ1=Rot3DSurf(Cyc1[j],i*ThetaZ,2)
    #         ax.plot_surface(CZ1[0],CZ1[1],CZ1[2],color=couleurCyclides1[j],antialiased=False,linewidth=0,alpha=alphaCyclides1)
        
    # if afficher_Cyclides2 :
    #     for j in np.arange(len(Cyc2)):
    #         CZ2=Rot3DSurf(Cyc2[j],i*ThetaZ,2)
    #         ax.plot_surface(CZ2[0],CZ2[1],CZ2[2],color=couleurCyclides2[j],antialiased=False,linewidth=0,alpha=alphaCyclides2)
            
    # if afficher_Cyclides1 :
    #     for j in np.arange(len(Cyc1)):
    #         CZ1=Rot3DSurf(Cyc1[j],i*ThetaZ,2)
    #         ax.plot_surface(CZ1[0],CZ1[1],CZ1[2])  
        
    if afficher_Cyclides2 :
        for j in np.arange(len(Cyc2Min)):
            CZ2Min=Rot3DSurf(Cyc2Min[j],i*ThetaZ,2)
            ax.plot_surface(CZ2Min[0],CZ2Min[1],CZ2Min[2],color=couleurCyclides2[j])
            # ax.plot_surface(-CZ2[0],-CZ2[1],-CZ2[2],color=couleurCyclides2[j])
            
        for j in np.arange(len(Cyc2Max)):
            CZ2Max=Rot3DSurf(Cyc2Max[j],i*ThetaZ,2)
            ax.plot_surface(CZ2Max[0],CZ2Max[1],CZ2Max[2],color=couleurCyclides2[j])
     
    # if afficher_Bord1 : 
    #     for CB in CycB1 :
    #         CBZ1=Rot3D(CB,i*ThetaZ,2)
    #         ax.plot(CBZ1[0],CBZ1[1],CBZ1[2],color=couleurBord1,linewidth=epaisseurBord)    
            
    # if afficher_Bord2 :
    #     for CB in CycB2 :
    #         CBZ2=Rot3D(CB,i*ThetaZ,2)
    #         ax.plot(CBZ2[0],CBZ2[1],CBZ2[2],color=couleurBord2,linewidth=epaisseurBord)

# ax.scatter(0, 0, 0, color='black')
ax.set(xlabel = 'x', ylabel = 'y', zlabel = 'z')
ax.set(xticklabels = [], yticklabels = [], zticklabels = [])
# ax.view_init(3,-83)
# ax.view_init(11,-31)
# ax.view_init(180,-31)
ax.view_init(-177,-90)
# ax.view_init(180,-90)

# ax.legend()

plt.show()    
