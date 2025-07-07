#Importations
import numpy as np
import matplotlib.pyplot as plt

#Définition des paramètres
T = 100 #Poussée (N)
R = 0.40 #Rayon (m)
Omega = 100 #Vitesse de rotation (rad/s)
Nb = 6 #Nombre de pales
rho = 1.225 #Masse volumique air
pi = 3.14159

#Définitions des lois variables
def Theta(r): #Theta angle de calage (degrés), on définit la loi de vrillage ici
    return pi/12

def c(r): #c corde de l'élément de pale (mm), on définit sa loi de variation 
    return 0.040

def Cz(alpha): #Cz en fonction de alpha, à partir de la polaire du profil
    return alpha

#Elément de pale
def Poussee_elmt(r):
    Beta = np.arccos(Omega*r/np.sqrt((Omega*r)**2 + Vi_moy**2)) #Angle d'arrivée du flux par rapport au plan de rotation
    alpha = Theta(r) - Beta #Angle d'attaque
    dF = 0.5*rho*c(r)*(Vi_moy**2 + (Omega*r)**2)*Cz(alpha)*np.cos(Beta) #Poussée élémentaire
    return dF 

def Trainee_elmt(r):
    Beta = np.arccos(Omega*r/np.sqrt((Omega*r)**2 + Vi_moy**2)) #Angle d'arrivée du flux par rapport au plan de rotation
    alpha = Theta(r) - Beta #Angle d'attaque
    dT = 0.5*rho*c(r)*(Vi_moy**2 + (Omega*r)**2)*Cz(alpha)*np.sin(Beta) #Poussée élémentaire
    return dT 

#Grandeur moyenne
Vi_moy = np.sqrt(T/(2*rho*pi*R**2)) #Vitesse induite moyenne (m/s)
P_moy = np.sqrt(T**3/(2*rho*pi*R**2)) #Puissance moyenne requise (W)
print("Vi_moy = ", Vi_moy,", Pi_moy = ", P_moy) #Affichage des valeurs

#Tracés

x = np.linspace(0.01, R, 500)
y = Poussee_elmt(x)

plt.plot(x, y, label="poussée élémentaire")
plt.show()