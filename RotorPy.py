#Importations
import math 

#Définition des paramètres
T = 100 #Poussée (N)
R = 0 #Rayon (m)
Vrot = 0 #Vitesse de rotation (tr/min)
Nb = 6 #Nombre de pales
rho = 1.225 #Masse volumique air
pi = 3.14159

#Définitions des lois variables
def Theta(r, omega): #Theta angle de calage (degrés), on définit la loi de vrillage ici
    return 8

def c(r, omega): #c corde de l'élément de pale (mm), on définit sa loi de variation 
    return 40

#Elément de pale


#Grandeur moyenne
def Vi_moy(T, R, rho):
    return sqrt(T/(2*rho*pi*R**2))

print(Vi_moy(T, R, rho))