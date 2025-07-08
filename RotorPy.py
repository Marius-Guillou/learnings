#Importations
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

#Définition des paramètres
pi = 3.14159
T = 100 #Poussée (N)
R = 0.75 #Rayon (m)
Vrot = 1050 #Rotations par minutes
Omega = Vrot*2*pi/60 #Vitesse de rotation (rad/s)
Nb = 6 #Nombre de pales
rho = 1.225 #Masse volumique air

#Grandeur moyenne
Vi_moy = np.sqrt(T/(2*rho*pi*R**2)) #Vitesse induite moyenne (m/s)
P_moy = np.sqrt(T**3/(2*rho*pi*R**2)) #Puissance moyenne requise (W)
print("Vi_moy = ", Vi_moy,", Pi_moy = ", P_moy) #Affichage des valeurs

#Définitions des lois variables
def Theta(r): #Theta angle de calage (degrés), on définit la loi de vrillage ici
    theta_racine = 0.28  # en radians
    theta_bout = 0.10  # ≈ 16° en radians
    return theta_racine + (theta_bout - theta_racine) * r / R

def c(r): #c corde de l'élément de pale (mm), on définit sa loi de variation 
    c_racine = 0.06  # 6 cm
    c_bout = 0.04  # 4 cm
    return c_racine + (c_bout - c_racine) * r / R

def Cz(alpha): #Cz en fonction de alpha, à partir de la polaire du profil
    return 2*pi*alpha

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

#Integration
def calc_poussee_totale(): # Fonction : poussée totale sur l'ensemble du rotor
    r_vals = np.linspace(0.01, R, 500)  # évite r=0
    dF_vals = Poussee_elmt(r_vals)
    F_par_pale = simpson(dF_vals, r_vals)  # intégrale sur une pale
    F_totale = Nb * F_par_pale
    return F_totale

def calc_couple_total(): # Fonction : couple total dû à la trainée
    r_vals = np.linspace(0.01, R, 500)
    dT_vals = Trainee_elmt(r_vals)
    dC_vals = r_vals * dT_vals  # moment élémentaire = r * dT
    C_par_pale = simpson(dC_vals, r_vals)
    C_total = Nb * C_par_pale
    return C_total

F_totale = calc_poussee_totale()
C_total = calc_couple_total()

print(f"Poussée totale estimée : {F_totale:.2f} N")
print(f"Couple total estimé : {C_total:.3f} N·m")

#Tracés

x = np.linspace(0.01, R, 500)
y = Poussee_elmt(x)
z = Trainee_elmt(x)

plt.figure(figsize=(10, 6))  # Taille plus large

plt.plot(x, y, label="Poussée élémentaire", linewidth=2)
plt.plot(x, z, label="Trainée élémentaire", linestyle='--', linewidth=2)

plt.title("Distribution de la poussée et de la trainée le long de la pale", fontsize=14)
plt.xlabel("Rayon r (m)", fontsize=12)
plt.ylabel("Force élémentaire (N/m)", fontsize=12)

plt.grid(True, which='both', linestyle=':', linewidth=0.5)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()

plt.show()


