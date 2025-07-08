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

def calc_Vi_P_moy(R, T=100, rho=1.225):
    pi = 3.14159
    Vi_moy = np.sqrt(T / (2 * rho * pi * R**2))
    P_moy = np.sqrt(T**3 / (2 * rho * pi * R**2))
    return Vi_moy, P_moy

# Grandeur moyenne
Vi_moy, P_moy = calc_Vi_P_moy(R, T, rho)
print("Vi_moy = ", Vi_moy, ", P_moy = ", P_moy)


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

#Mapage des couples (R, omega)

# Redéfinir Poussee_elmt pour dépendre de R, Omega, Vi_moy
def Poussee_elmt_param(r, R_test, Omega_test, Vi_moy):
    Beta = np.arccos(Omega_test * r / np.sqrt((Omega_test * r)**2 + Vi_moy**2))
    alpha = Theta(r) - Beta
    dF = 0.5 * rho * c(r) * (Vi_moy**2 + (Omega_test * r)**2) * Cz(alpha) * np.cos(Beta)
    return dF

# Grilles de R et Omega
R_vals = np.linspace(0.3, 1.0, 30)          # m
Omega_vals = np.linspace(30, 400, 5)    # rad/s
R_grid, Omega_grid = np.meshgrid(R_vals, Omega_vals)
F_totale_grid = np.zeros_like(R_grid)

# Remplissage de la grille avec ta fonction
for i in range(R_grid.shape[0]):
    for j in range(R_grid.shape[1]):
        R_test = R_grid[i, j]
        Omega_test = Omega_grid[i, j]
        Vi_moy = np.sqrt(100 / (2 * rho * pi * R_test**2))  # hypothèse constante

        r_vals = np.linspace(0.01, R_test, 300)
        dF_vals = Poussee_elmt_param(r_vals, R_test, Omega_test, Vi_moy)
        F_par_pale = simpson(dF_vals, r_vals)
        F_totale = Nb * F_par_pale
        F_totale_grid[i, j] = F_totale

# Tracé de la carte de poussée
plt.figure(figsize=(10, 6))
contourf = plt.contourf(R_grid, Omega_grid, F_totale_grid, levels=50, cmap="viridis")
contour = plt.contour(R_grid, Omega_grid, F_totale_grid, levels=[100], colors='red', linewidths=2)

cbar = plt.colorbar(contourf)
cbar.set_label("Poussée totale (N)")

plt.xlabel("Rayon R (m)")
plt.ylabel("Vitesse de rotation Ω (rad/s)")
plt.title("Carte de la poussée totale selon R et Ω (avec isocourbe à 100 N)")
plt.grid(True)
plt.tight_layout()
plt.show()