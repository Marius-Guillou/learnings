import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

# Définition des paramètres
pi = 3.14159
T = 100  # Poussée (N)
R_val = 0.75  # Rayon (m)
Vrot = 1050  # Rotations par minute
Omega_val = Vrot * 2 * pi / 60  # Vitesse de rotation (rad/s)
Nb = 6  # Nombre de pales
rho = 1.225  # Masse volumique de l'air

# Calcul des grandeurs moyennes
def calc_Vi_P_moy(R, T=100, rho=1.225):
    Vi_moy = np.sqrt(T / (2 * rho * pi * R**2))
    P_moy = np.sqrt(T**3 / (2 * rho * pi * R**2))
    return Vi_moy, P_moy

Vi_moy, P_moy = calc_Vi_P_moy(R_val, T, rho)

# Lois variables
def Theta(r, R):  # Loi de vrillage linéaire
    theta_racine = 0.28
    theta_bout = 0.10
    return theta_racine + (theta_bout - theta_racine) * r / R

def c(r, R):  # Corde variable linéairement
    c_racine = 0.06
    c_bout = 0.04
    return c_racine + (c_bout - c_racine) * r / R

def Cz(alpha):  # Cz réaliste pour NACA 0012
    Cz_max = 1.1
    Cz_lin = 2 * pi * alpha
    return np.clip(Cz_lin, -Cz_max, Cz_max)

# Élément de pale
def Poussee_elmt(r, R, Omega, Vi_moy):
    Beta = np.arccos(Omega * r / np.sqrt((Omega * r)**2 + Vi_moy**2))
    alpha = Theta(r, R) - Beta
    dF = 0.5 * rho * c(r, R) * (Vi_moy**2 + (Omega * r)**2) * Cz(alpha) * np.cos(Beta)
    return dF

def Trainee_elmt(r, R, Omega, Vi_moy):
    Beta = np.arccos(Omega * r / np.sqrt((Omega * r)**2 + Vi_moy**2))
    alpha = Theta(r, R) - Beta
    dT = 0.5 * rho * c(r, R) * (Vi_moy**2 + (Omega * r)**2) * Cz(alpha) * np.sin(Beta)
    return dT

# Intégration
def calc_poussee_totale(R, Omega, Vi_moy):
    r_vals = np.linspace(0.01, R, 500)
    dF_vals = Poussee_elmt(r_vals, R, Omega, Vi_moy)
    F_par_pale = simpson(dF_vals, r_vals)
    return Nb * F_par_pale

def calc_couple_total(R, Omega, Vi_moy):
    r_vals = np.linspace(0.01, R, 500)
    dT_vals = Trainee_elmt(r_vals, R, Omega, Vi_moy)
    dC_vals = r_vals * dT_vals
    C_par_pale = simpson(dC_vals, r_vals)
    return Nb * C_par_pale

F_totale = calc_poussee_totale(R_val, Omega_val, Vi_moy)
C_total = calc_couple_total(R_val, Omega_val, Vi_moy)

# Tracé poussée et trainée élémentaires
x = np.linspace(0.01, R_val, 500)
y = Poussee_elmt(x, R_val, Omega_val, Vi_moy)
z = Trainee_elmt(x, R_val, Omega_val, Vi_moy)

plt.figure(figsize=(10, 6))
plt.plot(x, y, label="Poussée élémentaire", linewidth=2)
plt.plot(x, z, label="Trainée élémentaire", linestyle='--', linewidth=2)
plt.title("Distribution de la poussée et de la trainée le long de la pale")
plt.xlabel("Rayon r (m)")
plt.ylabel("Force élémentaire (N/m)")
plt.grid(True, linestyle=':')
plt.legend()
plt.tight_layout()
plt.show()

# Mapping R-Omega
R_vals = np.linspace(0.3, 1.0, 30)
Omega_vals = np.linspace(30, 400, 5)
R_grid, Omega_grid = np.meshgrid(R_vals, Omega_vals)
F_totale_grid = np.zeros_like(R_grid)

for i in range(R_grid.shape[0]):
    for j in range(R_grid.shape[1]):
        R_test = R_grid[i, j]
        Omega_test = Omega_grid[i, j]
        Vi_moy_test, _ = calc_Vi_P_moy(R_test)
        r_vals = np.linspace(0.01, R_test, 300)
        dF_vals = Poussee_elmt(r_vals, R_test, Omega_test, Vi_moy_test)
        F_totale_grid[i, j] = Nb * simpson(dF_vals, r_vals)

# Tracé de la carte de poussée
plt.figure(figsize=(10, 6))
contourf = plt.contourf(R_grid, Omega_grid, F_totale_grid, levels=50, cmap="viridis")
contour = plt.contour(R_grid, Omega_grid, F_totale_grid, levels=[100], colors='red', linewidths=2)
cbar = plt.colorbar(contourf)
cbar.set_label("Poussée totale (N)")
plt.xlabel("Rayon R (m)")
plt.ylabel("Vitesse de rotation Ω (rad/s)")
plt.title("Carte de la poussée totale selon R et Ω (isocourbe à 100 N)")
plt.grid(True)
plt.tight_layout()
plt.show()

print(F_totale, C_total)
