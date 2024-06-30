# Bolt sizing code
# Adam Grendys
# 2023-10-14

# THIS PROGRAM IS A WORK IN PROGRESS! USE AT YOUR OWN RISK.

from math import *
bolt_count = 4 # driving parameter

d_nom = 0.164 # inches
TPI = 32 # threads per inch
preload_percent = 0.765 # 0.9 * 0.85
tensile_yield_strength = 170000 # start with 150 kpsi

A_t = (pi / 4) * (d_nom - (0.9743 / TPI))**2 # tensile stress area

F_PL = preload_percent * tensile_yield_strength * A_t # preload force
print("Preload Force: " + str(round(F_PL, 4)))

# parameters for torque (N/A, including equation)
#d_minor = d_nom - (1.08253175/TPI) # minor diameter
#r_t = (d_nom + d_minor) / 4 # mean thread radius
#r_c = 0.625 * d_nom # mean collar radius
#lambda = arctan(1 / (2*pi*TPI)) # lead angle
#((r_t / d_nom) * (tan(lambda) + f_t*sec(alpha)) / (1 - f_t * tan(lambda)*sec(alpha))) + (f_c * r_c / d_nom)

K_T = 0.2 # value for zinc-plated bolts, Shigley

T_PL = K_T * d_nom * F_PL # preload torque
print("Preload torque (lbf-ft): %.4f\n" %(T_PL / 12))

# calculate overall stress

chamber_pressure = 250 # psi, assuming MAWP
chamber_diameter = 3.75 # inches
chamber_length = 5.205 # inches
manifold_pressure = 342 #psi
chamber_cross_section = pi * (chamber_diameter / 2)**2 #chamber cross sectional area

F_app = (chamber_pressure) * 0.35**2 * pi #force applied, assuming manifold and chamber area are eqaul
F_sh = 0 #chamber_pressure * pi * chamber_diameter * chamber_length #shear force
M_b = 0 # assume no bending moment

A_s = pi / 4 * d_nom**2 # shear area, equal to nominal

sigma_PL = F_PL / (A_t * bolt_count) # preload stress
sigma_app = F_app / (A_t * bolt_count) # applied tensile stress
tau_sh = F_sh / (A_s * bolt_count) # shear stress
sigma_bnd = 32 * M_b / (pi * d_nom**3) # bending stress

load_factor = 2.4 # limited application of safety factor
sigma_VM = load_factor*(1/2*((-sigma_PL + (sigma_app + sigma_bnd)))**2 + 3*((tau_sh)**2))**(1/2) # von Mises stress

print("Shear Stress: " + str(round(tau_sh,4)))

print("Applied Tensile Stress: " + str(round(sigma_app, 4)))

print("Von Mises stress (psi): %.4f\n" %sigma_VM)