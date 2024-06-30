# Bolt sizing code
# Alejandro Diaz Contreras and Adam Grendys
# 2023-10-14
#
# This code can be used to size bolts for a pressure cap.
# In this case, it is used to size the bolts for the injector flange
#
# Revised calcs with:
# https://mechanicalc.com/reference/bolted-joint-analysis


from math import *
#import math

source = '' #Do not change values here!!!! Scroll down to "d" section below (marked with a comment)
bolt_count = 0
d_nom = 0
tpi = 0
preload_percent = 0
tensile_yield_strength = 0
chambers = 0
chamber_press = 0
chamber_diam = 0
F_app = 0

while(not(source == 'u' or source == 'd')):
    source = input('Enter "u" for user defined values or "d" for defaults: ')

    if(source == 'u'):
        bolt_count = int(input('Enter the number of bolts: '))
        d_nom = float(input('Enter nominal diameter (in): '))
        tpi = float(input('Enter Threads Per Inch: '))
        preload_percent = float(input('Enter preload percent'))
        tensile_yield_strength = float(input('Enter bolt yield stress (psi): '))
        
        chamber_press = float(input(f'Enter the chamber pressure (psi): '))
        chamber_diam = float(input(f'Enter the chamber diameter (in): '))

        F_sh = float(input('Enter shear force (lbf): '))
        M_b = float(input('Enter bending moment (lb*ft): '))

    elif(source == 'd'):  # default parameters, edit these
        bolt_count = 8 # driving parameter
        d_nom = 0.2175 # inches
        tpi = 20 # threads per inch
        preload_percent = 0.765 # 0.9 * 0.85
        tensile_yield_strength = 120000 # start with 150 kpsi

        chamber_press = 342 # = 342 psi
        chamber_diam = (3.75) # = 1.2307 * 2 inches

        F_sh = 0 #assuming zero shear force, otherwise --> chamber_pressure * pi * chamber_diameter * chamber_length
        M_b = 0 # assume no bending moment

    else:
        print('Command not recognized, try again.')


generateText = input('Would you like to generate a text file with the outputs? Enter y/n: ')



#Calculate useful values
A_t = (pi / 4) * (d_nom - 0.9743 / tpi)**2 # tensile stress area

F_preload = preload_percent * tensile_yield_strength * A_t * bolt_count # preload force

F_app = chamber_press * pi * ((chamber_diam / 2) ** 2)


# parameters for torque (N/A, including equation)
#d_minor = d_nom - (1.08253175/tpi) # minor diameter
#r_t = (d_nom + d_minor) / 4 # mean thread radius
#r_c = 0.625 * d_nom # mean collar radius
#lambda = arctan(1 / (2*pi*tpi)) # lead angle
#((r_t / d_nom) * (tan(lambda) + f_t*sec(alpha)) / (1 - f_t * tan(lambda)*sec(alpha))) + (f_c * r_c / d_nom)

K_T = 0.2 # value for zinc-plated bolts, Shigley

T_PL = K_T * d_nom * (F_preload / bolt_count) # preload torque


print("Preload torque (lbf-ft): %.4f\n" %(T_PL / 12 ))


"""
############################
# Calculate overall stress # 
############################


for i in range(chambers):
    chamber_cross_sections.append(pi * (chamber_diams[i] / 2)**2) #chamber cross sectional area
    F_app += chamber_press[i] * chamber_cross_sections[i] #This method has a good chance of being wrong...

A_s = pi / 4 * d_nom**2 # shear area, equal to nominal

sigma_PL = F_preload / (A_t * bolt_count) # preload stress
sigma_app = F_app / (A_t * bolt_count) # applied tensile stress
tau_sh = F_sh / (A_s * bolt_count) # shear stress
sigma_bnd = 32 * M_b / (pi * d_nom**3) # bending stress


sigma_VM = (((sigma_PL + (sigma_app + sigma_bnd)))**2 + 3*((tau_sh)**2))**(1/2) # von Mises stress 
#Revise Von mises, im not sure if preload stress should be added like that

"""


if(F_app > F_preload):
    print('Joint failed!\n')
    print(f'Preload force: {F_preload:.2f} lbf')
    print(f'Applied force: {F_app:.2f} lbf')
else:
    print(f'Safety factor: {(F_preload/F_app):.2f}\n')
    print(f'Preload force: {F_preload:.2f} lbf')
    print(f'Applied force: {F_app:.2f} lbf')





"""
print("Shear Stress: " + str(round(tau_sh,4)))

print("Applied Tensile Stress: " + str(round(sigma_app, 4)))

print("Von Mises stress (psi): %.4f\n" %sigma_VM)
"""
