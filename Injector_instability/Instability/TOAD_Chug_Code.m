%%%
% TOAD Chug Instability Code
% Using crocco time-lag method as outlined in NASA SP-194 
% Written by Jan Ayala
% Last Modified 5/1/2024
%%%
clear
close all
time_lag = 0.001;
chamber_res_time = 0.001;
mass_flow = 2.2;
cstar = 4000;
dcdr = 1;
OF = 1.2;
pc = 250;
fuel_dp = 0.30 * pc;
ox_dp = 0.30 * pc;
rho_fuel = 800;
line_volume = 0.05;
fuel_man_volume = 0.01;
bulk_fuel = 100;
fuel_manifold_res_time = 0.1;
fuel_area = 0.0001;
ox_area = 0.0001;
fuel_vel = 20;
ox_vel = 20;
length_line = 1;
area_line = 0.02;
line_vel = 20;
a_L = 900;
s = tf([1 0], 1);
%Fuel Injector
fuel_admittance = (fuel_vel / a_L) ^ 2 * s * rho_fuel / bulk_fuel * (line_volume + fuel_man_volume) * pc/mass_flow;
%Ox Injector
ox_admittance = (fuel_vel / a_L) ^ 2 * s * rho_fuel / bulk_fuel * (line_volume + fuel_man_volume) * pc/mass_flow;
%Stability Calcs
stab_equation = exp(-s*time_lag) / (1 + s*chamber_res_time) * (((1 + (1 + OF) / cstar) * dcdr) * ox_admittance ...
    + ((1 - OF * (1 + OF) / cstar) * dcdr) * fuel_admittance);
stab_equation_pade = tf(pade(stab_equation, 4));
lambda = zero(stab_equation_pade + 1);
nyquist(stab_equation)
figure
nyquist(stab_equation_pade)
