% SuperTadpole Chamber Sizing
% Adam Grendys
% 6/28/2025

%% ASSUMPTIONS
% Quasi 1-D 
% Isentropic
% Frozen 
% Stagnation Temperature = Chamber Temperature (Tc = To)
% Stagnation Pressure = Chamber Pressure (Pc = Po)
% Conical Nozzle

clear
clc
close all

%% INTILIATZATION
% Inputs
Po = 250; % Chamber Pressure / Stagnation Pressure (psia)
Ft = 550; % Engine Thrust (lbf)
L_star = 50; % Characteristic Lenght (in)
Dc = 3.75; % Chamber Diameter (in), Consistency with Tapdole Rev1 design
Patm = 14.7; % Atmoshperic Pressure (psia)
Pe = 17; % Exit Pressure (psia)
Convergence_angle = 37.5; % Convergence Angle (input)
Divergence_angle = 15; % Conical nozzel divergence half angle (same as Tadpole Rev1)

% Data from CEA Run 
OF = 1.2;
g_o = 32.174049; % ft/sec^2
R_prime = 1544; % universal gas constant (ft-lbf/lb-mol-◦R)
ER = 2.7704; % Expansion Ratio, Area of Exit / Area of Throat

% CEA Vectors Parameters in order of [Chamber, Throat, Exit]
C_star = [NaN, 6068.8, 6068.8]; % Characteristic Velocity FT/SEC 
Cf = [NaN, 0.6902, 1.3415]; % Thrust Coeffcient
Ivac = [NaN, 235.2, 288.6]; % LB-SEC/LB 
Isp = [NaN, 130.2, 253.0]; % Specific Impulse LB-SEC/LB
gamma = [1.2355, 1.2396, 1.2602]; % Specific heat ratio
Mach = [NaN, 1.0, 2.377]; % Gas Mach Number
M = [18.321, 18.321, 18.321]; % Gas Molecular weight
Son_vel = [1348.3, 1276.7, 1043.9] .* 3.281; % Sonic Velocity (speed of sound) (ft/s)
T = [3242.15, 2897.55, 1905.42] .* (9/5); % Gas Temperature (R)
P = [17.011, 9.4720, 1.1568] .* 14.696; % Pressure (psia) 

%% CALCULATIONS
R = R_prime / M(1);

% Effcicency Additons
eta_cf = 0.95;
eta_cstar = 0.92;

C_star = C_star .* eta_cstar;
Cf = Cf + (Pe / P(1) - Patm / P(1)) * ER;
Cf = Cf .* eta_cf;
Isp(3) = C_star(3) * Cf(3) / g_o;

% General
Ve = Isp(3) * g_o;
m_dot = Ft / Ve;
m_dot_fuel = m_dot / (1 + OF);
m_dot_ox = m_dot - m_dot_fuel;

% Geometry Calculations
At = m_dot * C_star(2) / P(1); % Throat Area (in^2)
Dt = sqrt(4 * At / pi); % Throat Diameter (in)
Vc = L_star * At; % Chamber Volume (in^3)
Ac = pi * Dc^2 / 4; % Chamber Area (in^2)
Lc = Vc / (1.1 * Ac); % Chamber Length (in)
CR = Ac / At; % Convergence Ratio
Ae = At * ER; % Exit Area (in^2)
De = sqrt(4 * Ae / pi); % Exit Diameter (in^2)
Lcon = (Dc - Dt) / (2 * tand(Convergence_angle)); % Length of Converging Section (in)
Ldiv = (De - Dt) / (2 * tand(Divergence_angle)); % Length of Diverging Section (in)

%% FORMATTED OUTPUT
fprintf("OF: %.2f", OF);
fprintf("\nIsp: %.3f s", Isp(3));
fprintf("\nExpansion Ratio: %.3f ", ER);
fprintf("\nConvergence Ratio: %.3f \n", CR);

fprintf("\nChamber Diameter: %.3f inches", Dc);
fprintf("\nThroat Diameter: %.3f inches", Dt);
fprintf("\nExit Diameter: %.3f inches\n", De);

fprintf("\nChamber Length: %.3f inches", Lc);
fprintf("\nConverging Length: %.3f inches", Lcon);
fprintf("\nDiverging Length: %.3f inches", Ldiv);
fprintf("\nTotal Length: %.3f inches\n", Lc + Lcon + Ldiv);

fprintf("\nTotal mdot: %.3f lbs/s", m_dot * g_o);
fprintf("\nFuel mdot: %.3f lbs/s", m_dot_fuel * g_o);
fprintf("\nOx mdot: %.3f lbs/s", m_dot_ox * g_o);


