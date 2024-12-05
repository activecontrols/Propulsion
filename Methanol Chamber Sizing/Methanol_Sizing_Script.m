% Tadpole Methanol Sizing Script
% Adam Grendys
% 11/21/2024

%% ASSUMPTIONS
% Quasi 1-D 
% Isentropic
% Equilibrium
% Stagnation Temperature = Chamber Temperature (Tc = To)
% Stagnation Pressure = Chamber Pressure (Pc = Po)
% 80% Bell Contour
% Theta_n = ~23 degrees (From NASA SP-8120)
% Theta_e = ~12.5 degrees (From NASA SP-8120)

clear
clc
close all

%% INTILIATZATION
% Inputs
Po = 250; % Chamber Pressure / Stagnation Pressure (psia)
Ft = 550; % Engine Thrust (lbf)
L_star = 45; % Characteristic Lenght (in)
Dc = 3.75; % Chmaber Diameter, Consistency with current tapdole infrastructure
Patm = 14.7; % Atmoshperic Pressure (psi)
Pe = 17; 

% Data from NASA SP-8120 for 80% Bell Contour
Rao_inlet = 23; % Parabolic Start angle (Throat End)
Rao_exit = 12.5; % Parabolic End angle (Exit)
Convergence_angle = 37.5; % Convergence Angle (input)


% Data from CEA Run 
OF = 1.05;
g_o = 32.174049; % ft/sec^2
R_prime = 1544; % universal gas constant (ft-lbf/lb-mol-◦R)
ER = 2.8936; % Expansion Ratio, Area of Exit / Area of Throat

% CEA Vectors Parameters in order of [Chamber, Throat, Exit]
C_star = [NaN, 5551.7, 5551.7]; % Characteristic Velocity FT/SEC 
Cf = [NaN, 0.6777, 1.3480]; % Thrust Coeffcient
Ivac = [NaN, 214.3, 266.6]; % LB-SEC/LB 
Isp = [NaN, 116.9, 232.6]; % Specific Impulse LB-SEC/LB
gamma = [1.1980, 1.2012, 1.2181]; % Specific heat ratio
Me = [NaN, 1.0, 2.372]; % Gas Mach Number
M = [21.259, 21.259, 21.259]; % Gas Molecular weight
Son_vel = [1201.0, 1146.8, 961.8] .* 3.281; % Sonic Velocity (speed of sound) (ft/s)
T = [3078.83, 2799.18, 1941.88] .* (9/5); % Gas Temperature (R)
P = [17.011, 9.5973, 1.1568] .* 14.696; % Pressure (psia) 

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
%Lcon = ((Dc - Dt) / 2) / tand(Convergence_angle); % Length of Converging Section (in)
%Ldiv = ((De - Dt) / 2) / tand(Rao_inlet); % Length of Diverging Section (in)
Ln = 0.8 * ((Dt / 2) * (sqrt(CR) - 1) + (1.5 * Dt / 2) * (secd(Rao_inlet) - 1)) / (tand(Rao_inlet));
disp(Ln)


Surface_Area = 2 * Lc * sqrt(pi * CR * At) + cscd(Convergence_angle) * (CR - 1) * At; % Chamber and Converging (in^2)

% Plotting Calculations
gamma_y = mean(gamma);


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


% % Chamber
% x = linspace(0, Lc, 5);
% Cy = [x, Dc / 2];
% % Converging Radius
% x = linspace(Lc, Lcon / 2, 5);
% R_con = 
% Cy = {Cy, };
% % Throat Inlet Radius
% x = linspace(Lc + Lcon / 2, Lc + Lcon, 5);
% R_inlet = 1.5 * Dt/2;
% Cy = {Cy, R_inlet};
% % Throat Exit Radius
% x = linespace(Lc + Lcon, Lc + Lcon + Ldiv / 2, 5);
% R_exit = 0.383 * Dt /2;
% Cy = {Cy, R_exit};
% % Exit Parabola
% x = linespcae(Lc + Lcon + Ldiv/2, Lc + Lcon + Ldiv, 5);
% P_exit = ;
% Cy = {Cy, P_exit};



