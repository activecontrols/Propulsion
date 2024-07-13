function [mdot, b_ratio] = sonicNozzleSizing(P1, P2, T1, d_throat, d_inlet, gamma, Cd, fluid)

% This function computes either the mdot produced by a sonic nozzle devic. 
% To perform calculations, this code uses an
% approach outlined in the Scottmonauts Zucrow document.

% Inputs:
% P1:       Upstream pressure of sonic nozzle, Units: PSIA
% P2:       Donwstream pressure of sonic nozzle. Used to help validate choking, Units: PSIA
% T1:       Upstream temperature of sonic nozzle, Units: F
% d_throat: Diameter of sonic nozzle throat, Units: in.
% d_inlet:  Diameter of sonic nozzle inlet tubing, Units: in.
% gamma:    Specific heat ratio of gas, Unitless
% Cd:       Approximate coefficient of discharge for sonic nozzle, specified by vendor, Unitless
% fluid:    String input representing flow gas. Must be inputted as a valid string for use in REFPROP calculations

% Outputs:
% mdot:     Mass flow rate through nozzle, Units: lbm/s

%% Computations

T1 = (5/9)*(T1 - 32) + 273.15; %Converts T1 to K
P2 = P2 * 6.89476; %Converts P2 to kPa
A_throat = (d_throat/2)^2 * pi * 0.0254^2; % Computes nozzle throat area
b_ratio = d_throat/d_inlet;
P1 = P1 * 6.89476; %P1 in kPA

if P2 <= P1*0.8
    rho = refpropm("D", "T", T1, "P", P1, fluid);
    mdot = Cd*A_throat*(2*rho*P1*1000)^0.5 * 2.20462; %Mdot in lbm/s
    b_ratio = d_throat/d_inlet;
else
    fprintf("\nUpstream nozzle pressure is not high enough to induce flow choking! Please check pressure inputs!\n")
end