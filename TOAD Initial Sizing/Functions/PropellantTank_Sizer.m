function [tank_mass] = PropellantTank_Sizer(prop_mass, OF)
% TOAD Tank Sizer
% Author: Jacob Metcalf
% Reviewer: Adam Grendys
% Last Edited: 2/16/2025

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing
% Ellipsoildal Tank Endcaps

%%%% ADJUSTABLE VARIABLES FOR TANK SIZING %%%%%%%%%%%
% Bechtel Limits us to max external diameter of 16".
wall_thickness = 1/4; % Tank wall thickness (in)
tank_OD = 12; % Tank outer diameters (in)
tank_pressure = 550; % Internal tank pressure (psia)
safety_factor = 2; % Tank design safety factor
weld_factor = 0.8; % Additional safety factor due to tank welds

%% INITIALZATION
s_y = (276*10^6) / (6894.76); % 6061-T6 Yield Stress (psi)
s_u = (310*10^6) / (6894.76); % 6061-T6 Ultimate Tensile Strength (psi)
density_ipa = 49.06838; % lbm/ft^3 @ 293.15K
density_lox = 71.798; % lbm/ft^3 @ 90K, 550psi
density_tanks = 168.5555; % lbm/ft^3 (Aluminium 6061-T6 @ 293.15K)
ullage_vol_coef = 1.1; % 10% Ullage volume
fuel_reserves_coef = 1.5; % 50% Fuel excess Fuel and Oxidizer mass

%% CALCULATIONS
% Propellant Calcs
ipa_mass = prop_mass / (1 + OF);
lox_mass = prop_mass - ipa_mass;

ipa_mass = fuel_reserves_coef * ipa_mass; % Total IPA mass (lbm)
lox_mass = fuel_reserves_coef * lox_mass; % Total LOX mass (lbm)
ipa_vol = ullage_vol_coef * (ipa_mass / density_ipa); % Volume (ft^3)
lox_vol = ullage_vol_coef * (lox_mass / density_lox); % Volume (ft^3)

% Stress Calcs
tank_radius = tank_OD / 2 - wall_thickness; % Tank INNER radii (in)

s_h = (tank_pressure * tank_radius) / wall_thickness; % Hoop stress (psi) (same for IPA and LOx)
s_a = (tank_pressure * tank_radius) / (2 * wall_thickness); % Axial stress (psi) (same for IPA and LOx)
s_v = sqrt(s_h^2 + s_a^2 - s_h * s_a); % Von Mises stress (psi) (same for IPA and LOx)

M_y = s_y / ((safety_factor + weld_factor) * s_h) - 1; % Margin to yield
M_u = s_u / ((safety_factor + weld_factor) * s_h) - 1; % Margin to ultimate

% Tank Dimensions
bulkhead_inner_vol = 4/3 * pi * tank_radius^2 * sqrt(tank_radius); % Ellipsoildal tank bulkhead volume (ft^3)
h_ipa = (ipa_vol * 12^3 - bulkhead_inner_vol) / (pi * tank_radius^2); % Cylindrical section height (in)
h_lox = (lox_vol * 12^3 - bulkhead_inner_vol) / (pi * tank_radius^2); % Cylindrical section height (in)

% Tank Mass
if h_lox < 0 || h_ipa < 0
    fprintf("\n\nTANK IS SHORTER THAN IT IS WIDE\n")
    output = 0;
else
    output = 1;

    bulkhead_vol = 4/3 * pi * (tank_radius + wall_thickness)^2 * sqrt(tank_radius + wall_thickness) - 4/3 * pi * tank_radius^2 * sqrt(tank_radius); % Tank material volume of bulkheads (in^3)
    clyinder_vol_ipa = pi * h_ipa * ((tank_radius + wall_thickness)^2 - (tank_radius)^2); % IPA tank material volume of cylinder section (in^3)
    clyinder_vol_lox = pi * h_lox * ((tank_radius + wall_thickness)^2 - (tank_radius)^2); % LOx tank material volume of cylinder section (in^3)

    ipa_total_tank = (bulkhead_vol + clyinder_vol_ipa) / 12^3; % IPA tank material volume (ft^3)
    lox_total_tank = (bulkhead_vol + clyinder_vol_lox) / 12^3; % LOx tank material volume (ft^3)
    ipa_tank_mass = ipa_total_tank * density_tanks; % IPA dry tank mass (lbm)
    lox_tank_mass = lox_total_tank * density_tanks; % LOx dry tank mass (lbm)

    tank_mass = ipa_tank_mass + lox_tank_mass;
end

%% FORMATTED OUTPUT
if output
    fprintf("\n\nTANK DIMENSIONS\n")
    fprintf("Wall Thickness: %.3f in\n", wall_thickness)
    fprintf("Tank Radii: %.3f in\n", tank_radius)
    fprintf("IPA Height: %.3f in\n", h_ipa + 2 * sqrt(tank_radius + wall_thickness))
    fprintf("LOx Height: %.3f in\n", h_lox + 2 * sqrt(tank_radius + wall_thickness))

    fprintf("\nTANK MASSES\n")
    fprintf("IPA Tank Volume: %.3f ft^3\n", ipa_vol)
    fprintf("LOx Tank Volume: %.3f ft^3\n", lox_vol)
    fprintf("IPA Tank Mass: %.3f lbs\n", ipa_tank_mass)
    fprintf("LOx Tank Mass: %.3f lbs\n", lox_tank_mass)    
    
    fprintf("\nTANK STRESS\n")
    fprintf("Safety Factor: %.2f\n", safety_factor)
    fprintf("Weld Factor: %.2f\n", weld_factor)
    fprintf("Von Mises: %.3f psi\n", s_v)
    fprintf("Margin to Yield: %.3f\n", M_y)
    fprintf("Margin to Ultimate: %.3f\n", M_u)
end
