function [total_mass] = tank_sizer(mdot, OF, flight_time)
% TOAD Tank Sizer
% Author: Jacob Metcalf
% Reviewer: Adam Grendys
% Last Edited: 2/16/2025

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing

%%%% ADJUSTABLE VARIABLES FOR TANK SIZING %%%%%%%%%%%
max_radius_ipa = 0.66; % ft
min_radius_ipa = 0.25; % ft
max_height_ipa = 5; % ft

max_radius_lox = 0.66; % ft
min_radius_lox = 0.25; % ft
max_height_lox = 5; % ft

tank_wall_thickness = 1/8 / 12; % Tank wall thickness (1/8" to ft)

%% INITIALZATION
density_ipa = 49.0684; % lbm/ft^3 @ 293.15K
density_lox = 71.23030301738; % lbm/ft^3 @ 90K
density_tanks = 168.5555; % lbm/ft^3 (Aluminium 6061-T6 @ 293.15K)
ullage_vol_coef = 1.1; % 10% Ullage volume
fuel_reserves_coef = 1.1; % 10% Fuel excess Fuel and Oxidizer mass

%% CALCULATIONS
mdot_IPA = mdot / (OF + 1);
mdot_LOX = mdot - mdot_IPA;

ipa_mass = fuel_reserves_coef * (mdot_IPA * flight_time); % Total IPA mass (lbm)
ipa_vol = ullage_vol_coef * (ipa_mass / density_ipa); % Volume (ft^3)
lox_mass = fuel_reserves_coef * (mdot_LOX * flight_time); % Total LOX mass (lbm)
lox_vol = ullage_vol_coef * (lox_mass / density_lox); % Volume (ft^3)

%% FORMATTED OUTPUT
fprintf("\nIPA Tank Volume: %.4f ft^3\n", ipa_vol);
fprintf("IPA mass: %.02f lbm\n", ipa_mass);
fprintf("LOX Tank Volume: %.4f ft^3\n", lox_vol);
fprintf("LOX mass: %.02f lbm", lox_mass);

% Define Radius Range
ipa_radius = linspace(min_radius_ipa, max_radius_ipa, 100);
lox_radius = linspace(min_radius_lox, max_radius_lox, 100);

% Compute Heights
height_ipa = ipa_vol ./ (pi * ipa_radius.^2);
height_lox = lox_vol ./ (pi * lox_radius.^2);

cap_height_ipa = ipa_radius/sqrt(2);
cap_height_lox = lox_radius/sqrt(2);
new_cyl_ipa_height = (ipa_vol) ./ ((pi * ipa_radius.^2)) - (4/3) * cap_height_ipa;
new_cyl_lox_height = (lox_vol) ./ ((pi * lox_radius.^2)) - (4/3) * cap_height_lox;
new_height_ipa = new_cyl_ipa_height + (2 * cap_height_ipa);
new_height_lox = new_cyl_lox_height + (2 * cap_height_lox);

% Ensure heights don’t exceed max_height
valid_indices_ipa = new_height_ipa <= max_height_ipa; % Logical mask for valid region
radius_valid_ipa = ipa_radius(valid_indices_ipa); % Keep only valid radii

valid_indices_lox = new_height_lox <= max_height_lox; % Logical mask for valid region
radius_valid_lox = lox_radius(valid_indices_lox); % Keep only valid radii

% Minimum Tank Mass Calculations
min_IPA_tank_vol = min(pi * ((radius_valid_ipa + 2 * tank_wall_thickness).^2 - (radius_valid_ipa.^2)));
min_LOX_tank_vol = min(pi * ((radius_valid_lox + 2 * tank_wall_thickness).^2 - (radius_valid_lox.^2)));
tank_mass_IPA = density_tanks * min_IPA_tank_vol;
tank_mass_LOX = density_tanks * min_LOX_tank_vol;

total_mass = tank_mass_LOX + lox_mass + tank_mass_IPA + ipa_mass;
