function [total_mass, prop_mass] = tank_sizer(mdot, OF, flight_time)
% TOAD Tank Sizer
% Author: Jacob Metcalf
% Reviewer: Adam Grendys
% Last Edited: 2/16/2025

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing

%%%% ADJUSTABLE VARIABLES FOR TANK SIZING %%%%%%%%%%%
% Bechtel Limits us to max external diameter of 16".
tank_wall_thickness = .125 / 12; % Tank wall thickness (1/2" to ft)

max_radius_ipa = 8/12 - tank_wall_thickness; % ft
min_radius_ipa = 1/12; % ft
max_height_ipa = 6; % ft

max_radius_lox = 8/12 - tank_wall_thickness; % ft
min_radius_lox = 1/12; % ft
max_height_lox = 6; % ft

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

% Define Radius and Height Ranges
ipa_radius = linspace(min_radius_ipa, max_radius_ipa, 100);
lox_radius = linspace(min_radius_lox, max_radius_lox, 100);
cyl_height_ipa = (ipa_vol - 4/3 * pi * ipa_radius.^3) ./ (pi * ipa_radius.^2); % Spherical Endcaps
cyl_height_lox = (lox_vol - 4/3 * pi * lox_radius.^3) ./ (pi * lox_radius.^2); % Spherical Endcaps

% Ensure heights don’t exceed max_height
valid_indices_ipa = (cyl_height_ipa + ipa_radius + tank_wall_thickness) <= max_height_ipa; % Logical mask for valid region
ipa_radius = ipa_radius(valid_indices_ipa); % Keep only valid radii
cyl_height_ipa = cyl_height_ipa(valid_indices_ipa);

valid_indices_lox = (cyl_height_lox + lox_radius + tank_wall_thickness) <= max_height_lox; % Logical mask for valid region
lox_radius = lox_radius(valid_indices_lox); % Keep only valid radii
cyl_height_lox = cyl_height_lox(valid_indices_lox);

% Minimum Tank Mass Calculations
[min_IPA_tank_vol, ipa_index] = min(pi * ((ipa_radius + tank_wall_thickness).^2 - ipa_radius.^2) .* cyl_height_ipa + 4/3  * pi * ((ipa_radius + tank_wall_thickness).^3 - ipa_radius.^3));
[min_LOX_tank_vol, lox_index] = min(pi * ((lox_radius + tank_wall_thickness).^2 - lox_radius.^2) .* cyl_height_lox + 4/3  * pi * ((lox_radius + tank_wall_thickness).^3 - lox_radius.^3));

tank_mass_IPA = density_tanks * min_IPA_tank_vol;
tank_mass_LOX = density_tanks * min_LOX_tank_vol;

prop_mass = lox_mass + ipa_mass;
total_mass = tank_mass_LOX + lox_mass + tank_mass_IPA + ipa_mass;

% Calculating Total Heights
total_height_ipa = cyl_height_ipa(ipa_index) + ipa_radius(ipa_index) + tank_wall_thickness;
total_height_lox = cyl_height_lox(ipa_index) + lox_radius(lox_index) + tank_wall_thickness;

%% FORMATTED OUTPUT
fprintf("\nIPA Tank Volume: %.4f ft^3\n", ipa_vol);
fprintf("LOX Tank Volume: %.4f ft^3\n", lox_vol);
fprintf("IPA mass: %.02f lbm\n", ipa_mass);
fprintf("LOX mass: %.02f lbm\n", lox_mass);

fprintf("\nWall Thickness (INPUT): %.4f ft (%.4f in.)\n", tank_wall_thickness, tank_wall_thickness * 12);
fprintf("IPA Internal Radius: %.4f ft (%.4f in.)\n", ipa_radius(ipa_index), ipa_radius(ipa_index) * 12);
fprintf("LOx Internal Radius: %.4f ft (%.4f in.)\n", lox_radius(lox_index), lox_radius(lox_index) * 12);
fprintf("IPA Height: %.2f ft (%.4f in.)\n", total_height_ipa, total_height_ipa * 12);
fprintf("LOx Height: %.2f ft (%.4f in.)\n", total_height_lox, total_height_lox * 12);


