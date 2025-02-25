%% TOAD Sizer Main
% Description: An iterative method to determine foundational TOAD
% parameters: Required Thrust, Required Chamber Pressure, TOAD
% mass. Be sure to add the "functions" folder to path.
% Author: Adam Grendys
% Last Edited: 2/16/2025

%% ASSUMPTIONS/DECISIONS
% IPA-LOx
% 6061-T6 Tank Material
% Rev1 Tadpole Capabilites and Geometry
% 10% Fuel/Ox Reserves by mass
% 10% Tank Ullage volume 
% Spherical Endcaps

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing

clear
clc
close all

%% USER INPUT
selection = NaN;
while selection ~= 1 && selection ~= 2
    selection = input("Default Parameters (1) or Custom (2): ");
    if selection == 1
        mdot = 2.2; % Total mass flow rate (lbm/s)
        OF = 1.2; % Desired OF Ratio
        flight_time = 60; % Flight Duration (s)
    elseif selection == 2
        mdot = input("Input Total Mdot (lbm/s): "); % Total mass flow rate (lbm/s)
        OF = input("\nInput OF Ratio: "); % Desired OF Ratio
        flight_time = input("\nFlight time (s): "); % Flight Duration (s)
    else
        fprintf("\nINVALID SELECTION\n", 2)
    end
end

%% INITIALIATION
% Engine Parameters (From Tadpole Rev1)
r_t =  1.4352 / 2; % Radius of throat (in) (From Tadpole CMM)
CF = 1.32; % Thrust Coefficient (95% ηcf)

prop_massFraction = 0.5613; % ESTIMATE, derived from Masten's XODIAC Lander

%% CALCULATIONS
[TankProp_mass, Prop_mass] = tank_sizer(mdot, OF, flight_time);
TOAD_mass = Prop_mass / prop_massFraction;

[thrust_req, Pc_req] = ThrustPc_Sizer(TOAD_mass, CF, r_t);

%% FORMMATED OUTPUT
fprintf("\nTOAD Mass: %.3f lbm", TOAD_mass)
fprintf("\nRequired Thrust: %.3f lbf", thrust_req)
fprintf("\nRequired Chamber Pressure: %.3f psia", Pc_req)
