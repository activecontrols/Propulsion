function [thrust_req, Pc_req] = ThrustPc_Sizer(toad_mass, CF, r_t)
% Engine Size from TOAD Mass
% Adam Grendys
% Last Edited: 2/16/2025

%% INITIALIZATION
g_o = 32.174049; % (ft/sec^2)
min_thrust = 500; % (lbf) (Req. 8.4)
TWR_min = 2; % SSR Slide 18

%% CALCULATIONS
A_t = pi * r_t^2;
thrust_req = max(min_thrust, TWR_min * toad_mass);
Pc_req = thrust_req / (CF * A_t);


