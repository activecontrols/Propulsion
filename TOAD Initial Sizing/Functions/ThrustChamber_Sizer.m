function [thrust_req, Pc_req] = ThrustChamber_Sizer(max_thrust, CF, r_t)
% Engine Size from Vehicle Mass
% Adam Grendys
% Last Edited: 2/16/2025

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing

%% INITIALIZATION
min_thrust = 100; % (lbf) (Req. 8.4)

%% CALCULATIONS
A_t = pi * r_t^2;
thrust_req = max(min_thrust, max_thrust);
Pc_req = thrust_req / (CF * A_t);
