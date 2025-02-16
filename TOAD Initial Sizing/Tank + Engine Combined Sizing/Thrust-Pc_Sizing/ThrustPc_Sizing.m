% Engine Size from TOAD Mass
% Adam Grendys
% 2/16/2025

%% ASSUMPTIONS
% Using IPA-LOx
% Rev1 Tadpole Capabilites

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing

clear
clc
close all

%% INITIALIZATION
g_o = 32.174049; % (ft/sec^2)
toad_mass = linspace(100, 800, 1000); % TOAD mass (lbm)
min_thrust = 500; % (lbf) (Req. 8.4)
TWR_min = 2; % SSR Slide 18

% Engine Parameters (From Tadpole Rev1)
r_t =  1.4352 / 2; % Radius of throat (in) (From Tadpole CMM)
CF = 1.32; % Thrust Coefficient (95% ηcf)

%% CALCULATIONS
A_t = pi * r_t^2;
thrust_req = max(min_thrust, TWR_min * toad_mass);
Pc_req = thrust_req / (CF * A_t);

%% FORMATTED OUTPUT
subplot(2,1,1)
plot(toad_mass, thrust_req)
xlabel("TOAD Mass (lbm)")
ylabel("Required Thrust (lbf)")
ylim("padded")
grid on

subplot(2,1,2)
plot(toad_mass, Pc_req)
xlabel("TOAD Mass (lbm)")
ylabel("Required Pc (psia)")
ylim("padded")
grid on

sgtitle("Required Thrust and Pc vs. TOAD Mass")
