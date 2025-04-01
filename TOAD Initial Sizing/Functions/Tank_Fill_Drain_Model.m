% Simulation of Tank Draining and Fill

% Author: Rugved Dikay
% Date: 12/24/2024

% Description: As part of the process of pressurization system design for
% the TOAD lander, this script is a simple representation of a tank drain
% and fill operation from a COPV to a theoretical fuel tank. This code is
% made to mimic a PSP Liquids Script, which is the basis for validating
% this script.

clear;
close all;

%% Param Definition

% Assumptions --> Let's get these reviewed!

% 1. Isenthalpic expansion assumed for pressurant into propellant tank
% (diffuser design condition)
% 2. Isentropic expansion assumed for fluids expanding into a tank volume
% 3. Inlet press gas temp modeled assuming isentropic expansion of helium
%    into prop. tank
% 4. Press gas in tank never mixes with propellant
% 5. Thermally insulated tanks
% 6. When in contact in prop. tank, press gas and prop fluids never mix
% 7. Only tank 1 fluid expands into tank 1 open volume. Tank 2 fluid only
% drains
% 8. Initial ullage of 10% assumed (standard, derived from industry and
% PSP-L)
% 9. Isentropic expansion and temperature drop of helium in COPV
% 10. Thermal mass averaging used for determining ullage temp as f(time)
% 11. Initial pressure of prop tank is desired pressure

% NOTE: All model calculations are done in metric units for simplicity and
% converted back to English/convenient units as part of post processing

%% User Inputs

% Fluids
press_name = 'helium'; %Name of press gas
prop_name = 'oxygen'; %Name of propellant

% Press System Properties
P1 = 4950; %Initial COPV pressure, psia
P_setpoint = 600; %Desired prop tank pressure, psia
T1 = 294.261; %Initial COPV gas temp, K
T2 = 90.14; %Initial propellant temp, K
V1 = 550; %Internal volume of COPV, in^3

% Propellant Properties
m2_tot = 100; %Total propellant mass, lb
mdot_prop = 1; %Mass flow rate exiting prop tank, lb/s

%% Model Setup

% Unit Conversion

P1 = P1*6894.76; %P1 in Pa
P_setpoint = P_setpoint*6894.76; %P_setpoint in Pa
V1 = V1*0.0254^3; %COPV Volume in m^3
mdot_prop = mdot_prop/2.205; %mdot_prop in kg/s
m2_tot = m2_tot/2.205; %Total prop mass in kg

% Fixed Quantities

gamma = gammaCalc(T1,P1,press_name); %Specific heat of press gas
R_spec = R_calc(T1,P1,press_name); %Specific gas law constant for Helium, J/kg-K

% Initializations

COPV_rho = densityCalc(T1,P1,press_name); %Density of gas in COPV, initially, kg/m^3
m1_tot = COPV_rho*V1; %Total press mass available, kg
ullage_rho = densityCalc(T2,P_setpoint,press_name); %Density of ullage initially, kg/m^3
prop_rho = densityCalc(T2, P_setpoint, prop_name); %Density of propellant, kg/m^3
V2_prop = m2_tot / prop_rho; %Total volume of propellant, m^3
V2 = V2_prop/0.9; %Total volume of propellant tank, m^3
V2_ullage = V2*0.1; %Initial ullage volume, m^3 (10% assumption)
m_ullage = ullage_rho*V2_ullage; %Mass of press gas in ullage, kg

% Data Vectors
P1_vec = []; P1_vec(end+1) = P1;
T1_vec = []; T1_vec(end+1) = T1;
m1_vec = []; m1_vec(end+1) = m1_tot;
T2_vec = []; T2_vec(end+1) = T2;
m2_vec = []; m2_vec(end+1) = m2_tot;
V2_ullage_vec = []; V2_ullage_vec(end+1) = V2_ullage;
m_ullage_vec = []; m_ullage_vec(end+1) = m_ullage;
ullage_rho_vec = []; ullage_rho_vec(end+1) = ullage_rho;
mdot_press = []; mdot_press(end+1) = 0;
press_check = []; press_check(end+1) = P_setpoint/6894.76; %Checks pressure of prop. tank for each time step and stores value
press_inlet_temp = []; press_inlet_temp(end+1) = T1;

%% Computations

dt = 0.01; %Timestep, sec.
timeline = 0:dt:40; %Timeline of tank sim. in seconds

for i = 2:length(timeline) %Runs system simulation for specified duration
    % Current System Properties
    P1 = P1_vec(end);
    T1 = T1_vec(end);
    m1 = m1_vec(end);
    V2_ullage = V2_ullage_vec(end);
    m_ullage = m_ullage_vec(end);
    T2 = T2_vec(end);
    m2 = m2_vec(end);

    % Volume & Pressure Change in Prop Tank
    dV_tank = mdot_prop*dt/prop_rho;
    V2_ullage_vec(end+1) = V2_ullage + dV_tank;

    % Calculating press mass needed to maintain pressure
    ullage_rho = m_ullage / V2_ullage_vec(end);
    press_gas_T = isentropic_temp(P1,P_setpoint,T1,gamma);
    mdot_press(end+1) = (densityCalc(T2,P_setpoint,press_name)*V2_ullage_vec(end) - m_ullage)/dt;

    % Ullage Temperature Change (Mass Averaged)
    T2_vec(end+1) = ((mdot_press(end)*dt*press_gas_T) + m_ullage*T2) / (m_ullage + mdot_press(end)*dt);

    % COPV Property Update
    m1_vec(end+1) = m1 - mdot_press(end)*dt;
    P1_vec(end+1) = refpropm('P','T',T1,'D',m1_vec(end)/V1,press_name)*1000;
    if P1_vec(end) <= 845*6894.76 %Ensures press gas flow is choked
        T1_vec(end+1) = isentropic_temp(P1,P1_vec(end),T1,gamma);

        % Prop Tank Property + General Updates
        m2_vec(end+1) = m2 - mdot_prop*dt;
        m_ullage_vec(end+1) = m_ullage + mdot_press(end)*dt;
        ullage_rho_vec(end+1) = ullage_rho;
        press_check(end+1) = (refpropm('P','T',T2,'D',m_ullage_vec(end)/V2_ullage(end),press_name)*1000)/6894.76;
        press_inlet_temp(end+1) = press_gas_T;
        break
    end
    T1_vec(end+1) = isentropic_temp(P1,P1_vec(end),T1,gamma);

    % Prop Tank Property + General Updates
    m2_vec(end+1) = m2 - mdot_prop*dt;
    m_ullage_vec(end+1) = m_ullage + mdot_press(end)*dt;
    ullage_rho_vec(end+1) = ullage_rho;
    press_check(end+1) = (refpropm('P','T',T2,'D',m_ullage_vec(end)/V2_ullage(end),press_name)*1000)/6894.76;
    press_inlet_temp(end+1) = press_gas_T;

end

%% Model Validation

end_idx = length(mdot_press);
mdot_est = trapz(timeline(1:end_idx), mdot_press);

% Error in mass calculations
mdot_error1 = (m1_tot - m1_vec(end)) - (m_ullage_vec(end) - m_ullage_vec(1));
mdot_error2 = mdot_est - (m1_tot - m1_vec(end));
mdot_error3 = mdot_est - (m_ullage_vec(end) - m_ullage_vec(1));

%% Outputs

clc;
fprintf("\n-------------------\n");
fprintf("Total COPV Press Gas Mass: %f lbf\n", m1_tot*2.205);
fprintf("Total Press Gas Removed: %f lbf\n", mdot_est * 2.205);
fprintf("Final COPV Pressure: %f psi\n", P1_vec(end)/6894.76);
fprintf("Final COPV Temperature: %f K\n\n", T1_vec(end));
fprintf("Final Ullage Mass: %f lbf\n", m_ullage_vec(end)*2.205);
fprintf("Mass Estimate Error: %f lbf\n", (mdot_error1^2 + mdot_error2^2 + mdot_error3^2)^0.5 * 2.205);
fprintf("-------------------\n");

figure(1)
plot(timeline(1:end_idx),P1_vec./6894.76);
xlabel("Time (sec.)");
ylabel("COPV Pressure (psia)");
title("COPV Pressure vs. Time");
grid on;

figure(2)
plot(timeline(1:end_idx),T1_vec);
xlabel("Time (sec.)");
ylabel("COPV Temperature (K)");
title("COPV Temperature vs. Time");
grid on;

figure(3)
plot(timeline(1:end_idx),V2_ullage_vec);
xlabel("Time (sec.)");
ylabel("Ullage Volume (m^3)");
title("Ullage Volume vs. Time");
grid on;

figure(4)
plot(timeline(1:end_idx),m_ullage_vec*2.205);
xlabel("Time (sec.)");
ylabel("Ullage Gas Mass (lb)");
title("Ullage Gas Mass vs. Time");
grid on;

figure(5)
plot(timeline(1:end_idx),mdot_press*2.205);
xlabel("Time (sec.)");
ylabel("Mass Flow Rate Required (lb/s)");
title("Required Press Mdot vs. Time");
grid on;

figure(6)
plot(timeline(1:end_idx),T2_vec);
xlabel("Time (sec.)");
ylabel("Ullage Volume Temperature (K)");
title("Ullage Volume Temperature vs. Time");
grid on;

figure(7)
plot(timeline(1:end_idx),press_check);
xlabel("Time (sec.)");
ylabel("Prop Tank Pressure (psia)");
title("Prop Tank Pressure vs. Time");
grid on;

%% Helper Functions

function rho = densityCalc(T,P,fluid_name)
% Computes density of fluid

P = P * (1/1000); % Converts pressure to kPa;
rho = refpropm('D', 'T', T, 'P', P, fluid_name);
end

function gamma = gammaCalc(T,P,fluid_name)
% Computes specific heat capacity ratio for a fluid

P = P * (1/1000); % Converts pressure to kPa
gamma = refpropm('K', 'T', T, 'P', P, fluid_name);
end

function R_specific = R_calc(T,P,fluid_name)
% Compute specific gas constant for gas species

P = P * (1/1000); %Converts pressure to kPa
molar_mass = refpropm('M','T',T, 'P', P,fluid_name);
R_specific = 8314/molar_mass;
end

function T2 = isentropic_temp(P1,P2,T1,gamma)
% Computes temperature of gas assuming isentropic expansion
T2 = (P2/P1)^(1-1/gamma) * T1;
end