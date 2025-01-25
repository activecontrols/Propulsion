%% Throttle Code
% PSP Active Control

clear;
clc;
fclose all;
CEA_input_name = 'throttle';


%% System Constants
breakpoints = 200;

R_t = 0.7335/39.37; % throat radius [m]
At =  pi*R_t^2; % Throat Area [meters]
exp_ratio = 2.8153; % Expansion ratio;
Pc_max = 250;   % Max Throttle Chamber Pressure [psi]
Pe_max = 17;    % Max Throttle Exit Pressure [psi]
cstar_eff = 0.92;   % C* efficency
cf_eff = 0.95;  % Cf efficiency;
throttle_pct = linspace(0.4, 1, breakpoints);
thrust_max = 2446.52;   % Max Thrust [N]
g = 9.81;   % Gravity [m/s^2]
g_imperial = 32.2 * 12;            % in/s^2
Pa = 14.7;  % Atmospheric Pressure [P]


%% Injector Geometry
injector_type = "pintle";   % injector type (i.e. pintle, impinging, coax shear, ect.)
pintle_center = "OX";       % which propellant is centered (injected through the holes on the pintle post)
hole_diameter = 0.031;  % in 
hole_number = 66;           % number of holes on pintle tip
A_OX = pi * hole_diameter ^ 2 / 4 * hole_number;              % ox orifice area [in^2]
A_FUEL = 0.040310068980444;             % fuel orifice area [in^2]
cd_OX = .7;                % ox orifice discharge coefficient [N/A]
cd_FUEL = .7;              % fuel orifice discharge coefficient [N/A]
annulus_width = 0.017;     % annulus width [in]


%% Propellant Values
fuel_temp = 293.15; % [K]
oxidizer_temp = 90.17; % [K]
fuel = 'C3H8O,2propanol'; % fuel definition
fuel_weight = 100; % 
oxidizer = 'O2(L)'; % oxidizer definition
OF = 1.2; % oxidizer/fuel ratio
mdot = 1.2566; % Propellant mass flow rate [kg/s]
rho_FUEL = 49.06838 / 12^3; % fuel density [lbm/in^3]
rho_OX = 56.34 / 12^3;      % oxidizer density [lbm/in^3]


%% Matrix Initialization

throttle_thrust_actual = zeros(1, breakpoints);
Pc_throttle_actual = zeros(1, breakpoints);
mdot_throttle_actual = zeros(1, breakpoints);
fuel_massflow_rate = zeros(1, breakpoints);
ox_massflow_rate= zeros(1, breakpoints);
P_OX_manifold = zeros(1, breakpoints);
P_FUEL_manifold = zeros(1, breakpoints);
isp_throttle = zeros(1, breakpoints);
Pe_throttle = zeros(1, breakpoints);
OX_stiffness = zeros(1, breakpoints);
FUEL_stiffness = zeros(1, breakpoints);
P_sep = zeros(1, breakpoints);


%% Throttle Iteration

for i=1:length(throttle_pct)
    
    Target_thrust = throttle_pct(i)*thrust_max;
    converged = 0;
    Pc_Max = 300;   % [psi]
    Pc_Min = 50;    % [psi]
    counter = 0;
    Pc_throttle_guess = Pc_max * throttle_pct(i);   % Inital guess[psi]
    
    while ~(converged)
        % Enter CEA Call
        [cstar_cea, cf_cea, ~, ~, ~, gamma, Pe_cea, T_chamber, ~, ~, ~, ~, ~, ~, ~] = throttleCEA(Pc_throttle_guess, 0, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, exp_ratio, 2, 1, 0, CEA_input_name);
        Pc_throttle_guess_SI = Pc_throttle_guess * 6895; % [Pa]
        Pe_cea = Pe_cea / 6895; % Convert the CEA Pa to PSI
        
        % Adjust perfectly expanded parameters 
        cf_cea = cf_cea + (Pe_cea - Pa) / Pc_throttle_guess * exp_ratio; % Adjust the Cf from perfectly expanded to desired nozzle geometry
        cstar_guess_actual = cstar_cea * cstar_eff; % Adjusting C* with efficiency
        cf_guess_actual = cf_cea * cf_eff; % Adjusting Cf with efficiency;
        isp_actual = cstar_guess_actual * cf_guess_actual / g; % Adjust the isp from perfectly expanded to desired nozzle geometry
        mdot_guess = (Pc_throttle_guess_SI) * At / cstar_guess_actual; % Propellant mass flow guess
        thrust_guess = cf_guess_actual*(Pc_throttle_guess_SI)*At;
    
        % Check for convergence
        if abs(thrust_guess - Target_thrust) > 1 && counter < 250 % check for tolerance
            % convergence loop
            if thrust_guess - Target_thrust > 0
                Pc_Max = Pc_throttle_guess;
            else 
                Pc_Min = Pc_throttle_guess;
            end 
            Pc_throttle_guess = (Pc_Max+Pc_Min) / 2;
    
            counter = counter + 1;
        else
            converged = 1;
    
        end
    end
    
    % Throttle results
    throttle_thrust_actual(i) = thrust_guess;   % N
    Pc_throttle_actual(i) = Pc_throttle_guess;   % Psi
    Pe_throttle(i) = Pe_cea;
    mdot_guess = mdot_guess * 2.20462; % lbm/s
    mdot_throttle_actual(i) = mdot_guess;  % lbm/s
    fuel_massflow_rate(i) = mdot_guess / (1 + OF);  % lbm/s
    ox_massflow_rate(i)= mdot_guess - fuel_massflow_rate(i); % lbm/s
    isp_throttle(i) = isp_actual; % s
    cstar_throttle(i) = cstar_cea * cstar_eff;
    cf_throttle(i) = cf_cea * cf_eff;
    gamma_throttle(i) = gamma;
    t_throttle(i) = T_chamber;

    % Injector pressure calculations
    P_OX_manifold(i) = (ox_massflow_rate(i) / (cd_OX * A_OX)) ^ 2 / (2 * rho_OX * g_imperial) + Pc_throttle_actual(i);   % Psi
    P_FUEL_manifold(i) = (fuel_massflow_rate(i) / (cd_FUEL * A_FUEL)) ^ 2 / (2 * rho_FUEL * g_imperial) + Pc_throttle_actual(i); % Psi
    OX_stiffness(i) = (P_OX_manifold(i) - Pc_throttle_actual(i)) / Pc_throttle_actual(i); 
    FUEL_stiffness(i) = (P_FUEL_manifold(i) - Pc_throttle_actual(i)) / Pc_throttle_actual(i);
     % Seperation Pressure
    P_sep(i) = 2/3 *(Pc_throttle_actual(i)/Pa)^-0.2 * Pa;

end

[cstar_min, ~, ~, ~, ~, ~, ~, Tc_ns_min, ~, ~, ~, ~, ~, ~, ~] = throttleCEA(Pc_throttle_actual(1), Pe_throttle(1), fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF, 0, 0, 1, 1, 0, CEA_input_name);

%% FIGURES

% Throttle results for chamber
f=figure('Name', 'Throttle Chamber Results');

subplot(2,2,1)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, throttle_thrust_actual * 0.224809, "Linewidth", 2)
title("Thrust","Interpreter","latex")
xlabel("Throttle (%)","Interpreter","latex")
ylabel("Thrust (lbf)","Interpreter","latex")

subplot(2,2,2)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, Pe_throttle, "Linewidth", 2)
plot(throttle_pct*100,P_sep, "--", "Linewidth", 2)
legend('$$P_{e}$$','$$P_{sep}$$','Interpreter','latex','Location','northwest')
title("Exit Pressure","Interpreter","latex")
xlabel("Throttle (%)","Interpreter","latex")
ylabel("Pressure (Psi)","Interpreter","latex")

subplot(2,2,3)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, Pc_throttle_actual, "Linewidth", 2)
title("Chamber Pressure","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")
ylabel("Pressure [psi]","Interpreter","latex")

subplot(2,2,4)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, isp_throttle, "Linewidth", 2)
title("Isp","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")
ylabel("Isp [s]","Interpreter","latex")

exportgraphics(f,'Chamber_throttle.png','Resolution',600)

% Throttle results for injector
f=figure('Name', 'Injector Throttle Results');

subplot(2,2,1)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, P_OX_manifold, "Linewidth", 2)
title("Oxidizer Manifold Pressure",'Interpreter','latex')
xlabel("Throttle [%]","Interpreter","latex")
ylabel("Pressure [psi]","Interpreter","latex")

subplot(2,2,2)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, P_FUEL_manifold, "Linewidth", 2)
title("Fuel Manifold Pressure [psi]","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")
ylabel("Pressure [psi]","Interpreter","latex")

subplot(2,2,3)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, OX_stiffness*100, "Linewidth", 2)
title("Oxidizer Stiffness","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")
ylabel("Stiffness [%]","Interpreter","latex")

subplot(2,2,4)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, FUEL_stiffness*100, "Linewidth", 2)
title("Fuel Stiffness","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")
ylabel("Stiffness [%]","Interpreter","latex")

exportgraphics(f,'Injector_throttle.png','Resolution',600)

% Mass flows
f=figure('Name', 'Throttle Massflows');

subplot(2,2,[1 2])
hold on 
grid on
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, mdot_throttle_actual, "Linewidth", 2)
title("Total Mass Flow Rate","Interpreter","latex")
ylabel("Mass Flow [lbm/s]","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")

subplot(2,2,3)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, fuel_massflow_rate, "Linewidth", 2)
title("Fuel Massflow Rate","Interpreter","latex")
ylabel("Mass Flow [lbm/s]","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")

subplot(2,2,4)
hold on 
grid on 
set(gcf,'color','w')
hAxes.TickLabelInterpreter = 'latex';
plot(throttle_pct*100, ox_massflow_rate, "Linewidth", 2)
title("Oxidizer Massflow Rate","Interpreter","latex")
ylabel("Mass Flow [lbm/s]","Interpreter","latex")
xlabel("Throttle [%]","Interpreter","latex")
exportgraphics(f,'MassFlow_throttle.png','Resolution',600)

fclose all;