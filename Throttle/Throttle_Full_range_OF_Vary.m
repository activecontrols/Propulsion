% %% PSP Active Control Throttle Code 
% % Authors: Yash Amin, Alex Suppiah
% % First Created: 1/24/2024
% % Last Updated: 

   %{ 
    Description:
    Throttle calculations for different thrust percentages and OF values

    Program Methodology
    - 
    
    Inputs:
    - 
    
    Outputs: 
    - 
    
    Assumptions:
    - Steady state
    - Frozen chemistry

   %}


%% INITIALIZATION
clear;
clc;
close all;
fclose all;
CEA_input_name = 'throttle';


%% SYSTEM CONSTANTS
breakpoints = 7;
g = 9.81;   % Gravity [m/s^2]
g_imperial = 32.2 * 12; % Gravity [in/s^2]


%% CHAMBER PARAMETERS
R_t = 0.7335/39.37; % throat radius [m]
At =  pi*R_t^2; % Throat Area [meters]
exp_ratio = 2.8153; % Expansion ratio;
Pc_max = 250;   % Max Throttle Chamber Pressure [psi]
Pe_max = 17;    % Max Throttle Exit Pressure [psi]
cstar_eff = 0.92;   % C* efficency
cf_eff = 0.95;  % Cf efficiency;
throttle_pct = linspace(0.4, 1, breakpoints);
thrust_max = 2446.52;   % Max Thrust [N]
Pa = 14.7;  % Atmospheric Pressure [P]


%% INJECTOR GEOMETRY
injector_type = "pintle";   % injector type (i.e. pintle, impinging, coax shear, ect.)
pintle_center = "OX";       % which propellant is centered (injected through the holes on the pintle post)
hole_diameter = 0.031;  % in 
hole_number = 66;           % number of holes on pintle tip
A_OX = pi * hole_diameter ^ 2 / 4 * hole_number;              % ox orifice area [in^2]
A_FUEL = 0.040310068980444;             % fuel orifice area [in^2]
cd_OX = .7;                % ox orifice discharge coefficient [N/A]
cd_FUEL = .7;              % fuel orifice discharge coefficient [N/A]
annulus_width = 0.017;     % annulus width [in]
shaft_length = 0.7500;     % Shaft length [in]
shaft_radius = 0.3750;     % Shaft radius [in]
chamber_diameter = 3.75;   % Chamber Diameter [in]

%% PROPELLANT VALUES
fuel_temp = 293.15; % [K]
oxidizer_temp = 90.17; % [K]
fuel = 'C3H8O,2propanol'; % fuel definition
fuel_weight = 100; % 
oxidizer = 'O2(L)'; % oxidizer definition
OF = linspace(0.8, 2, breakpoints); % oxidizer/fuel ratio
mdot = 1.2566; % Propellant mass flow rate [kg/s]
rho_FUEL = 49.06838 / 12^3; % fuel density [lbm/in^3]
rho_OX = 56.34 / 12^3;      % oxidizer density [lbm/in^3]


%% MATRIX INITIALIZATION
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
OF_throttle = zeros(1,breakpoints);
cstar_cea = zeros(1,breakpoints);
cf_cea = zeros(1,breakpoints);
Pe_cea = zeros(1,breakpoints);
Tc_ns_throttle = zeros(1,breakpoints);
annulus_velocity = zeros(1,breakpoints);
radial_velocity = zeros(1,breakpoints);
tmr = zeros(1,breakpoints);
lmr = zeros(1,breakpoints);
spray_angle = zeros(1,breakpoints);
impingement_point = zeros(1,breakpoints);


%% THROTTLE ITERATION
% Iterate over OF and throttle percent
 for j = 1:length(OF)   
    for i=1:length(throttle_pct)
        
        Target_thrust = throttle_pct(i)*thrust_max;
        converged = 0;
        Pc_Max = 300;   % [psi]
        Pc_Min = 50;    % [psi]
        counter = 0;
        Pc_throttle_guess = Pc_max * throttle_pct(i);   % Inital guess[psi]
        
        while ~(converged)
            % Enter CEA Call
           
            [cstar_cea, cf_cea, ~, ~, ~, ~, Pe_cea, Tc_ns, ~, ~, ~, ~, ~, ~, ~] = throttleCEA(Pc_throttle_guess, 0, fuel, fuel_weight, fuel_temp, oxidizer, oxidizer_temp, OF(j), 0, exp_ratio, 2, 1, 0, CEA_input_name);
            
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
        throttle_thrust_actual(i,j) = thrust_guess * 0.224809;   % [lbf]
        Pc_throttle_actual(i,j) = Pc_throttle_guess;   % [Psi]
        Pe_throttle(i,j) = Pe_cea;
        mdot_guess = mdot_guess * 2.20462; % [lbm/s]
        mdot_throttle_actual(i,j) = mdot_guess;  % [lbm/s]
        fuel_massflow_rate(i,j) = mdot_guess / (1 + OF(j));  % [lbm/s]
        ox_massflow_rate(i,j)= mdot_guess - fuel_massflow_rate(i,j); % [lbm/s]
        isp_throttle(i,j) = isp_actual; % [s]
        Tc_ns_throttle(i,j) = (Tc_ns - 273.15) * 9/5 + 32; % [F]
    
        % Injector pressure calculations
        P_OX_manifold(i,j) = (ox_massflow_rate(i,j) / (cd_OX * A_OX)) ^ 2 / (2 * rho_OX * g_imperial) + Pc_throttle_actual(i,j);   % [Psi]
        P_FUEL_manifold(i,j) = (fuel_massflow_rate(i,j) / (cd_FUEL * A_FUEL)) ^ 2 / (2 * rho_FUEL * g_imperial) + Pc_throttle_actual(i,j); % [Psi]
        OX_stiffness(i,j) = (P_OX_manifold(i,j) - Pc_throttle_actual(i,j)) / Pc_throttle_actual(i,j); 
        FUEL_stiffness(i,j) = (P_FUEL_manifold(i,j) - Pc_throttle_actual(i,j)) / Pc_throttle_actual(i,j);
        
        % Separation Pressure
        P_sep(i,j) = 2/3 *(Pc_throttle_actual(i,j)/Pa)^-0.2 * Pa; % [Psi]

        % Injection velocity
        annulus_velocity(i,j) = fuel_massflow_rate(i,j) / (rho_FUEL * A_FUEL);
        radial_velocity(i,j) = ox_massflow_rate(i,j) / (rho_OX * A_OX); % in/s

        % Momentum ratios
        tmr(i,j) = radial_velocity(i,j) * ox_massflow_rate(i,j) / (annulus_velocity(i,j) * fuel_massflow_rate(i,j)); % Total Momentum Ratio (radial / annulus)
        lmr(i,j) = (rho_OX * radial_velocity(i,j)^2 * A_OX / hole_number) / (rho_FUEL * annulus_velocity(i,j)^2 * hole_diameter * annulus_width);

        % Sprau calculations
        spray_angle(i,j) = (0.67 * atan(1.8 * lmr(i,j)) * 180 / pi) + 20; % blakely, freeberg, and hogge (2019) spray formation from pintle injector system
        impingement_point(i,j) = ((chamber_diameter / 2 - shaft_radius) / tan(spray_angle(i,j) / 180 * pi)) + shaft_length; % based only on simple trig

    end
 end

% Save workspace and data for post processing 
save('data_40%_OF_08_2') % percent throttle and OF range for this run


%% FIGURES

% OF vs flow rate 
f=figure('Name', 'Throttle OF and Flow Rate Trend');
set(gcf,'color','w')
colormap jet
hAxes.TickLabelInterpreter = 'latex';
contourf(fuel_massflow_rate, ox_massflow_rate, ox_massflow_rate./fuel_massflow_rate)
hbar = colorbar;
title("Throttle OF and Flow Rate Trend",'Interpreter','latex')
xlabel("Fuel Flow Rate [lb/s]",'Interpreter','latex')
ylabel("Ox Flow Rate [lb/s]",'Interpreter','latex')
ylabel(hbar, "OF Ratio",'Interpreter','latex')
exportgraphics(f,'OF_FlowRate_contour.png','Resolution',600)


% Thrust vs flow rate 
f=figure('Name', 'Throttle Thrust and Flow Rate Trend');
set(gcf,'color','w')
colormap parula
hAxes.TickLabelInterpreter = 'latex';
contourf(fuel_massflow_rate, ox_massflow_rate, throttle_thrust_actual)
hbar = colorbar;
title("Throttle Thrust and Flow Rate Trend",'Interpreter','latex')
xlabel("Fuel Flow Rate [lb/s]",'Interpreter','latex')
ylabel("Ox Flow Rate [lb/s]",'Interpreter','latex')
ylabel(hbar, "Thrust [lbf]",'Interpreter','latex')
exportgraphics(f,'Thrust_FlowRate_contour.png','Resolution',600)


% Temperature vs flow rate 
f=figure('Name', 'Chamber Temperature and Flow Rate Trend');
set(gcf,'color','w')
colormap hot
hAxes.TickLabelInterpreter = 'latex';
contourf(fuel_massflow_rate, ox_massflow_rate, Tc_ns_throttle)
hbar = colorbar;
title("Chamber Temperature and Flow Rate Trend",'Interpreter','latex')
xlabel("Fuel Flow Rate [lb/s]",'Interpreter','latex')
ylabel("Ox Flow Rate [lb/s]",'Interpreter','latex')
ylabel(hbar, "Chamber Temperature [F]",'Interpreter','latex')
exportgraphics(f,'Temp_FlowRate_contour.png','Resolution',600)


% Injector Stiffness vs flow rate 
f=figure('Name', 'Oxidizer Stiffness and Flow Rate Trend');
set(gcf,'color','w')
colormap hot
hAxes.TickLabelInterpreter = 'latex';
contourf(OX_stiffness, FUEL_stiffness,tmr)
hbar = colorbar;
title("Injector Stiffness and Momentum Ratio",'Interpreter','latex')
xlabel("Oxidizer Stiffness",'Interpreter','latex')
ylabel("Fuel Stiffness",'Interpreter','latex')
ylabel(hbar, "Total Mixture Ratio",'Interpreter','latex')
exportgraphics(f,'Stiffness_Momentum_contour.png','Resolution',600)

% Momentum Ratio vs flow rate 
f=figure('Name', 'Total Momentum Ratio and Flow Rate Trend');
set(gcf,'color','w')
colormap parula
hAxes.TickLabelInterpreter = 'latex';
contourf(fuel_massflow_rate, ox_massflow_rate, tmr)
hbar = colorbar;
title("Total Momentum Ratio and Flow Rate Trend",'Interpreter','latex','FontSize',12)
xlabel("Fuel Flow Rate [lb/s]",'Interpreter','latex')
ylabel("Ox Flow Rate [lb/s]",'Interpreter','latex')
ylabel(hbar, "Total Mixture Ratio",'Interpreter','latex','FontSize',11.5)
exportgraphics(f,'MassFlow_Momentum_contour.png','Resolution',600)

eqn_best_fit=polyfit(throttle_thrust_actual(:,3),mdot_throttle_actual(:,3),1);
syms T Mdot
Mdot_Polyfit =Mdot == eqn_best_fit(1)*T + eqn_best_fit(2);
disp((Mdot_Polyfit))
thrust_range = linspace(200,500);
Mdot_fit = eqn_best_fit(1).*(thrust_range) + eqn_best_fit(2);

% Momentum Ratio vs flow rate 
f=figure('Name', 'Massflow and Thrust (Constant OF = 1.2)');
set(gcf,'color','w')
hold on
hAxes.TickLabelInterpreter = 'latex';
plot(mdot_throttle_actual(:,3), throttle_thrust_actual(:,3))
plot(Mdot_fit,thrust_range)
legend('Actual Thrust', 'Polyfit', 'Interpreter','latex')
title("Total Mass Flow and Thrust",'Interpreter','latex','FontSize',12)
xlabel("Total Flow Rate [lb/s]",'Interpreter','latex')
ylabel("Thrust [lbf]",'Interpreter','latex')
exportgraphics(f,'MassFlow_Thrust.png','Resolution',600)




% % Throttle results for chamber @ OF 1.2
% f=figure('Name', 'Throttle Chamber Results');
% 
% subplot(2,2,1)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, throttle_thrust_actual * 0.224809)
% title("Thrust",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Thrust (lbf)",'Interpreter','latex')
% 
% subplot(2,2,2)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, Pe_throttle)
% plot(throttle_pct*100,P_sep)
% legend('$$P_{e}$$','$$P_{sep}$$','Interpreter','latex','Location','northwest')
% title("Exit Pressure",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Pressure (Psi)",'Interpreter','latex')
% 
% subplot(2,2,3)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, Pc_throttle_actual)
% title("Chamber Pressure",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Pressure (psi)",'Interpreter','latex')
% 
% subplot(2,2,4)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, isp_throttle)
% title("Isp",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Isp (s)",'Interpreter','latex')
% 
% exportgraphics(f,'Chamber_throttle.png','Resolution',600)
% 
% % Throttle results for injector @ OF 1.2
% f=figure('Name', 'Injector Throttle Results');
% 
% subplot(2,2,1)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, P_OX_manifold)
% title("Oxidizer Manifold Pressure",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Pressure (psi)",'Interpreter','latex')
% 
% subplot(2,2,2)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, P_FUEL_manifold)
% title("Fuel Manifold Pressure (psi)",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Pressure (psi)",'Interpreter','latex')
% 
% subplot(2,2,3)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, OX_stiffness*100)
% title("Oxidizer Stiffness",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Stiffness $$(\%)$$",'Interpreter','latex')
% 
% subplot(2,2,4)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, FUEL_stiffness*100)
% title("Fuel Stiffness",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Stiffness $$(\%)$$",'Interpreter','latex')
% 
% exportgraphics(f,'Injector_throttle.png','Resolution',600)
% 
% % Massflows @ OF 1.2
% f=figure('Name', 'Throttle Massflows');
% 
% subplot(2,2,[1 2])
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, mdot_throttle_actual)
% title("Total Mass Flow Rate",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Mass Flow (lbm/s)",'Interpreter','latex')
% 
% subplot(2,2,3)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, fuel_massflow_rate)
% title("Fuel Massflow Rate",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Mass Flow (lbm/s)",'Interpreter','latex')
% 
% subplot(2,2,4)
% hold on 
% set(gcf,'color','w')
% hAxes.TickLabelInterpreter = 'latex';
% plot(throttle_pct*100, ox_massflow_rate)
% title("Oxidizer Massflow Rate",'Interpreter','latex')
% xlabel("Throttle $$(\%)$$",'Interpreter','latex')
% ylabel("Mass Flow (lbm/s)",'Interpreter','latex')
% exportgraphics(f,'MassFlow_throttle.png','Resolution',600)
% 
fclose all;