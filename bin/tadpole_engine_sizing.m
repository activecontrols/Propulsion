%% Tadpole Chamber Sizing Code
% Author: Andrew Radulovich
% Last Updated: 
%{ 
Description: This Code is used to determine the performance parameters of 
the tadpole engine given relevant design parameters for its mission. 
These calculations include throttle performance to a set percentage of 
nominal thrust.Optionally this code can call a seperate function to 
generate the chamber geometry.

Inputs: 
    - Nominal Thrust 
    - Chamber Pressure
    - Propellant Combination
    - Oxidizer Fuel Ratio
    - Nominal Exit Pressure
    - Atmospheric Pressure
    - Minimum Throttle thrust percentage
    - Cstar efficiency 
    - Cf efficiency 
    - L star 
    - Frozen or Equilibrium CEA Option
Outputs: 
    - Nominal Performance Parameters
    - Throttle Performance Parameters (Optional)
    - Chamber geometry Parameters (Optional) 
    - Chamber contour output (Optional) 

Assumptions: Uses NASA CEA for chemical equilibrium calculations thus
assumptions stated in the User Manual apply. Steady State, 
characteristics are uniform at any point along an axial slice of
the engine. 
%}

clear; clc; close all
%% inputs
% File directory
filename = "Tadpole_Sizing_Parameters.txt";
% Calculation Options
throttlebool = 0; % 1: Perform Throttle Calculations
geobool = 1; 
numpoints = 5;
tolerance = .01; % Can run into convergence issues if tolerance is lowered
maxiter = 50; 
% Engine Desin Input Parameters
Tf = 550; % Thrust [lbf]
minthrottle_percent = 40; % Minimum thrust percentage of nominal %
Pc = 250; % Nominal Chamber Pressure [psi]
Pe = 16.5; % Nominal Exit Pressure [psi]
%fuel = 'C3H8O,2propanol';
fuel = 'CH3OH'; % fuel
ox = 'O2(L)'; % Oxidizer
Lstar = 45; % Characteristic Length [in]
OF = .8; % OF ratio
cstareff = .92; % Cstar efficiency 
cfeff = .95; % Cf efficiency
Pa = 14.7; % Atmospheric Pressure
R_fillets = [1.5,1.5,.5];

%% initalization and precalculations
minthrottle_percent = minthrottle_percent * .01;
converge_vector = ones(1,numpoints);
PcPe = Pc/Pe;
g = 9.80665;
cstar_vector = ones(1, numpoints); 
mdot_vector = ones(1,numpoints); 
cf_vector = ones(1,numpoints);
Pc_vector = ones(1,numpoints);
Pe_vector = ones(1,numpoints);
gamma_vector = ones(1,numpoints);
thrust_vector = linspace(Tf, Tf*minthrottle_percent,numpoints);


%% Nominal (Max) Thrust Sizing
[out2] = callCEA('fr', Pc, 'psi', 'o/f', OF,'pip',PcPe,fuel,'K',440,100,ox,'K',90.18,100);
cstar_out = squeeze(out2('cstar'));
cstar = cstar_out(1) * cstareff;
cf_out = squeeze(out2('cf'));
cf_cea = cf_out(2); 
CEA_isp = squeeze(out2('isp'));
ue = CEA_isp(2); 
epsilon_out = squeeze(out2('ae/at'));
gamma_out = squeeze(out2('gammas'));
gamma = gamma_out(1);
%gamma_t = gamma_out(2);
%gamma_t = (gamma + gamma_out(2))/2;
epsilon = epsilon_out(2);
%epsilon_kamon = 1/(((gamma_t+1)/2)^(1/(gamma_t-1)) * (Pe/Pc)^(1/gamma_t) * ((gamma_t+1)/(gamma_t-1)*(1-(Pe/Pc)^((gamma_t-1)/gamma_t)))^.5);
%cf = (sqrt(((2*gamma^2)/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1-(Pe/Pc)^((gamma-1)/gamma))) + (Pe-Pa)*epsilon/Pc) * cfeff;
cf = (cf_cea + (Pe-Pa)*epsilon/Pc) * cfeff;
Isp = cstar * cf / g;
mdot = Tf/Isp;
At = (cstar*mdot)/ (g * Pc);
Dt = 2*(At/pi())^(1/2);
rt = Dt/2; 
Ae = At * epsilon;
De = 2*(Ae/pi())^(1/2);
% Engine Geometry 
Vc = Lstar * At; 



%% Throttle Thrust Sizing
if throttlebool
    Pc_guess = Pc;
    j = 1;
    for thrust = thrust_vector
        Pc_max = 1000;
        Pc_min = 0;
        convergence = 0;
        i = 1;
        while convergence == 0
            [out1] = callCEA('fr',Pc_guess,'psi','o/f',OF,'sup', epsilon, fuel,'K', 293.15, 100, ox, 'K', 90.18, 100);
            cf_out1 = squeeze(out1('cf')); 
            cf_cea1 = cf_out1(2);
            p_out = squeeze(out1('p')); 
            Pe_vector(j) = p_out(2) * 0.0001450377;
            %cf_vector(j) = (sqrt(((2*gamma_vector(j)^2)/(gamma_vector(j)-1)) * (2/(gamma_vector(j)+1))^((gamma_vector(j)+1)/(gamma_vector(j)-1)) * (1-(Pe/Pc)^((gamma_vector(j)-1)/gamma_vector(j)))) + (Pe-Pa)*epsilon/Pc_guess) * cfeff;
            cf_vector(j) = (cf_cea1 + (Pe_vector(j)-Pa)*epsilon/Pc_guess) * cfeff;
            thrust_guess = Pc_guess * cf_vector(j) * At;
            a = thrust_guess - thrust
            if abs(thrust_guess - thrust) < tolerance || i == maxiter
                convergence = 1;
                if i == maxiter
                    converge_vector(1,i) = 0;
                end
            elseif (thrust_guess - thrust) > 0 
                Pc_max = Pc_guess;
                Pc_guess = (Pc_guess + Pc_min)/2
            else 
                Pc_min = Pc_guess;
                Pc_guess = (Pc_guess + Pc_max)/2
            end
            i = i + 1;
        end
        p_out = squeeze(out1('p')); 
        Pe_vector(j) = p_out(2) * 0.0001450377;
        Pc_vector(j) = Pc_guess;
        cstar_out = squeeze(out1('cstar'));
        cstar_vector(j) = cstar_out(1) * cstareff;
        mdot_vector(j) = g*Pc_vector(j)*At/(cstar_vector(j));
        j = j + 1
    end
    
    Isp_vector = thrust_vector ./ mdot_vector;
    Psep = .667.*((Pc_vector./Pa).^(-.2)).*Pa;
    
    
    figure(1)
    plot(thrust_vector, mdot_vector)
    title("mass flow rate v.s. thrust")
    figure(2)
    plot(thrust_vector, Isp_vector)
    title("Isp v.s thrust")
    figure(3)
    plot(thrust_vector, Pc_vector)
    title("chamber pressure vs thrust"); 
    figure(4)
    plot(thrust_vector, Pe_vector, thrust_vector, Psep)
    title("Exit pressure vs thrust")
    figure(5) 
    subplot(2,2,1)
    plot(thrust_vector,Pc_vector, linewidth=2)
    title("Chamber Pressure vs. Thrust")
    ylabel("Chamber Pressure [psi]")
    xlabel("Thrust [lbf]")
    xlim([220,550]);
    subplot(2,2,2)
    plot(thrust_vector,mdot_vector, linewidth=2)
    title("Total Mass Flow Rate vs. Thrust")
    ylabel("Total Mass Flow [lbm/s]")
    xlabel("Thrust [lbf]")
    xlim([220,550]);
    subplot(2,2,3)
    plot(thrust_vector,Pe_vector, linewidth=2)
    hold on
    plot(thrust_vector,Psep, 'r--' , linewidth=2)
    title("Exit Pressure vs. Thrust")
    ylabel("Exit Pressure [psi]")
    xlabel("Thrust [lbf]")
    legend("Designed Exit Pressure", "Kalt-Bendall Seperation Criterion")
    xlim([220,550]);
    
    subplot(2,2,4)
    plot(thrust_vector,Isp_vector, linewidth=2)
    title("Isp vs. Thrust")
    ylabel("Isp [sec]")
    xlabel("Thrust [lbf]")
    xlim([220,550]);
end

if geobool 
    R_t = Dt/2;
    Ac = (pi*((3.75/2)^2));
    con_ratio = Ac/At;
    conv_angle = 37;
    theta_i = 15;
    [x_contour, r_contour, L_c, L_total, L_converging, L_diverging, L_seg] = engineContour("conical", .8, R_t, epsilon, con_ratio, conv_angle, theta_i, 0,Lstar, R_fillets, 1,100);
end






