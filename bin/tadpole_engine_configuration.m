%% Tadpole Configuration Calculations
% Author: Andrew Radulovich
clear; clc; close all

% Inputs
%fuel = 'CH3OH'; % fuel
fuel = 'C3H8O,2propanol';
ox = 'O2(L)'; % Oxidizer
mdot = 2.72; % Total mass flow rate lb/s
OF = 1.2; % OF ratio
epsilon = 2.904; % Supersonic area ratio
cstareff = .92;
cfeff = .95;
At = 1.6178; % Area of throat [in^2]
numpoints = 100;
Pa = 14.7;
throttlebool = 1; % Do throttle analysis?

% Constants
g = 9.80665;
gi = 32.17405;

% Initialization
Pc_guess = 250
Pc_max = 1000;
Pc_min = 0;
convergence = 0;
minthrottle_percent = 40 * .01;
tolerance = .00005 % Can run into convergence issues if tolerance is lowered
maxiter = 50; 
tolerance2 = .00001

% Nominal Thrust Calculations

while convergence == 0
    [out1] = callCEA('fr',Pc_guess,'psi','o/f',OF,'sup', epsilon, fuel,'K', 293.15, 100, ox, 'K', 90.18, 100);
    cstar = squeeze(out1('cstar'));
    mdot_guess = g*Pc_guess*At/(cstar(1)*cstareff)
    if abs(mdot_guess - mdot) < tolerance
        convergence = 1;
    elseif (mdot_guess - mdot) > 0 
        Pc_max = Pc_guess;
        Pc_guess = (Pc_guess + Pc_min)/2
    else 
        Pc_min = Pc_guess;
        Pc_guess = (Pc_guess + Pc_max)/2
    end
end
Pc = Pc_guess;
P = squeeze(out1('p'));
Pe = P(2) * 0.0001450377;
%[out1] = callCEA('fr',Pc1,'psi','o/f',OF,'sup', aeat1, fuel,'K', 298.15, 100, ox, 'K', 90.18, 100);
CEA_isp = squeeze(out1('isp')) ;
%cstar1 = squeeze(out1('cstar')) * 3.28084;
%cstar2 = 9.80665*Pc1 * At1/ (mdot1 * cstareff)* 3.28084
cstar = cstar(1) * cstareff;
cf_out = squeeze(out1('cf'));
cf_cea = cf_out(2); 
CEA_isp = squeeze(out1('isp'));
ue = CEA_isp(2); 
epsilon_out = squeeze(out1('ae/at'));
gamma_out = squeeze(out1('gammas'));
gamma = gamma_out(1);
epsilon = epsilon_out(2);
cf = (cf_cea + (Pe-Pa)*epsilon/Pc) * cfeff;
Isp = cstar * cf / g;
Tf = Isp * mdot_guess;
Tf1 = Pc * cf * At;
mdot2 = Tf/Isp;


thrust_vector = linspace(Tf1, Tf*minthrottle_percent,numpoints);


% Throttle Point Calculations 
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
            if abs(thrust_guess - thrust) < tolerance2 || i == maxiter
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