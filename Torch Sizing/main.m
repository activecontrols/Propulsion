tic
%%design criteria:
%Pc Between--- under 250 and above 90
%%mass flow is 1%
%%OF between 1.5-3
%%throat diameter is more than 0.0065
%%[data_return] = callCEA(type, Pc, PcUnits, mixture_type, mix_val, exit_cond_type, exit_cond(area ratio), fuel, fuel_temp_units, fuel_t, fuel_wt, ox, ox_temp_units, ox_t, ox_wt)
%%change units of throat area to inches
%%change pc units
%%call Cea
clear; close all; clc
%OF=linspace(1,3.5,10);
OF = [1, 1.5,2, 2.5,3, 3.5];
mdot = linspace(.0272, .06, 15);
PC = zeros(length(OF), length(mdot));
H = zeros(length(OF), length(mdot));
Tg = zeros(length(OF), length(mdot));
for j = 1:length(OF)
    a = OF(j)
    for i = 1:length(mdot)
        out = otherguessingmachine(.15,.1925,mdot(i),OF(j));
        Cp(j,i) = out(2);
        PC(j,i) = out(1);
        Tg(j,i) = out(3);
        hg(j,i) = out(4);
        
    end 
    H(j,:) = Cp(j,:) .*Tg(j,:) .* mdot .* 0.4535924;
end
h_tadpole = 6.1896e+06; %J/kg
tadpole_mdot = 2.72;
H_tadpole = h_tadpole * tadpole_mdot * 0.4535924;
H_tadpole_vector = ones(length(H)).* H_tadpole;
%PC= otherguessingmachine(0.15,0.1925,0.0285, OF);
%%EXTRACT OUTPUTS
% CEAexit = cell(6, 1);c    
% for i = 1:length(OF)
%     CEAexit{i} = callCEA('fr', PC(i), 'psi', 'o/f', OF(i), 'sup', 1, 'H2', 'K', 297, 100, 'O2', 'K', 297, 100);
%     throatcond{i}=[OF(i) CEAexit(i)];
% end
%% plot
% figure(1)
% 
% plot(OF,out(1,1,:),"Color","c","LineStyle","-",LineWidth=1.5)
% title('AC Torch Chamber Pressure vs. O/F Ratio (Dthroat=0.15in,Mdot=0.0272lb/s) -Rahaf Adi')
% xlabel('O/F Raio')
% ylabel('Chamber Pressure(psi)')
% grid on
% 


figure(2) 
plot(mdot, Tg(1,:),mdot, Tg(2,:), mdot, Tg(3,:), mdot, Tg(4,:), mdot, Tg(5,:), mdot,Tg(6,:) );
title("Combustion temperature vs Mass flow and OF ratio");
xlabel("Mass flow rate [lbm/s]");
ylabel("Temperature [K]");
legend("OF = 1", "OF = 1.5", "OF = 2", "OF = 2.5", "OF = 3", "OF = 3.5");
figure(3) 
plot(mdot, H(1,:),mdot, H(2,:), mdot, H(3,:), mdot, H(4,:), mdot, H(5,:), mdot, H(6,:), mdot, H_tadpole_vector .* .01 );
title("Total Enthalpy vs Mass flow and OF ratio");
xlabel("Mass flow rate [lbm/s]");
ylabel("Temperature [K]");
legend("OF = 1", "OF = 1.5", "OF = 2", "OF = 2.5", "OF = 3", "OF = 3.5","1% tadpole total enthalpy")
figure(4) 
plot(mdot, PC(1,:),mdot, PC(2,:), mdot, PC(3,:), mdot, PC(4,:), mdot, PC(5,:), mdot, PC(6,:) );
title("Chamber Pressure vs Mass flow and OF ratio");
xlabel("Mass flow rate [lbm/s]");
ylabel("Temperature [K]");
legend("OF = 1", "OF = 1.5", "OF = 2", "OF = 2.5", "OF = 3", "OF = 3.5");
figure(5)
plot(mdot, hg(1,:),mdot, hg(2,:), mdot, hg(3,:), mdot, hg(4,:), mdot, hg(5,:), mdot, hg(6,:));
title("Film Coeff (Bartz) vs Mass flow and OF ratio");
xlabel("Mass flow rate [lbm/s]");
ylabel("Temperature [K]");
legend("OF = 1", "OF = 1.5", "OF = 2", "OF = 2.5", "OF = 3", "OF = 3.5");



% figure(2) 
% plot(mdot, Tg(1,:),mdot, Tg(2,:), mdot, Tg(3,:) );
% title("Combustion temperature vs Mass flow and OF ratio");
% xlabel("Mass flow rate [lbm/s]");
% ylabel("Temperature [K]");
% legend("1.5","2.3","2")
% figure(3) 
% plot(mdot, H(1,:),mdot, H(2,:), mdot, H(3,:));
% title("Total Enthalpy vs Mass flow and OF ratio");
% xlabel("Mass flow rate [lbm/s]");
% ylabel("Temperature [K]");
% legend("1.5","2.3","2")
% figure(4) 
% plot(mdot, PC(1,:),mdot, PC(2,:), mdot, PC(3,:) );
% title("Chamber Pressure vs Mass flow and OF ratio");
% xlabel("Mass flow rate [lbm/s]");
% ylabel("Temperature [K]");
% legend("1.5","2.3","2")
% figure(5)
% plot(mdot, hg(1,:),mdot, hg(2,:), mdot, hg(3,:));
% title("Film Coeff (Bartz) vs Mass flow and OF ratio");
% xlabel("Mass flow rate [lbm/s]");
% ylabel("Temperature [K]");
% legend("1.5","2.3","2")

%% extract outputs
%Pc=squeeze(CEAdata('p'))
%exit_mach=squeeze(CEAdata('mach'));
%cstar = squeeze(CEAdata('cstar'));
%gamma= squeeze(CEAdata('gammas'))
%gamma=gamma(2)
%cstar = cstar(1)
%tempExit = temp(end)
%% throat equation equation
%mdot=0.033
%At = mdot*cstar / Pc
%Ae = epsilon * At;clc
%dt = 2*sqrt(At/pi);
%de = 2*sqrt(Ae/pi);
%dt_inch = meter2inch(dt)
%de_inch = meter2inch(de)
%has context menu