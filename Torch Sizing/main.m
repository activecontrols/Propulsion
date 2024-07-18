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
OF=linspace(1,3, 20);
PC= otherguessingmachine(0.15,0.1925,0.0272, OF);
 
%%EXTRACT OUTPUTS
CEAexit = cell(6, 1);
for i = 1:length(OF)
    CEAexit{i} = callCEA('fr', PC(i), 'psi', 'o/f', OF(i), 'sup', 1, 'H2', 'K', 297, 100, 'O2', 'K', 297, 100);
    throatcond{i}=[OF(i) CEAexit(i)];
end
%% plot
figure

plot(OF,PC,"Color","c","LineStyle","-",LineWidth=1.5)
title('AC Torch Chamber Pressure vs. O/F Ratio (Dthroat=0.1875in,Mdot=0.0272lb/s) -Rahaf Adi')
xlabel('O/F Raio')
ylabel('Chamber Pressure(psi)')
grid on
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