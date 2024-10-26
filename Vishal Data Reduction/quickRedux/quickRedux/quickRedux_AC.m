clear
clc
close all

%% FILE INTIALIZATION
testid = 1;
saveData = true;
tdmsLoc = '\Propulsion\Vishal Data Reduction'; %where you are saving directory
addpath(cd)

dataFileNameMB = 'DataLog_2024-1004-1004-07_PSPAC_Data_Wiring.tdms'; %file name

TOR = 'AC Test 1'; %Label of test

LFmatFilename = sprintf('Test_%d_Data',testid);

        LFMB = reduceTDMS(dataFileNameMB,1,nan);
 
   % Package and Save
    if saveData % && ZeroSave
        fprintf('Saving Low Frequency Data...\n')
        save([pwd,'\',LFmatFilename],'LFMB')
        fprintf('Low Frequency Data Saved.\n\n')
    end

allData.LFMB = LFMB;


%% RAW DATA PLOTTING, hardcoded :( 
figure()
title([ TOR ' Ox Supply Pressure'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_pdox_083.Value,'DisplayName',LFMB.pt_pdox_083.Name)

figure()
title([ TOR ' N2 Supply Pressure'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_n2_026.Value,'DisplayName',LFMB.pt_n2_026.Name)

figure()
title([ TOR ' H2 Supply Pressure'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_ch_058.Value,'DisplayName',LFMB.pt_ch_058.Name)

figure()
title([ TOR ' Purge Pressure'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_n2_076.Value,'DisplayName',LFMB.pt_n2_076.Name)

figure()
title([ TOR ' O2 Orifice Upstream'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igox_04.Value,'DisplayName',LFMB.pt_igox_04.Name)

figure()
title([ TOR ' H2 Orifice Upstream'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igfu_02.Value,'DisplayName',LFMB.pt_igfu_02.Name)

figure()
title([ TOR ' H2 Sleeve'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igfu_05.Value,'DisplayName',LFMB.pt_igfu_05.Name)

figure()
title([ TOR ' H2 Core'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igfu_06.Value,'DisplayName',LFMB.pt_igfu_06.Name)

figure()
title([ TOR ' Igniter Ox'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igox_07.Value,'DisplayName',LFMB.pt_igox_07.Name)

figure()
title([ TOR ' H2 Temp'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Temperature (deg F)')
plot(LFMB.time.Value,LFMB.tc_igfu_01.Value,'DisplayName',LFMB.tc_igfu_01.Name)

figure()
title([ TOR ' O2 Temp'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Temperature (deg F)')
plot(LFMB.time.Value,LFMB.tc_igox_03.Value,'DisplayName',LFMB.tc_igox_03.Name)

% computed values
R_u = 8.314; % J/mol/K

%% OX mdot
Cd = 0.79;
A = pi*(0.032)^2/4; % put in in^2
gamma = 1.4; 
P0 = LFMB.pt_igox_07.Value;
MW = 0.032; % kg/mol
R = R_u/MW; % J/kg/K
T0 = LFMB.tc_igox_03.Value;
mdot_ox = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of ox orifice, assuming ideal gas


% %% H2 mdot
% Cd = 0.85;
% A = pi*(0.032)^2/4; % put in in^2
% gamma = 1.4;
% P0 = LFMB.pt_igfu_02.Value;
% MW = 0.002; % kg/mol
% R = R_u/MW; % J/kg/K
% T0 = LFMB.tc_igfu_01.Value; 
% mdot_h2 = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of h2 core orifice, assuming ideal gas

%% H2 Core mdot
Cd = 0.85;
A = pi*(0.01)^2/4; % put in in^2
gamma = 1.4;
P0 = LFMB.pt_igfu_06.Value;
MW = 0.002; % kg/mol
R = R_u/MW; % J/kg/K
T0 = LFMB.tc_igfu_01.Value; 
mdot_h2_core = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of h2 core orifice, assuming ideal gas

%% H2 Sleeve mdot
Cd = 0.75;
A = pi*(0.063)^2/4;
gamma = 1.4;
P0 = LFMB.pt_igfu_05.Value;
MW = 0.002; % kg/mol
R = R_u/MW; % J/kg/K
T0 = LFMB.tc_igfu_01.Value; 
mdot_h2_sleeve = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of h2 sleeve orifice, assuming ideal gas

%mdot_h2_sleeve = mdot_h2 - mdot_h2_core;


%% Plotting mass flow rates
t = LFMB.time.Value;
figure()
plot(t, mdot_ox, t, mdot_h2_sleeve, 'r', t, mdot_h2_core, 'b') % , t, mdot_h2
title('Mass flow rates')
grid on
legend on
xlabel('Time (s)')
ylabel('Mass flow rate (kg/s)')
legend('O2','H2 sleeve','H2 core', 'H2 Upstream')

%Plot OF ratio
t = LFMB.time.Value;
figure()
OF_tot = mdot_ox ./ (mdot_h2_sleeve + mdot_h2_core);
OF_core = mdot_ox ./ mdot_h2_core;
plot(t, OF_tot, 'r') % , total OF
title('Total OF Ratio')
grid on
legend on
xlabel('Time (s)')
ylabel('OF Ratio')

figure()
plot(t, OF_core, 'b') % , t, core OF
title('Total OF Ratio')
grid on
legend on
xlabel('Time (s)')
ylabel('OF Ratio')

function mdot = compute_mdot(Cd, A, gamma, P0, R, T0)

    % Cd     [unitless]
    % A      [in^2]
    % gamma  [unitless]
    % P0     [psia]
    % R      [J/kg/K]
    % T0     [deg F]

    A = A * (0.0254)^2; % m^2
    P0 = P0 * 6894.76; % Pa
    T0 = (T0 - 32)*5/9 + 273.15; % K 

    C_star_i = sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1))); % unitless
    
    mdot = zeros(length(P0),1);

    for i = 1:length(mdot)
        G_i = C_star_i*(P0(i)/sqrt(R*T0(i))); % kg/s/m^2
        mdot_i = Cd*A*G_i; % kg/s

        mdot(i) = mdot_i;
    end
end