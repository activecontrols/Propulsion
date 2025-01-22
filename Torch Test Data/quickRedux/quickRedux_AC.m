%AC Data Redux Code, based on Zucrow ZL3 Data redux Code
%Authors: Vishal Ravi, Jan Ayala
%TODO: Add a way to mass convert data
%      Automatic upload to sift

clear all; close all

testid = 1;
saveData = true;
TDMSLoc = 'C:\Users\Andrew Radulovich\Documents\GitHub\Propulsion\Torch Test Data\hotfire4'; %where you are saving  directory
saveData =true;
excelsheet = true; % Export data to excel sheet
excelfilename = 'hotfire4_data.xlsx';
TDMSLoc = '..\hotfire4'; %where you are saving  directory
addpath(cd)

dataFileNameMB = 'DataLog_2024-1215-0026-08_PSPAC_Data_Wiring.tdms'; %file name

TOR = 'AC Test 1'; %Label of test   

LFmatFilename = sprintf('Test_%d_Data.mat',testid);

        LFMB = reduceTDMS([TDMSLoc,'\',dataFileNameMB],1,nan);
 
   % Package and Save
    if saveData% && ZeroSave
        fprintf('Saving Low Frequency Data...\n')
        dataNames = fieldnames(LFMB);
        for i = 1:numel(dataNames)
            dataMatrix(:, i) = LFMB.(dataNames{i}).Value;
        end
        dataTable = array2table(dataMatrix, 'VariableNames', dataNames);
        writetable(dataTable, [dataFileNameMB, '.csv']) 
        fprintf('Low Frequency Data Saved.\n\n')
    end
    % 
    % if excelsheet 
    %     time = LFMB.time.Value;
    %     pt_pdox_083 = LFMB.pt_pdox_083.Value;
    %     pt_n2_026 = LFMB.pt_n2_026.Value;
    %     pt_ch_058 = LFMB.pt_ch_058.Value;
    %     pt_n2_076 = LFMB.pt_n2_076.Value;
    %     pt_igox_04 = LFMB.pt_igox_04.Value;
    %     pt_igfu_02 = LFMB.pt_igfu_02.Value;
    %     pt_igfu_05 = LFMB.pt_igfu_05.Value;
    %     pt_igfu_06 = LFMB.pt_igfu_06.Value;
    %     pt_igox_07 = LFMB.pt_igox_07.Value;
    %     pv_pdfu_059_fb = LFMB.pv_pdfu_059_fb.Value;
    %     pv_pdox_087_fb = LFMB.pv_pdox_087_fb.Value;
    %     pv_n2_077_fb = LFMB.pv_n2_077_fb.Value;
    %     tc_igfu_01 = LFMB.tc_igfu_01.Value;
    %     tc_igox_03 = LFMB.tc_igox_03.Value;
    %     T = table(pt_pdox_083,pt_n2_026,pt_ch_058,pt_n2_076,pt_igox_04,pt_igfu_02,pt_igfu_05,pt_igfu_06,pt_igox_07,pv_pdfu_059_fb,pv_pdox_087_fb,pv_n2_077_fb,tc_igfu_01,tc_igox_03);
    %     writetable(T,excelfilename,'Sheet',1,'Range','A1')
    % end 

allData.LFMB = LFMB;


% raw data plots

% figure()
% title([ TOR ' Ox Supply Pressure'])
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Pressure (psia)')
% plot(LFMB.time.Value,LFMB.pt_pdox_083.Value,'k','DisplayName',LFMB.pt_pdox_083.Name)
% 
% figure()
% title([ TOR ' N2 Supply Pressure'])
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Pressure (psia)')
% plot(LFMB.time.Value,LFMB.pt_n2_026.Value,'k','DisplayName',LFMB.pt_n2_026.Name)

% figure()
% title([ TOR ' H2 Supply Pressure'])
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Pressure (psia)')
% plot(LFMB.time.Value,LFMB.pt_ch_058.Value,'k','DisplayName',LFMB.pt_ch_058.Name)
% 
% figure()
% title([ TOR ' Purge Pressure'])
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Pressure (psia)')
% plot(LFMB.time.Value,LFMB.pt_n2_076.Value,'k','DisplayName',LFMB.pt_n2_076.Name)
%%
figure()
title([ TOR ' O2 Upstream Timings'])
hold on
grid on
legend on; xlim([0 3])
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igox_04.Value,'b','DisplayName',LFMB.pt_igox_04.Name); hold on

figure()
title([ TOR ' H2 Upstream Timings'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left; xlim([0 3])
ylabel('Pressure (psia)')
plot(LFMB.time.Value,LFMB.pt_igfu_02.Value,'b','DisplayName',LFMB.pt_igfu_02.Name); hold on
%
figure()
title([ TOR ' H2 Sleeve'])
hold on
grid on
fieldnames(LFMB);legend on
xlabel('Time (s)'); xlim([0 5])
yyaxis left 
ylabel('Pressure (psia)'); xlim([0 5])
plot(LFMB.time.Value,LFMB.pt_igfu_05.Value,'k','DisplayName',LFMB.pt_igfu_05.Name)
plot(LFMB.time.Value,LFMB.pt_igfu_06.Value,'b','DisplayName',LFMB.pt_igfu_06.Name)
plot(LFMB.time.Value,LFMB.pt_igox_07.Value,'r','DisplayName',LFMB.pt_igox_07.Name)
%%
% 
% figure()
% title([ TOR ' H2 Temp'])
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Temperature (deg F)')
% plot(LFMB.time.Value,LFMB.tc_igfu_01.Value,'k','DisplayName',LFMB.tc_igfu_01.Name)
% 
% figure()
% title([ TOR ' O2 Temp'])
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Temperature (deg F)')
% plot(LFMB.time.Value,LFMB.tc_igox_03.Value,'k','DisplayName',LFMB.tc_igox_03.Name)

% computed values
R_u = 8.314; % J/mol/K

% ox mdot
Cd = 0.76;
A = pi*(0.04)^2/4; % put in in^2
gamma = 1.4; 
P0 = LFMB.pt_igox_04.Value;
MW = 0.032; % kg/mol
R = R_u/MW; % J/kg/K
T0 = LFMB.tc_igox_03.Value;
mdot_ox = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of ox orifice, assuming ideal gas

% h2 core mdot
Cd = 0.73;
A = pi*(0.016)^2/4;
gamma = 1.4;
P0 = LFMB.pt_igfu_02.Value;
MW = 0.002; % kg/mol
R = R_u/MW; % J/kg/K
T0 = LFMB.tc_igfu_01.Value; 
mdot_h2_core = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of h2 core orifice, assuming ideal gas

% h2 sleeve mdot
Cd = 0.75;
A = pi*(0.055)^2/4;
gamma = 1.4;
P0 = LFMB.pt_igfu_02.Value;
MW = 0.002; % kg/mol
R = R_u/MW; % J/kg/K
T0 = LFMB.tc_igfu_01.Value; 
mdot_h2_sleeve = compute_mdot(Cd, A, gamma, P0, R, T0); % mdot using sensors upstream of h2 sleeve orifice, assuming ideal gas

% plotting mass flow rates
t = LFMB.time.Value;
figure()
plot(t, mdot_ox,'k', t, mdot_h2_sleeve, 'r', t, mdot_h2_core, 'b')
title('Mass flow rates'); xlim([0 3])
grid on
legend on
xlabel('Time (s)')
ylabel('Mass flow rate (kg/s)')
legend('O2','H2 sleeve','H2 core')

t = LFMB.time.Value;
figure()
plot(t, mdot_ox./(mdot_h2_sleeve + mdot_h2_core))
title('OF Ratio'); xlim([0 3])
grid on
legend on
xlabel('Time (s)')
ylabel('OF')
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
        mdot_i = Cd*A*G_i* 2.20462; % lb/s

        mdot(i) = mdot_i;
    end
end
