%RocketLab Data Reduction Code

%%%%%%TDMS Reading%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note: TDMS Reading can be commented out after TDMS Data Struc is saved

testid = 1;
saveData =true;
% cd 'Z:\data\Rocket Lab\Testing';
cd 'G:\rocketlab\Loging\quickRedux'
pathname = genpath('Data');
addpath(pathname)

dataFileNameMB = 'DataLog_2023-0324-1410-27_Config_MB_Fac_Fenix_20230303.tdms';

LFmatFilename = sprintf('Test_%d_Data',testid);

        LFMB = reduceTDMS(dataFileNameMB,1,nan);
 
        
%Package and Save
    if saveData% && ZeroSave
        fprintf('Saving Low Frequency Data...\n')
        save([pwd,'\',LFmatFilename],'LFMB')
        fprintf('Low Frequency Data Saved.\n\n')
    end

allData.LFMB = LFMB;
%% 

%%%%%%%Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
title('LOX')
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
% ylabel('Pressure (psi)')
% plot(LFMB.time.Value,LFMB.pt_lox_mb_02.Value,'g','DisplayName',LFMB.pt_lox_mb_02.Name)
% plot(LFMB.time.Value,LFMB.ptf_n2_04.Value,'r','DisplayName',LFMB.ptf_n2_04.Name)
% yyaxis right
% ylabel('Temperature (F)')
% plot(LFMB.time.Value,LFMB.tcg_m_o.Value,'k','DisplayName',LFMB.tcg_m_o.Name)
% plot(LFMB.time.Value,LFMB.tc_lox_mb_01.Value,'b','DisplayName',LFMB.tc_lox_mb_01.Name)
% plot(LFMB.time.Value,LFMB.tcf_ox_04.Value,'m','DisplayName',LFMB.tcf_ox_04.Name)
plot(LFMB.time.Value,LFMB.cv_lox_mb_01_fb.Value,'g','DisplayName',LFMB.cv_lox_mb_01_fb.Name)
plot(LFMB.time.Value,LFMB.cv_lng_mb_01_fb.Value,'r','DisplayName',LFMB.cv_lng_mb_01_fb.Name)




%% 

% figure(2)
% title('LNG')
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis left 
% ylabel('Pressure (psi)')
% plot(LFMB.time.Value,LFMB.pt_lng_mb_02.Value,'g','DisplayName',LFMB.pt_lng_mb_02.Name)
% plot(LFMB.time.Value,LFMB.pt_lng_mb_06.Value,'k','DisplayName',LFMB.pt_lng_mb_06.Name)
% yyaxis right
% ylabel('Temperature (F)')
% plot(LFMB.time.Value,LFMB.tc_lng_mb_05.Value,'r','DisplayName',LFMB.tc_lng_mb_05.Name)
% plot(LFMB.time.Value,LFMB.tcg_m_f.Value,'m','DisplayName',LFMB.tcg_m_f.Name)
% plot(LFMB.time.Value,LFMB.tc_lng_mb_02.Value,'color',[0.8500 0.3250 0.0980],'DisplayName',LFMB.tc_lng_mb_02.Name)

% figure(3)
% title('TEA TEB')
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% yyaxis right
% ylabel('Valve Feedback (V)')
% plot(LFMB.time.Value,LFMB.pv_tt_ifc_01_fb.Value,'DisplayName',LFMB.pv_tt_ifc_01_fb.Name)
% yyaxis left
% ylabel('Pressure (psi)')
% plot(LFMB.time.Value,LFMB.pt_tt_ifc_03.Value,'DisplayName',LFMB.pt_tt_ifc_03.Name)
% % plot(LFMB.time.Value,LFMB.pt_ta_vent_01.Value,'DisplayName',LFMB.pt_ta_vent_01.Name)

% figure(4)
% title('Chamber Pressures')
% hold on
% grid on
% legend on
% xlabel('Time (s)')
% ylabel('Pressure (psi)')
% plot(LFMB.time.Value,LFMB.sp_m_o.Value,'DisplayName',LFMB.sp_m_o.Name)
% plot(LFMB.time.Value,LFMB.sp_m_f.Value,'DisplayName',LFMB.sp_m_f.Name)
% plot(LFMB.time.Value,LFMB.sp_cc_a1.Value,'r-','DisplayName',LFMB.sp_cc_a1.Name)
% plot(LFMB.time.Value,LFMB.sp_ch_a1.Value,'m-','DisplayName',LFMB.sp_ch_a1.Name)
