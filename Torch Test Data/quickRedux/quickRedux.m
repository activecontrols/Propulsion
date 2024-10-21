close all

testid = 1;
saveData =true;
tdmsLoc = 'Y:\rocketlab\Loging\quickRedux\'; %where you are saving  directory
addpath(cd)

dataFileNameMB = 'DataLog_2024-0412-1610-22_Config_MB_Fac_RPS_Preburner_20240223.tdms'; %file name

TOR = 'HF19'; %Label of test

LFmatFilename = sprintf('Test_%d_Data',testid);

        LFMB = reduceTDMS(dataFileNameMB,1,nan);
 
   % Package and Save
    if saveData% && ZeroSave
        fprintf('Saving Low Frequency Data...\n')
        save([pwd,'\',LFmatFilename],'LFMB')
        fprintf('Low Frequency Data Saved.\n\n')
    end

allData.LFMB = LFMB;

%%%%%%%Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure()
title([ TOR ' LOX - Pressures'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psi)')
plot(LFMB.time.Value,LFMB.ptf_n2_04.Value,'k','DisplayName',LFMB.ptf_n2_04.Name)
plot(LFMB.time.Value,LFMB.pt_lox_mb_02.Value,'DisplayName',LFMB.pt_lox_mb_02.Name)
plot(LFMB.time.Value,LFMB.pt_ta_01.Value,'g-','DisplayName',LFMB.pt_ta_01.Name)
plot(LFMB.time.Value,LFMB.pt_ta_03.Value,'r-.','DisplayName',LFMB.pt_ta_03.Name)
saveas(gcf, [TOR  '-LOXpress'] );

figure()
title([TOR ' LNG - Pressures'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psi)')
plot(LFMB.time.Value,LFMB.pt_lng_mb_06.Value,'k','DisplayName',LFMB.pt_lng_mb_06.Name)
plot(LFMB.time.Value,LFMB.pt_lng_mb_02.Value,'DisplayName',LFMB.pt_lng_mb_02.Name)
plot(LFMB.time.Value,LFMB.pt_ta_02.Value,'g-','DisplayName',LFMB.pt_ta_02.Name)
plot(LFMB.time.Value,LFMB.pt_ta_03.Value,'r-.','DisplayName',LFMB.pt_ta_03.Name)
saveas(gcf, [ TOR '-LNGpress'] );

figure()
title([TOR ' Differential - Pressures'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left 
ylabel('Pressure (psi)')
plot(LFMB.time.Value,(LFMB.pt_ta_01.Value - LFMB.pt_ta_03.Value),'b','DisplayName','Ta-01 - Ta-03')
plot(LFMB.time.Value,(LFMB.pt_ta_02.Value - LFMB.pt_ta_03.Value),'r','DisplayName','Ta-02 - Ta-03')
saveas(gcf, [TOR '-Diffpress'] );

figure()
title([TOR ' LNG - Temperatures'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left
ylabel('Temperatures (F)')
plot(LFMB.time.Value,LFMB.tc_lng_mb_02.Value,'DisplayName',LFMB.tc_lng_mb_02.Name)
plot(LFMB.time.Value,LFMB.tc_ta_02.Value,'g','DisplayName',LFMB.tc_ta_02.Name)
plot(LFMB.time.Value,LFMB.tc_ta_03.Value,'r','DisplayName',LFMB.tc_ta_03.Name)
saveas(gcf, [TOR  '-LNGtemp']);

figure()
title([TOR ' LOX - Temperatures'])
hold on
grid on
legend on
xlabel('Time (s)')
yyaxis left
ylabel('Temperatures (F)')
plot(LFMB.time.Value,LFMB.tc_lox_mb_01.Value,'DisplayName',LFMB.tc_lox_mb_01.Name)
plot(LFMB.time.Value,LFMB.tc_lox_mb_02.Value,'DisplayName',LFMB.tc_lox_mb_02.Name)
plot(LFMB.time.Value,LFMB.tc_ta_01.Value,'g','DisplayName',LFMB.tc_ta_01.Name)
plot(LFMB.time.Value,LFMB.tc_ta_03.Value,'r','DisplayName',LFMB.tc_ta_03.Name)
saveas(gcf, [TOR '-LOXtemps']);

figure
title([TOR ' Valve Feedback'])
hold on
grid on
legend on
xlabel('Time (s)')
ylabel('Feedback')
plot(LFMB.time.Value,LFMB.cv_lox_mb_01_fb.Value,'DisplayName',LFMB.cv_lox_mb_01_fb.Name)
plot(LFMB.time.Value,LFMB.cv_lng_mb_02_fb.Value,'DisplayName',LFMB.cv_lox_mb_01_fb.Name)


tc_lng_mb_02 = LFMB.tc_lng_mb_02.Value;
tc_ta_02 = LFMB.tc_ta_02.Value;
tc_ta_03 =LFMB.tc_ta_03.Value;
