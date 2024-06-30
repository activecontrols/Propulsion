% This function is accessed only if the 'Use Input File' button is pressed.
% It reads the text file, converts the strings to numbers, and write the
% values into the edit windows in the GUI.
function inparams1(hObject, eventdata, s1)
global fileinp
fileinp=1; % Flag indicating inputs are taken from input file.
%%% read floating point value from Excel
Dc=xlsread('CSCinput.xls','Sheet1','D2');
Lc=xlsread('CSCinput.xls','Sheet1','D3');
Dt=xlsread('CSCinput.xls','Sheet1','D4');
Lt=xlsread('CSCinput.xls','Sheet1','D5');
pod=xlsread('CSCinput.xls','Sheet1','D7');
mo=xlsread('CSCinput.xls','Sheet1','D8');
a=xlsread('CSCinput.xls','Sheet1','D9');
pfd=xlsread('CSCinput.xls','Sheet1','D10');
mf=xlsread('CSCinput.xls','Sheet1','D11');
b=xlsread('CSCinput.xls','Sheet1','D12');
tau_o=xlsread('CSCinput.xls','Sheet1','D14');
pc=xlsread('CSCinput.xls','Sheet1','D15');
cstar=xlsread('CSCinput.xls','Sheet1','D16');
MWg=xlsread('CSCinput.xls','Sheet1','D17');
tau_f=xlsread('CSCinput.xls','Sheet1','D18');
256
Tc=xlsread('CSCinput.xls','Sheet1','D19');
dcdMR=xlsread('CSCinput.xls','Sheet1','D20');
LoO=xlsread('CSCinput.xls','Sheet1','D22');
CmO=xlsread('CSCinput.xls','Sheet1','D23');
LoF=xlsread('CSCinput.xls','Sheet1','D24');
CmF=xlsread('CSCinput.xls','Sheet1','D25');
fmax=xlsread('CSCinput.xls','Sheet1','D27');
xmax=xlsread('CSCinput.xls','Sheet1','D28');
df=xlsread('CSCinput.xls','Sheet1','D29');
zmax=xlsread('CSCinput.xls','Sheet1','D30');
%%% Convert number to string so it can be displayed in GUI
Dc_str=num2str(Dc);
Lc_str=num2str(Lc);
Dt_str=num2str(Dt);
Lt_str=num2str(Lt);
pod_str=num2str(pod);
mo_str=num2str(mo);
a_str=num2str(a);
pfd_str=num2str(pfd);
mf_str=num2str(mf);
b_str=num2str(b);
tau_o_str=num2str(tau_o);
pc_str=num2str(pc);
cstar_str=num2str(cstar);
MWg_str=num2str(MWg);
tau_f_str=num2str(tau_f);
Tc_str=num2str(Tc);
dcdMR_str=num2str(dcdMR);
LoO_str=num2str(LoO);
CmO_str=num2str(CmO);
LoF_str=num2str(LoF);
CmF_str=num2str(CmF);
fmax_str=num2str(fmax);
xmax_str=num2str(xmax);
df_str=num2str(df);
zmax_str=num2str(zmax);
%%% Rename Editable Text Areas on the GUI
% % Inputs taken from input file; disable edit value capability
% Gas MW (lbm/lbmol)
set(s1.edit1, 'string', MWg_str,'Enable', 'off');
% c* (ft/sec)
set(s1.edit2, 'string', cstar_str,'Enable', 'off');
% Chamber Pressure (psia)
set(s1.edit3, 'string', pc_str,'Enable', 'off');
% Oxidizer Time Lag (msec)
set(s1.edit4, 'string', tau_o_str,'Enable', 'off');
% Slope of c*(MR) (ft/sec-MR)
set(s1.edit5, 'string', dcdMR_str,'Enable', 'off');
% Chamber Temperature (deg R)
set(s1.edit6, 'string', Tc_str,'Enable', 'off');
% Fuel Time Lag (msec)
set(s1.edit7, 'string', tau_f_str,'Enable', 'off');
% Oxidizer Exponent, a (mo=Cd*Dpo^a)
set(s1.edit8, 'string', a_str,'Enable', 'off');
% Oxidizer Flow Rate (lbm/sec)
set(s1.edit9, 'string', mo_str,'Enable', 'off');
% Oxidizer Pressure Drop (psia)
set(s1.edit10, 'string', pod_str,'Enable', 'off');
% Fuel Exponent, b (mf=Cd*Dpf^b)
set(s1.edit11, 'string', b_str,'Enable', 'off');
% Fuel Flow Rate (lbm/sec)
set(s1.edit12, 'string', mf_str,'Enable', 'off');
% Fuel Pressure Drop (psia)
set(s1.edit13, 'string', pfd_str,'Enable', 'off');
% Cylinder Length (in)
set(s1.edit14, 'string', Lc_str,'Enable', 'off');
% Chamber Diameter (in)
set(s1.edit15, 'string', Dc_str,'Enable', 'off');
% Convergent Section Length (in)
set(s1.edit16, 'string', Lt_str,'Enable', 'off');
% Throat Diameter (in)
set(s1.edit17, 'string', Dt_str,'Enable', 'off');
% Dpo/Pc Axis Max
set(s1.edit18, 'string', xmax_str,'Enable', 'off');
% Maximum Frequency (Hz)
set(s1.edit19, 'string', fmax_str,'Enable', 'off');
% Dpf/Pc Axis Max
set(s1.edit20, 'string', zmax_str,'Enable', 'off');
% Frequency Resolution (Hz)
set(s1.edit21, 'string', df_str,'Enable', 'off');
% Ox. Injector Inertance
set(s1.edit22, 'string', LoO_str,'Enable', 'off');
% Ox. Injector Compliance
257
set(s1.edit23, 'string', CmO_str,'Enable', 'off');
% Fuel Injector Inertance
set(s1.edit24, 'string', LoF_str,'Enable', 'off');
% Fuel Injector Compliance
set(s1.edit25, 'string', CmF_str,'Enable', 'off');