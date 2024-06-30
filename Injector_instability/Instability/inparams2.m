% This function takes the string information from the GUI and converts
% to a real number. It also performs basic calculations.
function inparams2(hObject, eventdata, s1)
global a b cstar dcdMR m mf mo MR f fmax pc pod pfd
global tau_f tau_o theta_g xmax zmax gc df
global LoO CmO LoF CmF
%%% Define variables from GUI inputs
in_str1 = get(s1.edit1, 'string');
MWg = str2num(in_str1);
in_str2 = get(s1.edit2, 'string');
cstar = str2num(in_str2)*12; % put in in/sec
in_str3 = get(s1.edit3, 'string');
pc = str2num(in_str3);
in_str4 = get(s1.edit4, 'string');
tau_o = str2num(in_str4)/1000; % put in seconds
in_str5 = get(s1.edit5, 'string');
dcdMR = str2num(in_str5)*12; % put in in/sec-MR
in_str6 = get(s1.edit6, 'string');
Tc = str2num(in_str6);
in_str7 = get(s1.edit7, 'string');
tau_f = str2num(in_str7)/1000; % put in seconds
in_str8 = get(s1.edit8, 'string');
a = str2num(in_str8);
in_str9 = get(s1.edit9, 'string');
mo = str2num(in_str9);
in_str10 = get(s1.edit10, 'string');
pod = str2num(in_str10);
in_str11 = get(s1.edit11, 'string');
b = str2num(in_str11);
in_str12 = get(s1.edit12, 'string');
mf = str2num(in_str12);
in_str13 = get(s1.edit13, 'string');
pfd = str2num(in_str13);
in_str14 = get(s1.edit14, 'string');
Lc = str2num(in_str14);
in_str15 = get(s1.edit15, 'string');
Dc = str2num(in_str15);
in_str16 = get(s1.edit16, 'string');
Lt = str2num(in_str16);
in_str17 = get(s1.edit17, 'string');
Dt = str2num(in_str17);
in_str18 = get(s1.edit18, 'string');
xmax = str2num(in_str18);
in_str19 = get(s1.edit19, 'string');
fmax = str2num(in_str19);
in_str20 = get(s1.edit20, 'string');
zmax = str2num(in_str20);
in_str21 = get(s1.edit21, 'string');
df = str2num(in_str21);
in_str22 = get(s1.edit22, 'string');
LoO = str2num(in_str22);
in_str23 = get(s1.edit23, 'string');
CmO = str2num(in_str23);
in_str24 = get(s1.edit24, 'string');
LoF = str2num(in_str24);
in_str25 = get(s1.edit25, 'string');
CmF = str2num(in_str25);
%%% Constants
gc=386.0874; % Grav. Constant of Proportionality [lbm*in/lbf-s^2]
Ru=18544.8; % Universal Gas Constant [in*lbf/lbmol-degR]
%%% Basic Calculations
m=mo+mf; % Total Flow Rate [lbm/sec]
MR=mo/mf; % Mixture Ratio [-]
At=pi*Dt^2/4; % Throat Area [in^2]
Vc=pi*Dc^2*Lc/4+pi*Lt*(Dt^2+Dc^2+Dc*Dt)/12; % Chamber Volume [in^3]
Rg=Ru/MWg; % Gas Constant for Combustion Gas
theta_g=Vc*cstar/(At*gc*Rg*Tc); % Gas Residence Time [sec]
%%% Define Variable value ranges
f=-fmax:df:fmax; % Frequency Range
%%% Define Input File Layout and save into excel
inpLayout={'Variable' 'Description' 'Units' 'Value';
'Dc' 'Chamber Diameter' 'in' Dc;
'Lc' 'Chamber Length' 'in' Lc;
'Dt' 'Throat Diameter' 'in' Dt;
'Lt' 'Convergent Section Length' 'in' Lt;
' ' ' ' ' ' ' ';
'pod' 'Ox. Pressure Drop' 'psia' pod;
'mo' 'Ox. Flow Rate' 'lbm/sec' mo;
'a' 'Ox. Exponent' '--' a;
'pfd' 'Fuel Pressure Drop' 'psia' pfd;
'mf' 'Fuel Flow Rate' 'lbm/sec' mf;
'b' 'Fuel Exponent' '--' b;
' ' ' ' ' ' ' ';
'tau_o' 'Ox. Time Lag' 'msec' tau_o*1000;
'pc' 'Chamber Pressure' 'psia' pc;
'cstar' 'c*' 'ft/sec' cstar/12;
'MWg' 'Gas MW' 'lbm/lbmol' MWg;
'tau_f' 'Fuel Time Lag' 'msec' tau_f*1000;
'Tc' 'Chamber Temperature' 'deg R' Tc;
'dcdMR' 'Slope of c*(MR)' 'ft/sec-MR' dcdMR/12;
' ' ' ' ' ' ' ';
'LoO' 'Ox. Injector Inertance' 'lbf*s^2/lbm-in^2' LoO;
'CmO' 'Ox. Injector Compliance' 'lbm*in^2/lbf' CmO;
'LoF' 'Fuel Injector Inertance' 'lbf*s^2/lbm-in^2' LoF;
'CmF' 'Fuel Injector Compliance' 'lbm*in^2/lbf' CmF;
' ' ' ' ' ' ' ';
'fmax' 'Max Frequency' 'Hz' fmax;
'xmax' 'Dpo/Pc Axis Max' '--' xmax;
'df' 'Frequency Resolution' 'Hz' df;
'zmax' 'Dpf/Pc Axis Max' '--' zmax;
};
xlswrite('CSCinput.xls',inpLayout,'Sheet1','A1');
mainfunc(hObject, eventdata, s1) % Call mainfunc