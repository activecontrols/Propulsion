%
% Concentrated Combustion Stability Code v1.3
%
% This program solves Crocco's concentrated combustion model in dimensional
% form. The solution was derived from chapter 3 of Crocco's "Theory of
% Combustion Instability in Liquid Propellant Rocket Motors" book. The
% program also solves a concentrated combustion model commonly used in
% software tools today. The program also solves a concentrated combustion
% model that incorporates injector admittance boundary conditions.
% The program also solves a new concentrated combustion model that
% incorporates damping aspects.
%
%%%% Model Descriptions %%%%
% Crocco Model: Derived from Chapter 3 of Crocco's 'Rocket Instability'
% Book, 1956
%
% Contemporary Model: Crocco Model with injection response and double
% time lag
%
% Modified Crocco Model: Contemporary Model with injector boundary
% conditions in chamber response -ref Casiano PhD dissertation
%
% Damping Model: Contemporary Model with injector boundary conditions in
% chamber response and chamber damping -ref Casiano PhD dissertation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % Clear all variables from workspace
format('long')
%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%
%%% Constants
gc=386.0874; % Grav. Constant of Proportionality [lbm*in/lbf-s^2]
%%% Define Variable value ranges
df=1.005; % Frequency Resolution [Hz]
fmin=0; % Minimum Frequency [Hz]
fmax=1000; % Maximum Frequency [Hz]
f=fmin:df:fmax; % Frequency Range [Hz]
omega=2*pi*f; % Circular Frequency [rad/sec]
s=i*omega; % Laplace Variable (s=a+i*omega, a=0)[rad/sec]
% Inputs for Concentrated Combustion Models
Lc=59.9; % Length of Chamber Cylindrical Section [in]
Dc=3.5; % Chamber Diameter [in]
Ac=pi*Dc^2/4; % Chamber x-sec. area [in2]
pc=1119.7; % Chamber Pressure [psia]
mO=6.638; % Fuel Flow Rate [lbm/sec]
mF=9.064; % Oxidizer Flow Rate [lbm/sec]
m=mO+mF; % Total Flow Rate [lbm/sec]
% Additional Parameters for Crocco Concentrated Combustion
c=5030*12; % Speed of Sound [in/s]
u=872.2*12; % Gas Velocity [in/s]
rho=0.00016729; % Gas Density [lbm/in^3]
rhog=rho/gc; % Gas Density [lbf-s^2/in^4]
gam=1.3778; % Specific heat ratio []
Me=u/c; % Mach number []
% Additional Parameters for Damping Model
k=omega/c; % Wave number [rad/in]
lambda=2*rho*2554*omega.*omega.^(-1.19); % Damping Parameter [lbm/in3-s]
% Time lag Model Inputs
tauO=0.002; % Total Time lag - oxidizer [sec]
tauF=0.002; % Total Time lag - fuel [sec]
taus=0.0013; % Sensitive time lag [sec]
n=2; % Interaction index []
psi=0.05*Lc; % Concentrated combustion axial location [in]
%%%%%%%%%%%%%%%%%%%%% INECTOR ADMITTANCE %%%%%%%%%%%%%%%%%%
% Manifold
rhoMO=68.233/1728; % Ox Manifold Fluid Density [lbm/in^3]
aMO=2753.5*12; % Ox Manifold Fluid Soundspeed [in/sec]
VMO=7.37; % Ox Manifold Volume [in^3]
KMO=aMO^2*rhoMO/gc; % Ox Manifold Fluid Bulk Modulus [lbf/in^2]
CmO=rhoMO*VMO/KMO; % Ox Manifold Compliance [lbm*in^2/lbf]
Mo=52.1*12/aMO; % Ox Injector Fluid Mach # []
rhoMF=3.5812/1728; % Fuel Manifold Fluid Density [lbm/in^3]
aMF=3558.6*12; % Fuel Manifold Fluid Soundspeed [in/sec]
VMF=50.85; % Fuel Manifold Volume [in^3]
KMF=aMF^2*rhoMF/gc; % Fuel Manifold Fluid Bulk Modulus [lbf/in^2]
CmF=rhoMF*VMF/KMF; % Fuel Manifold Compliance [lbm*in^2/lbf]
Mf=897*12/aMF; % Fuel Injector Fluid Mach # []
% Injector
pox=1264; % Ox Manifold Pressure [psia]
DpO=pox-pc; % Ox Injector Pressure Drop [lbf/in^2]
loO=1.253; % Ox Injector Length [in]
AoO=0.137892; % Ox Injector Total Open Area [in^2]
RoO=2*DpO/mO; % Ox Injector Resistance [lbf*s/lbm-in^2]
LoO=loO/AoO/gc; % Ox Injector Inertance [lbf*s^2/lbm-in^2]
pf=1566.1; % Fuel Manifold Pressure [psia]
DpF=pf-pc; % Fuel Injector Pressure Drop [lbf/in^2]
loF=0.593; % Fuel Injector Length [in]
AoF=0.47416; % Fuel Injector Total Open Area [in^2]
RoF=2*DpF/mF; % Fuel Injector Resitance [lbf*s/lbm-in^2]
LoF=loF/AoF/gc; % Fuel Injector Inertance [lbf*s^2/lbm-in^2]
GO=-1./(RoO+LoO*s+1./(CmO*s)); % Ox Inj. Admittance (no rho') [lbm*in^2/lbf-s]
GF=-1./(RoF+LoF*s+1./(CmF*s)); % Fuel Inj. Admittance (no rho') [lbm*in^2/lbf-s]
%%%%%%%%%%%%%%%%%%%%% NOZZLE ADMITTANCE %%%%%%%%%%%%%%%%%%
Ae=Ac; % Cham. Area at Converging Entrance [in^2]
Ge=(gam+1)/(2*gam)*m/pc; % Nozzle Adm. (m'/p') [lbm*in^2/lbf-s]
alpha_e=-(Ge*pc-m)/(2*Ge*pc-m); % Rearrange equation for Crocco model
%%%%%%%%%%%%%%%%%%%%% MODELS %%%%%%%%%%%%%%%%%%
%%% Characteristic Equation for Crocco Model
B=(alpha_e*Me+1)/(alpha_e*Me-1);
% F0=((B*exp(2*s*(psi-Lc)/c/(Me^2-1))+1)./(B*exp(2*s*(psi-Lc)/c/(Me^2-1))-1)...
% +(-exp(2*s*psi/c/(Me^2-1))+1)./(exp(2*s*psi/c/(Me^2-1))+1)...
% +Me*(1-c^2*rhog/pc*n*(1-exp(-s*taus))))*pc/u/rhog/c;
%%% Characteristic Equation for Contemporary Model
B=(alpha_e*Me+1)/(alpha_e*Me-1);
% F1=((B*exp(2*s*(psi-Lc)/c/(Me^2-1))+1)./(B*exp(2*s*(psi-Lc)/c/(Me^2-1))-1)...
% +(-exp(2*s*psi/c/(Me^2-1))+1)./(exp(2*s*psi/c/(Me^2-1))+1)...
% +Me*(1-c^2*rhog/pc*n*(1-exp(-s*taus))-u*rhog*c*pc/m*(GO.*exp(-s*tauO)...
% +GF.*exp(-s*tauF))))*pc/u/rhog/c;
%%% Characteristic Equation for Modified Crocco Model
Bi=((Mo/aMO*AoO*rhoMF+Mf/aMF*AoF*rhoMO-rhoMO*rhoMF/c/rho*Ac)*gc-rhoMF*GO-rhoMO*GF)...
./((Mo/aMO*AoO*rhoMF+Mf/aMF*AoF*rhoMO+rhoMO*rhoMF/c/rho*Ac)*gc-rhoMF*GO-rhoMO*GF);
Be=(m*c/gam/pc-Ae*gc-Ge*c)./(m*c/gam/pc+Ae*gc-Ge*c);
% F2=((Be.*exp(2*s*(psi-Lc)/c/(Me^2-1))+1)./(Be.*exp(2*s*(psi-Lc)/c/(Me^2-1))-1)...
% -(Bi.*exp(2*s*psi/c/(Me^2-1))+1)./(Bi.*exp(2*s*psi/c/(Me^2-1))-1))*pc/u/rhog/c...
% +(pc/rhog/c^2-n*(1-exp(-s*taus))-pc/m*(GO.*exp(-s*tauO)+GF.*exp(-s*tauF)));
%%% Characteristic Equation for Damping Model
aK=Me^2-1;
bK=-(3/2*Me*lambda.*k/rho./s+2*Me*k);
cK=k.^2+k.^2.*lambda./s/rho;
K1=(-bK+sqrt(bK.^2-4*aK*cK))/2/aK;
K2=(-bK-sqrt(bK.^2-4*aK*cK))/2/aK;
A2=1/rhog/c*(lambda*Me/2+rho*K1.*s./k)./(lambda+rho*(1-Me*K1./k).*s);
B2=1/rhog/c*(lambda*Me/2+rho*K2.*s./k)./(lambda+rho*(1-Me*K2./k).*s);
Hi=-((Mo/aMO*AoO*rhoMF+Mf/aMF*AoF*rhoMO+rhoMO*rhoMF/gc*Ac*B2)*gc-rhoMF*GO-rhoMO*GF)...
./((Mo/aMO*AoO*rhoMF+Mf/aMF*AoF*rhoMO+rhoMO*rhoMF/gc*Ac*A2)*gc-rhoMF*GO-rhoMO*GF);
He=-(m/gam/pc+rho*Ae*B2-Ge)./(m/gam/pc+rho*Ae*A2-Ge).*exp(Lc./omega.*s.*(K1-K2));
% F3=(((He.*A2.*exp(-psi./omega.*s.*(K1-K2))+B2)./(He.*exp(-psi./omega.*s.*(K1-K2))+1)...
% -(Hi.*A2.*exp(-psi./omega.*s.*(K1-K2))+B2)./(Hi.*exp(-psi./omega.*s.*(K1-K2))+1)...
% +(Me/rhog/c-u*n/pc*(1-exp(-s*taus))-u/m*(GO.*exp(-s*tauO)+GF.*exp(-s*tauF))))*rhog*c)*pc/u/rhog/c;
%%%%%%%%%%%%%%%%%%% SYSTEM ADMITTANCES %%%%%%%%%%%%%%%%%%
%Crocco Model
Gc0=((B*exp(2*s*(psi-Lc)/c/(Me^2-1))+1)./(B*exp(2*s*(psi-Lc)/c/(Me^2-1))-1)...
+(-exp(2*s*psi/c/(Me^2-1))+1)./(exp(2*s*psi/c/(Me^2-1))+1))*pc/u/rhog/c+Me*pc/u/rhog/c;
Gj0=0;
Gb0=n*(1-exp(-s*taus));
%Contemporary Model
Gc1=((B*exp(2*s*(psi-Lc)/c/(Me^2-1))+1)./(B*exp(2*s*(psi-Lc)/c/(Me^2-1))-1)...
+(-exp(2*s*psi/c/(Me^2-1))+1)./(exp(2*s*psi/c/(Me^2-1))+1))*pc/u/rhog/c+Me*pc/u/rhog/c;
Gj1=pc/m*(GO.*exp(-s*tauO)+GF.*exp(-s*tauF));
Gb1=n*(1-exp(-s*taus));
%Modified Crocco Model
Gc2=((Be.*exp(2*s*(psi-Lc)/c/(Me^2-1))+1)./(Be.*exp(2*s*(psi-Lc)/c/(Me^2-1))-1)...
-(Bi.*exp(2*s*psi/c/(Me^2-1))+1)./(Bi.*exp(2*s*psi/c/(Me^2-1))-1))*pc/u/rhog/c+pc/rhog/c^2;
Gj2=pc/m*(GO.*exp(-s*tauO)+GF.*exp(-s*tauF));
Gb2=n*(1-exp(-s*taus));
%Damping Model
Gc3=((He.*A2.*exp(-psi./omega.*s.*(K1-K2))+B2)./(He.*exp(-psi./omega.*s.*(K1-K2))+1)...
-(Hi.*A2.*exp(-psi./omega.*s.*(K1-K2))+B2)./(Hi.*exp(-psi./omega.*s.*(K1-K2))+1))*pc/u+pc/rhog/c^2;
Gj3=pc/m*(GO.*exp(-s*tauO)+GF.*exp(-s*tauF));
Gb3=n*(1-exp(-s*taus));
%%%%%%%%%%%%%%%%%%% SYSTEM GAIN AND PHASE %%%%%%%%%%%%%%%%%%
Gain0=abs(-(Gj0+Gb0)./Gc0);
Phase0=angle(-(Gj0+Gb0)./Gc0)*180/pi;
Gain1=abs(-(Gj1+Gb1)./Gc1);
Phase1=angle(-(Gj1+Gb1)./Gc1)*180/pi;
Gain2=abs(-(Gj2+Gb2)./Gc2);
Phase2=angle(-(Gj2+Gb2)./Gc2)*180/pi;
Gain3=abs(-(Gj3+Gb3)./Gc3);
Phase3=angle(-(Gj3+Gb3)./Gc3)*180/pi;
%%%%%%%%%%%%%%%%%%% BODE PLOTS %%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1);
h1=plot(f,Gain0,'k-',f,Gain1,'b--',f,Gain2,'r:',f,Gain3,'g:','LineWidth',3);
axis([0 fmax 0 5])
grid on
title('System Bode Plot')
xlabel('Frequency [Hz]')
ylabel('Open-Loop Gain [-]')
legend('Crocco Model','Contemporary Model','Modified Crocco Model','Damping Model','Location','NorthWest')
subplot(2,1,2);
h2=plot(f,Phase0,'k-',f,Phase1,'b--',f,Phase2,'r:',f,Phase3,'g:','LineWidth',3);
xlabel('Frequency [Hz]')
ylabel('Open-Loop Phase [deg]')
grid on