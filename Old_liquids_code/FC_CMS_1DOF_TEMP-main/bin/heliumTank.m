function [mHelium, volumeHeTank] = heliumTank(pressureFuTank, pressureOxTank, volumeOxTank, volumeFuTank)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Helium Sizing
Matt Currie

Input variables:  pressureFuTank [psi]
                  pressureOxTank [psi]
                  volumeOxTank [in^3]
                  volumeFuTank [in^3]

Output variables: mHelium [lbm]
                  volumeHeTank [in^3]
%}

% Intialized Constants
oversizeFactor = 1.1;
heSpecificHeatRatio = 1.66;
isentropicExponent = (heSpecificHeatRatio - 1) / heSpecificHeatRatio;

% Conversion Factors
cubicMToCubicIn = 61024;
cubicInToCubicM = 1 / cubicMToCubicIn;
psiTokPa = 6.89476;
kilogramToPoundMass = 2.2046;

% Conversiions
heGasConstant = 2.0771; % [kJ/kg-K]
pressureHelium = 4500 * psiTokPa; % PSI to kPa
pressureFuTank = pressureFuTank * psiTokPa; % PSI to kPa
pressureOxTank = pressureOxTank * psiTokPa; % PSI to kPa
volumeOxTank = volumeOxTank * cubicInToCubicM; % in^3 to m^3
volumeFuTank = volumeFuTank  * cubicInToCubicM; % in^3 to m^3

% Define end temperatures
ambTemp = 273.15 + 22; % Ambient temp [K]
tempEndOxTank = (ambTemp) / ((pressureHelium / pressureOxTank) ^ isentropicExponent); % [K]
tempEndFuTank = (ambTemp) / ((pressureHelium / pressureFuTank) ^ isentropicExponent); % [K]
tempEndHe = (ambTemp) / ((pressureHelium / pressureFuTank) ^ isentropicExponent); % [K]

% Calculate mass of helium in each tank
mHeliumOxTank = (pressureOxTank * volumeOxTank) / (heGasConstant * tempEndOxTank); % Calculate mass of helium in the LOX tank
mHeliumFuTank = (pressureFuTank * volumeFuTank) / (heGasConstant * tempEndFuTank); % Calculate mass of helium in the RP1 tank

% Calculate volume of helium tank and total mass of helium
volumeHeTank = oversizeFactor * heGasConstant * (mHeliumOxTank + mHeliumFuTank) / ((pressureHelium / ambTemp) - (pressureFuTank / tempEndHe)); % [m^3]
mHelium = (pressureHelium * volumeHeTank) / (heGasConstant * ambTemp); % [kg]

volumeHeTank = volumeHeTank  * cubicMToCubicIn; % m^3 to in^3
mHelium = mHelium * kilogramToPoundMass; % kg to lbs
end