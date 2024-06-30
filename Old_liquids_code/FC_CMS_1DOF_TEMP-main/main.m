%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Main Script
Matt Currie, Jonah Fouts, Cameron Williams, Tyler Trostle
%}

clear;
clc;

%% Adds paths to search functions for
functionPath = append(pwd, '/bin');
CEAPath = append(pwd, '/PSP_CEA_function_wrapper');
INPPath = append(pwd, '/INP_OUT');
RefPropPath = append(pwd, '/REFPROP');

addpath(functionPath);
addpath(CEAPath);
addpath(INPPath);
addpath(RefPropPath);

% Number of outputs
outputNum = 44;

%% Runs Loop Function
dataReturn = frootLoops(outputNum);

%% Condenses Data into its own vectors
outerDiameter = dataReturn(:, 1); % [in]
mDry = dataReturn(:, 2); % [lbm]
mWet = dataReturn(:, 3); % [lbm]
mTanks = dataReturn(:, 4); % [lbm]
burnTime = dataReturn(:, 5); % [s]
expectedIsp = dataReturn(:, 6); % [s]
maxAlt = dataReturn(:, 7); % [ft]
railSpeed = dataReturn(:, 8); % [ft/s]
maxMach = dataReturn(:, 9);
maxAccel = dataReturn(:, 10); % [g]
pressureChamber = dataReturn(:, 11); % [psi]
OF = dataReturn(:, 12);
mDot = dataReturn(:, 13); % [lbm/s]
thrust = dataReturn(:, 14); % [lbf]
rocketHeight = dataReturn(:, 15); % [ft]
AR = dataReturn(:, 16);
heightFu = dataReturn(:, 17); % [ft]
heightLOX = dataReturn(:, 18); % [ft]
burnoutAlt = dataReturn(:, 19); % [ft]
tubeIDFu = dataReturn(:, 20); % [in]
tubeIDOx = dataReturn(:, 21); % [in]
expectedCStar = dataReturn(:, 22); % [m/s]
throatAblation = dataReturn(:, 23); % [in]
exitRadius = dataReturn(:, 24); % [in]
twr = dataReturn(:, 25);
rpTankVolume = dataReturn(:, 26); % [in^3]
loxTankVolume = dataReturn(:, 27); % [in^3]
heTankVolume = dataReturn(:, 28); % [L]
exitPresHalfBurn = dataReturn(:, 29); % [psi]
pressureFuelTank = dataReturn(:, 30); % [psi]
pressureOxTank = dataReturn(:, 31); % [psi]
velOx = dataReturn(:, 32); % [ft/s]
velFu = dataReturn(:, 33); % [ft/s]
runLineOx = dataReturn(:, 34); 
runLineFu = dataReturn(:, 35);
pressureExit = dataReturn(:, 36); % [psi]
combustionTemperature = dataReturn(:,37); % [K]
altAtTenSeconds = dataReturn (:, 38); % [ft]
mCopv = dataReturn(:, 39);
radiusThroatInit = dataReturn(:, 40);
dataIndex = dataReturn(:, 41); % index that is for all rockets (good and bad)
fcPercents = dataReturn(:, 42);
fcDiams = dataReturn(:, 43);
fcNHoles = dataReturn(:, 44);

rocketNum = transpose(1 : 1 : size(dataReturn, 1));

%% Clears uneeded points
clear CEAPath currentPath functionPath INPPath

%% Resets settings
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED

%% Saves Workspace
save('Data_Reduc');
