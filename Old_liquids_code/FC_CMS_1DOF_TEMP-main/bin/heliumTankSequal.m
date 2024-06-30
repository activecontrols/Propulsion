function [mHelium, bestCopvIdx] = heliumTankSequal(pressureFuTank, pressureOxTank, mDotFu, mDotOx, burnTime, pressureService, copvVolume, mCopv)
% function [mHelium, volumeHeTank] = heliumTankSequal(pressureFuTank, pressureOxTank, volumeOxTank, volumeFuTank, mDotFu, mDotOx, burnTime)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Helium Sizing
Matt Currie

Input variables:  pressureFuTank [psi]
                  pressureOxTank [psi]
                  volumeOxTank [in^3]
                  volumeFuTank [in^3]
                  mDotFu [lbm/s]
                  mDotOx [lbm/s]
                  burnTime [s]

Output variables: mHelium [lbm]
                  volumeHeTank [in^3]
%}

% Intialized Constants
oversizeFactor = 1.1;
heSpecificHeatRatio = 1.66;
isentropicExponent = (heSpecificHeatRatio - 1) / heSpecificHeatRatio;
densLOXMetric = 1.014; % [g/cm^3]
% densJetAMetric = 0.804; % density of JetA [g/cm^3]
% densJetAMetric = 0.85131; % density of Ethanol (25% H2O) [g/cm^3]
densJetAMetric = 0.8; % density of Ethanol (5% H2O) [g/cm^3]
deltaT = 0.01;

% Conversion Factors
cubicMToCubicIn = 61024;
cubicInToCubicM = 1 / cubicMToCubicIn;
cubicCmToCubicM = 1 / (100^3);
cubicMtoCubicCm = 100^3;
gramToKg = 1 / 1000;
kgToGram = 1000;
psiTokPa = 6.89476;
kilogramToPoundMass = 2.2046;
litersToCubicM = 0.001;

% Conversions
heGasConstant = 2.0771; % [kJ/kg-K]
pressureFuTank = pressureFuTank * psiTokPa; % PSI to kPa
pressureOxTank = pressureOxTank * psiTokPa; % PSI to kPa
% volumeOxTank = volumeOxTank * cubicInToCubicM; % in^3 to m^3
% volumeFuTank = volumeFuTank  * cubicInToCubicM; % in^3 to m^3
densLOX = densLOXMetric * gramToKg / cubicCmToCubicM; % LOX density [kg/m^3]
densJetA = densJetAMetric * gramToKg / cubicCmToCubicM; % Jet-A density [kg/m^3]
mDotFu = mDotFu / kilogramToPoundMass; % Mass flow rate of fuel [kg/s]
mDotOx = mDotOx / kilogramToPoundMass; % Mass flow rate of oxidizer [kg/s]

% Define end temperatures
ambTemp = 273.15 + 22; % Ambient temp [K]

% Calculate volumetric flow rate out of tanks
volumeDotFu = mDotFu / densJetA; % Volumetric Flow Rate Fuel [m^3/s]
volumeDotOx = mDotOx / densLOX; % Volumetric Flow Rate Oxidizer [m^3/s]
volumeDotHe = volumeDotOx + volumeDotFu;

bestCopvIdx = 0;
mBestCopv = intmax;
bestCopvVolume = intmax;

for tankIdx = 1 : length(pressureService)
    % Set COPV volume and pressure based on available COPVs
    pressureHelium = pressureService(tankIdx) * psiTokPa; % PSI to kPa
    volumeHeTankMetric = copvVolume(tankIdx) * litersToCubicM;

    % Calculate the mass of the helium in the COPV
    mHelium = (pressureHelium * volumeHeTankMetric) / (heGasConstant * ambTemp);

    currentHeTankPressure = pressureHelium; % [kPa]
    previousHeTankPressure = pressureHelium; % [kPa]

    currentMHe = mHelium;
    previousMHe = mHelium;

    currentSpecificVolumeHe = volumeHeTankMetric / mHelium;
    previousSpecificVolumeHe = volumeHeTankMetric / mHelium;

    currentHeTankTemp = ambTemp;
    previousHeTankTemp = ambTemp;

    currentHeTankTemp = ambTemp;
    time = 0; % Current time of the simulation [s]

    tankHeTempFu = [];
    tankHeTempOx = [];

    currentTankHeTempFu = 0;
    currentTankHeTempOx = 0;

    % Iterates through burn time
    while currentHeTankPressure > pressureFuTank && time < burnTime

        % Caluclates temp of Helium out of the tank
        HeTempOutToFu = (currentHeTankTemp) / ((currentHeTankPressure / pressureFuTank) ^ isentropicExponent);
        HeTempOutToOx = (currentHeTankTemp) / ((currentHeTankPressure / pressureOxTank) ^ isentropicExponent);

        % Calculates mass flow rate out of helium tank
        mDotHetoFu = (pressureFuTank * volumeDotFu) / (heGasConstant * HeTempOutToFu); % [kg/s]
        mDotHetoOx = (pressureOxTank * volumeDotOx) / (heGasConstant * HeTempOutToOx); % [kg/s]
        mDotHe = mDotHetoFu + mDotHetoOx;

        % Calculates new mass inside of helium tank
        exitHeMass = mDotHe * deltaT; % [kg]
        currentMHe = previousMHe - exitHeMass; % [kg]
        currentSpecificVolumeHe = volumeHeTankMetric / currentMHe; % [m^3/kg]

        % Calculates Fuel and Oxidizer Temps
        %         tankHeTempFu(end + 1) = HeTempOutToOx;
        %         tankHeTempOx(end + 1) = HeTempOutToFu;
        %
        %         currentTankHeTempFu = mean(tankHeTempFu);
        %         currentTankHeTempOx = mean(tankHeTempOx);

        % Isentropic Expansion equations
        currentHeTankTemp = previousHeTankTemp * (previousSpecificVolumeHe / currentSpecificVolumeHe) ^ (heSpecificHeatRatio - 1); % [K]
        currentHeTankPressure = previousHeTankPressure * (previousSpecificVolumeHe / currentSpecificVolumeHe) ^ heSpecificHeatRatio; % [kPa]

        % Sets new pressures to previous and such
        previousSpecificVolumeHe = currentSpecificVolumeHe;
        previousHeTankPressure = currentHeTankPressure;
        previousHeTankTemp = currentHeTankTemp;
        previousMHe = currentMHe;
        time = time + deltaT;
    end

    if (currentHeTankPressure > (2 * pressureFuTank)) && (currentHeTankPressure > (2 * pressureOxTank)) && (copvVolume(tankIdx) < bestCopvVolume)
        bestCopvIdx = tankIdx;
        bestCopvVolume = copvVolume(tankIdx);
    end
end

mHelium = mHelium * kilogramToPoundMass; % Converts COPV mass from kg to lbm
end