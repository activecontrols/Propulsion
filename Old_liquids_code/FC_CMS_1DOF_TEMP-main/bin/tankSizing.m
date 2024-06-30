function [heightOx, heightFu, mTanks, mFluidSys, mFluids, volumeFuTank, volumeOxTank, ...
    tubeIDFu, tubeIDOx, pressureFuTank, pressureOxTank, velOx, velFu, runLineOx, runLineFu, bestCopvIdx] = ...
    tankSizing(burnTime, mDotFu, OF, pressureChamber, isRegen, availThicknesses, tubeOD, pressureService, copvVolume, mCopv, lowerPlumbingDrop)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Tank Sizing
Jonah Fouts, Matt Currie, Cameron Williams, Tyler Trostle, Erik Anderson, Jose Lara

Input variables:  burnTime [s] (Burn Time)
                  mDotFu [lb/s] (mass flow fuel)
                  OF (O/F Ratio)
                  pressureChamber [psi] (chamber pressure)
                  isRegen [boolean] (is it a regen engine)
                  availThicknesses [in] (current tank thickness)
                  tubeOD [in] (tube outer diameter)
                  lowerPlumbingDrop [% or -1] (pressure drop setting)

Output variables: heightOx [in] (height of LOX tank)
                  heightFu [in] (height of Fuel tank)
                  heightHe [in] (height og Helium tank)
                  mTanks [lbm] (mass of propellant tanks)
                  mFluidSys [lbm] (mass of fluid system)
                  mFluids [lbm] (mass of the fluids)
                  volumeFuTank [in^3] (volume of the fuel tank)
                  volumeOxTank [in^3] (volume of the ox tank)
                  volumeHeTank [in^3] (volume of the he tank)
                  tubeIDFu [in] (fuel tank inner diameter)
                  tubeIDOx [in] (ox tank inner diamter)
                  pressureFuTank [psi] (fuel tank pressure)
                  pressureOxTank [psi] (ox tank pressure) 
                  velOx [ft/s] (ox runline velocity)
                  velFu [ft/s] (fu runline velocity)
                  velOxDashSixteen [ft/s] (ox runline velocity for -16 tube)                
%}

freedomDensityToMetricDensity = 27.68; % lbm/in^3 to g/ml
metricDensityToFreedomDensity = 1 / freedomDensityToMetricDensity; % g/ml to lbm/in^3
inchesToMeters = 1 / 39.37; % in to m
metersToInches = 39.37; % m to in
paToPSI = 1 / 6895; % n/m^2 to psi
centimetersToMeters = 0.01; % cm to m
gramsToKilograms = 0.001; % g to kg

% densities
densLOXMetric = 1.014; % [g/cm^3]
% densJetAMetric = 0.804; % density of JetA [g/cm^3]
% densJetAMetric = 0.85131; % density of Ethanol (25% H2O) [g/cm^3]
densJetAMetric = 0.8; % density of Ethanol (5% H2O) [g/cm^3]
densAluminum = 2.7; % [g/cm^3]
mCOPVperHeVol = 0.429; % [g/cm^3] (mass of COPV per cubic cm of He)

% conversions
densLOX = densLOXMetric * metricDensityToFreedomDensity; % LOX density [lb/in^3]
densJetA = densJetAMetric * metricDensityToFreedomDensity; % RP1 density [lb/in^3]
densAluminum = densAluminum * metricDensityToFreedomDensity; % Aluminum 6061 density [lb/in^3]
mCOPVperHeVol = mCOPVperHeVol * metricDensityToFreedomDensity; % Composite tank density [lb/in^3]
yieldStrengthAL = 40000; % 6061 Aluminum yield strength [psi]

% Flow velocities
tubeAreaDashEight = pi / 4 * (0.5 - (2 * 0.049)) ^ 2; % [in^2]
tubeDiamDashEight = 0.5 - (2 * 0.049); %[in]

tubeAreaDashTwelve = pi / 4 * (0.75 - (2 * 0.065)) ^ 2; % [in^2]
tubeDiamDashTwelve = 0.75 - (2 * 0.065); % [in]

tubeAreaDashSixteen = pi / 4 * (1 - (2 * 0.065)) ^ 2; % [in^2]
tubeDiamDashSixteen = 1 - (2 * 0.065); % [in]

FrictionFactorOx = 0.0125; % Colebrook w/ absolute roughness 1.5e-3 [mm], assumes 600 [psi], 110 [K] in runlines 
FrictionFactorFu = 0.017; % assumes 500 [psi], 294 [K] in runlines

mDotOx = mDotFu * OF; % [lbm/s]

velOxDashEight = mDotOx / (densLOX * tubeAreaDashEight) / 12; % [ft/s]
velOxDashTwelve = mDotOx / (densLOX * tubeAreaDashTwelve) / 12; % [ft/s]
velOxDashSixteen = mDotOx / (densLOX * tubeAreaDashSixteen) / 12; % [ft/s]

velFuDashEight = mDotFu / (densJetA * tubeAreaDashEight) / 12; % [ft/s]
velFuDashTwelve = mDotFu / (densJetA * tubeAreaDashTwelve) / 12; % [ft/s]
velFuDashSixteen = mDotFu / (densJetA * tubeAreaDashSixteen) / 12; % [ft/s]

oxVelocityLimit = 40; %[ft/s] (Zucrow fluid velocity limit)
fuVelocityLimit = 50; 

if velOxDashSixteen <= oxVelocityLimit %[ft/s]
    runLineOx = 8;
    velOx = velOxDashEight;
    diamOx = tubeDiamDashEight;
    if velOxDashEight >= oxVelocityLimit %[ft/s]
        runLineOx = 12;
        velOx = velOxDashTwelve;
        diamOx = tubeDiamDashTwelve;
    end
    if velOxDashTwelve >= oxVelocityLimit %[ft/s]
            runLineOx = 16;
            velOx = velOxDashSixteen;
            diamOx = tubeDiamDashSixteen;
    end
else 
    runLineOx = 0.000001;
    velOx = 0.000001;
    diamOx = 0.000001;
end

% if velFuDashSixteen <= fuVelocityLimit %[ft/s]
%     runLineFu = 16;
%     velFu = velFuDashSixteen;
%     diamFu = tubeDiamDashSixteen;
%     if velFuDashTwelve <= fuVelocityLimit %[ft/s]
%         runLineFu = 12;
%         velFu = velFuDashTwelve;
%         diamFu = tubeDiamDashTwelve;
%         if velFuDashEight <= fuVelocityLimit %[ft/s]
%             runLineFu = 8;
%             velFu = velFuDashEight;
%             diamFu = tubeDiamDashEight;
%         end
%     end
% else 
%     runLineFu = 0.000001;
%     velFu = 0.000001;
%     diamFu = 0.000001;
% end

  runLineFu = 16;
    velFu = velFuDashSixteen;
    diamFu = tubeDiamDashSixteen;

finCanLength = 20 * inchesToMeters;
arbitraryFuTankLength = 53 * inchesToMeters; % an overestimate of the length of the fuel tank [m]

% Fuel Mass
mFuel = burnTime * mDotFu; % [lbm]

% System pressure drops
injectorDrop = 0.2; % 20% pressure drop from the injector
regenDrop = 0.4; % 40% of chamber pressure drop from regen

if lowerPlumbingDrop == -1
    plumbingDropOx = paToPSI  * (FrictionFactorOx * ((finCanLength + arbitraryFuTankLength) / (diamOx * inchesToMeters)) ... 
                    * (1000 * densLOXMetric) * ((velOx * 12 * inchesToMeters) ^ 2) / 2); % Pressure drop over lower ox plumbing [psi]
    plumbingDropFu = paToPSI * (FrictionFactorFu * (finCanLength / (diamFu * inchesToMeters)) ...
                    * (1000 * densJetAMetric) * ((velFu * 12 * inchesToMeters) ^ 2) / 2); % Pressure drop over lower plumbing (fu) [psi]
else
    plumbingDropFu = lowerPlumbingDrop / 100; % Pressure drop over lower plumbing [%]
    plumbingDropOx = lowerPlumbingDrop / 100; % Pressure drop over lower plumbing [%]
end

venturiDrop = 0.20; % 20% pressure drop from venturi

kfactorDropFu = paToPSI * 0.35 * (1000 * densJetAMetric) * ((velFu * 12 * inchesToMeters) ^ 2) / 2; % K-factor pressure drop over 45 degree bend, K = 0.35 [psi]
kfactorDropOx = paToPSI * 0.35 * (1000 * densLOXMetric) * ((velOx * 12 * inchesToMeters) ^ 2) / 2; % K-factor pressure drop over 45 degree bend, K = 0.35 [psi]

% Pressure Head
% pressureHeadFu = paToPSI * (densJetAMetric * 9.81 * finCanLength * gramsToKilograms / (centimetersToMeters)^3); % head pressure from plumbing in fin can in Fu run line [psi]
% pressureHeadOxPlumbing = paToPSI * (densLOXMetric * 9.81 * finCanLength * gramsToKilograms / (centimetersToMeters)^3); % head pressure from plumbing in fin can in Ox run line [psi]
% pressureHeadOxPassThru = paToPSI * (densLOXMetric * 9.81 * arbitraryFuTankLength * gramsToKilograms / (centimetersToMeters)^3); % head pressure from Fu Tank pass-through [psi]

% Fuel Tank Pressures
pressureFuTank = pressureChamber / (1 - injectorDrop); % Calculates injector inlet pressure (20%)
if isRegen
    pressureFuTank = pressureFuTank + (pressureChamber * regenDrop); % Calculates regen inlet pressure (30%)
end
if lowerPlumbingDrop == -1
    pressureFuTank = pressureFuTank + plumbingDropFu + (2 * kfactorDropFu); % Calculates venturi outlet pressure/lower plumbing inlet (friction factor derived)
else
    pressureFuTank = pressureFuTank / (1- plumbingDropFu) + (2 * kfactorDropFu); % Calculates venturi outlet pressure/lower plumbing inlet (given percentage)
end
pressureFuTank = pressureFuTank / (1 - venturiDrop); % Calculates final fuel tank pressure/venturi inlet (20%)

% Ox Tank Pressures
pressureOxTank = pressureChamber / (1 - injectorDrop); % Calculates injector inlet pressure (20%)
if lowerPlumbingDrop == -1
    pressureOxTank = pressureOxTank + plumbingDropOx + (2 * kfactorDropOx); % Calculates venturi outlet pressure/lower plumbing inlet (friction factor derived)
else
    pressureOxTank = pressureOxTank / (1- plumbingDropOx) + (2 * kfactorDropOx); % Calculates venturi outlet pressure/lower plumbing inlet (given percentage)
end
pressureOxTank = pressureOxTank / (1 - venturiDrop); % Calculates venturi inlet pressure

% Tank Volumes
ullageFactor = 1.1; % Ullage volume factor
volumeOxTank = mFuel * OF / densLOX * ullageFactor; % [in^3]
volumeFuTank = mFuel / densJetA * ullageFactor; % [in^3]

% Helium Tank Sizing
[mHelium, bestCopvIdx] = heliumTankSequal(pressureFuTank, pressureOxTank, mDotFu, mDotOx, burnTime, pressureService, copvVolume, mCopv); % Calls Helium Tank function

% Tank Thicknesses
safetyFactor = 1.3; % tank safety factor
qualityFactor = 1; %quality factor for seamless aluminum tube - B31.3
BoardmanFactor = 0.4; %Boardman correction factor for circumferential stress on the ID of a cylinder
weldFactor = 0.80; %weld and heat affected zone knockdown

thickRequirementFu = safetyFactor * pressureFuTank * tubeOD / (2 * (yieldStrengthAL * qualityFactor * weldFactor + pressureFuTank * BoardmanFactor)); % (in)
thickRequirementOx = safetyFactor * pressureOxTank * tubeOD / (2 * (yieldStrengthAL * qualityFactor * weldFactor + pressureOxTank * BoardmanFactor)); % (in)

% Find max fuel tank ID
tankBoolFu = availThicknesses >= thickRequirementFu;
tubeThickFu = availThicknesses(find(tankBoolFu == true, 1)); % [in]
tubeIDFu = tubeOD - (2 * tubeThickFu); % [in]

% Find max ox tank ID
tankBoolOx = availThicknesses >= thickRequirementOx;
tubeThickOx = availThicknesses(find(tankBoolOx == true, 1)); % [in]
tubeIDOx = tubeOD - (2 * tubeThickOx); % [in]
tubeODHe = tubeOD - 0.25; % [in]
tubeIDHe = tubeODHe - 0.35; % [in]

if isempty(tubeIDOx) || isempty(tubeIDFu) || bestCopvIdx == 0
    % Set all returned values to intmax to check if a run needs to be skipped
    heightOx = intmax; % [fucktons]
    heightFu = intmax; % [fucktons]
    heightHe = intmax; % [fucktons]

    mFluidSys = intmax; % [fucktons]
    mFluids = intmax; % [fucktons]
    mTanks = intmax; % [fucktons]
else
    % Tank Heights
    heightOx = (volumeOxTank - (4 / 3) * pi * (tubeIDOx / 2) ^ 3) / (pi * (tubeIDOx / 2) ^ 2) + tubeIDOx; % [in]
    heightFu = (volumeFuTank - (4 / 3) * pi * (tubeIDFu / 2) ^ 3) / (pi * (tubeIDFu / 2) ^ 2) + tubeIDFu; % [in]
    heightHe = 60; %(volumeHeTank - (4 / 3) * pi * (tubeIDHe / 2) ^ 3) / (pi * (tubeIDHe / 2) ^ 2) + tubeIDHe; % [in]
    
    % Tank Cylinder Volumes
    materialLoxTankVol = pi * (heightOx - tubeOD) * ((tubeOD / 2) ^ 2 - (tubeIDOx / 2) ^ 2); % [in^3]
    materialFuTankVol = pi * (heightFu - tubeOD) * ((tubeOD / 2) ^ 2 - (tubeIDFu / 2) ^ 2); % [in^3]
    
    % Tank Cylinder Masses
    mOxCyl = materialLoxTankVol * densAluminum; % [lbm]
    mFuCyl = materialFuTankVol * densAluminum; % [lbm]
    mHeTank = mCopv(bestCopvIdx); %volumeHeTank * mCOPVperHeVol; % [lbm]

    % Prop Bulkhead Masses (Assume split tanks, 4x bulkheads)
    volTankBulkheadFu = (2 / 3) * pi * (tubeOD ^ 3 - tubeIDFu ^ 3) / 8; % [in^3]
    mTankBulkheadFu = 2 * densAluminum * volTankBulkheadFu * safetyFactor; % [lbm]
    volTankBulkheadOx = (2 / 3) * pi * (tubeOD ^ 3 - tubeIDOx ^ 3) / 8; % [in^3]
    mTankBulkheadOx = 2 * densAluminum * volTankBulkheadOx * safetyFactor; % [lbm]
    
    % Total tank masses
    mOxTank = mOxCyl + mTankBulkheadOx; % [lbm]
    mFuTank = mFuCyl + mTankBulkheadFu; % [lbm]
    
    % Plumbing Mass
    mPlumbing = 30; % [lbm]

    % Fluid System Mass
    mTanks = mOxTank + mFuTank; % [lbm]
    mFluidSys = mTanks + mHeTank + mPlumbing; % [lbm]

    % Fluids Mass
    mFluids = mFuel * (1 + OF) + mHelium; % [lbm]
end