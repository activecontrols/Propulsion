function[mStructures, rocketHeight] = structures(tubeOD, tankHeight, heightHe)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Structures Function
Alex Porter, Jonah Fouts, Julian Petrillo, Matthew Lewton

Input variables:  tubeOD [in]
                  tankHeight [in]
                  heightHe [in]

Output variables: mStructuresMass [lbm]
                  rocketHeight [in]
%}

% Givens
% bzbOD = 6.5; % [in] not using bzb estimates anymore
% mBZBNosecone = 2.28; % [lbm] not using bzb estimates anymore
engineLength = 24; % [in]
fincanLength = 18; % [in]
midPlumLength = 20; % [in]
% recoveryLength = 9; % [in] no recovery tube
avionicsLength = 6; % [in]
noseconeLength = tubeOD * 5; % [in]
mTip = 1.2; % lbm of nosecone tip 
couplerEstimate = 5; % lbm for aluminum tube couplers

rocketHeight = tankHeight + engineLength + fincanLength + midPlumLength + avionicsLength + noseconeLength + heightHe; % [in]

% Tube Size perameters
% odRatioBZB = tubeOD / bzbOD;

% Layup Properties
LayerCt = 9; % 4 on each side of the sandwich, plus 1 layer to approx sandwich core.
layerThickness = .025/3; % [in], based on previous vaccum bagged tubes
cfDensity = .052421; % [lb/in^3]

% Nosecone Calculations
noseconeSA = pi * (tubeOD / 2) * (tubeOD * sqrt((tubeOD / 2) ^ 2 + noseconeLength ^ 2)); % Surface area based on cone estimate
mNosecone = (noseconeSA * layerThickness * 7 * cfDensity) + mTip; % [lbm] 

% Helium Mid Calculations
mMidHelium = heightHe * pi * tubeOD * LayerCt * layerThickness * cfDensity + (2 * couplerEstimate); % mass of tube and couplers around He tank

% Mid Plumbing mass calcs
mMidPlumbing = (tubeOD * pi * midPlumLength * 3 * layerThickness * cfDensity) + 6; % assumes 3 layer panels and 3 1x1x20 inch struts 

%fincan mass calcs
mFincan = (tubeOD * pi * fincanLength * 3 * layerThickness * cfDensity) + 9 + 6; % assumes 3 layer panels and 3 1.5x1x20 inch struts (6 lb total fin mass out of matt's ass)

% Total Mass of Airframe
mRecovery = 15; % [lbm] (Colin's estimate on 8/28/28)
mStructures = mNosecone + mMidHelium + mRecovery + mMidPlumbing + mFincan; % [lbm]
end