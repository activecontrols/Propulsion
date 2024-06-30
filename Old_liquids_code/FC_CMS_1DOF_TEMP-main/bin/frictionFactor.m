%% Friction Factor Calculations

%% Constants
freedomDensityToMetricDensity = 27.68; % lbm/in^3 to g/ml
metricDensityToFreedomDensity = 1 / freedomDensityToMetricDensity; % g/ml to lbm/in^3
inchesToMeters = 1 / 39.37; % in to m
metersToInches = 39.37; % m to in
pascalsToPSI = 1 / 6895; % n/m^2 to psi
psiToPascals = 6895; % psi to Pa
centimetersToMeters = 0.01; % cm to m
gramsToKilograms = 0.001; % g to kg

homeDir = pwd;

%% Parameters
oxTubeOD = 1; % oxygen run line tube outer diameter [in]
oxTubeThick = 0.083; % oxygen run line tube thickness [in]
fuTubeOD = 1; % fuel run line tube outer diameter [in]
fuTubeThick = 0.083; % fuel run line tube thickness [in]  
tubeRoughness = .01; % roughness (drawn tubing = 0.0015, commercial steel = 0.046) [mm]

oxPressure = 600; % oxygen run line pressure [psi]
oxTemperature = 110; % oxygen run line bulk fluid temperature [K]
oxVelocity = 50; % oxygen bulk fluid velocity [ft/s]

fuPressure = 500; % fuel run line pressure [psi]
fuTemperature = 273 + 21; % fuel run line bulk fluid temperature [K]
fuVelocity = 50; % fuel bulk fluid velocity [ft/s]

%% Calculations
oxVelocity = oxVelocity * 12 * inchesToMeters; % [m/s]
fuVelocity = fuVelocity * 12 * inchesToMeters; % [m/s]

oxPressure = oxPressure * psiToPascals / 1000; % [kPa]
fuPressure = fuPressure * psiToPascals / 1000; % [kPa]

oxtubeArea = pi/4 * (oxTubeOD^2 - (2* oxTubeThick)^2); % tube inner area [in^2]
oxDia = oxTubeOD - 2 * oxTubeThick; % tube inner diameter [in]
oxDia = oxDia * inchesToMeters; % [m]
oxRelativeRoughness = tubeRoughness / (oxDia * 1000);

cd('T:\PSP\Liquids\Fluid Systems\REFPROP');
addpath("refprop\");
for index = 1:length(oxPressure)
    oxDensity(index) = refpropm('D', 'T', oxTemperature, 'P', oxPressure(index), 'oxygen'); % [kg/m^3]
    oxViscosity(index) = refpropm('V', 'T', oxTemperature, 'P', oxPressure(index), 'oxygen'); % [Pa * s]
end
cd(homeDir);

oxReynolds = oxDensity .* oxVelocity .* oxDia ./ oxViscosity; 
% fuReynolds = fuDensity * fuVelocity * dia / fuViscosity;oxRelativeRoughness

for index = 1:length(oxReynolds)
    if oxReynolds(index) <= 2300
        oxFrictionFactor(index) = 64 / oxReynolds(index);
    else
        syms r R f
        eq = f^(-1/2) == -2 * log10(r / 3.7 + 2.51 / (R * sqrt(f)));
        eq = subs(eq, r, oxRelativeRoughness);
        eq = subs(eq, R, oxReynolds(index));
        oxFrictionFactor(index) = eval(solve(eq, f));
    end
end

% if fuReynolds <= 2300
%     fuFrictionFactor = 64 / fuReynolds;c
% else
%     r = fuRelativeRoughness;
%     R = fuReynolds;
%     fuFrictionFactor = -2 * (1/(2*wrightOmega((500*R*r)/9287 - log(251/(50*R))) - (1000*R*r)/9287)^2);
% end