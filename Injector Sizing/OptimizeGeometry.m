%Variables
clc
clear

lOpen = 0.0022045; %(m)
rhoOx = 1141; %(kg/m^3)
Ufuel = 13.1064; %(m/s)
Uox = 11.2776; %(m/s)
surfaceTensionFuel = 0.023; %(N/m)

radialDensity = 786; % kg/m^3
annularDensity = 1141; % kg/m^3
radialVelocity = 13.104; % m/s
annularVelocity = 11.27; % m/s

pintleDiameter = 0.63*0.0254; % m

numberOfSlots = 12;
slotWidth = 0.056568*0.0254; % m
slotHeight = 0.025249*0.0254:0.01*0.0254:(0.09671573+.025249)*0.0254; % m

radialArea = slotWidth * slotHeight * numberOfSlots; % m^2

throttleLevel = 0.4:0.065:1; % (%)
massFlowRateRadial = 0.56064017*throttleLevel; % kg/s
massFlowRateAnnular = 0.6731310*throttleLevel; % kg/s

radialVelocity = massFlowRateRadial / (radialArea * radialDensity); % m/s


blockageFactor = numberOfSlots * slotWidth / (pintleDiameter * pi);




optimizeMomentumRatios(radialDensity, annularDensity, radialVelocity, massFlowRateAnnular, massFlowRateRadial, throttleLevel);




