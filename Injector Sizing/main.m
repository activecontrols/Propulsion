clc, clear

%% VARIABLES
lOpen = 0.0022045; %(m)
densityOx = 1141; %(kg/m^3)
Ufuel = 13.1064; %(m/s)
Uox = 11.2776; %(m/s)
surfaceTensionFuel = 0.023; %(N/m)
densityFuel = 786; %(kg/m^3)

pintleTipAngles = 0:5:85; % Range of pintle tip angles tested (0, 5, 10, ..., 90)
pintleDiam = 0.0127;%(m)
annularThickness = 0.001; %0.0006858;%(m)
sleeveThickness = 0.000508;%(m)


%% SMD calculation for best atomization
D32 = SMDdesignVar(lOpen, densityOx, Ufuel, Uox, surfaceTensionFuel, pintleTipAngles);

%% Spray Angle calculations
pintleExitArea = lOpen * pi * (pintleDiam - (2 * sleeveThickness));
annularExitArea = pi * ((annularThickness + (pintleDiam / 2))^2 - (pintleDiam / 2)^2);

% mass flow rate at exit
mDotFuel = massFlowRt(pintleExitArea, Ufuel, densityFuel);
mDotOx = massFlowRt(annularExitArea, Uox, densityOx);
%final calc
sprayAngles = sprayAngleDesignVar(pintleTipAngles, mDotFuel, Ufuel, mDotOx, Uox);
    
%% Output
subplot(2, 1, 1);
plot(pintleTipAngles, D32, 'o')
grid on
hold on
xlabel("Pintle Tip Angle (degrees)")
ylabel("Sauter Diameter(m)")
title('Sauter Diameter(m) vs. Pintle Tip Angle (deg)')
hold off

subplot(2, 1, 2);
plot(pintleTipAngles, sprayAngles, 'o')
grid on
hold on
xlabel("Pintle Tip Angle (degrees)")
ylabel("Spray Angle(degrees)")
title('Spray Angle (deg) vs. Pintle Tip Angle (deg)')
hold off

%%fprintf(pintleTipAngle);