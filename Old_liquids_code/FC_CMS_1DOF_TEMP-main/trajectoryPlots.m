function trajectoryPlots (varargin)
%{
Purdue Space Program - Liquids
CMS 1DoF - Trajectory Temporal Plots
Talal Zaim

Input variables:  varargin
        If no inputs are given, then the plots for the chosen rocket are
        produced. Plots for a different rocket must require a fileName to
        be inputted.

        fileName (name of file containing the trajectory array)

Output variables: n/a

%}

% Import data
if isempty(varargin)
    fileName = "PickedRocket.mat";
    fprintf("\nCreating plots for the chosen CMS rocket.\n");
end

load(strcat('Trajectory Models\', fileName));

% Separate into variables
time = trajectoryArray(:, 1); % time [s]
alt = trajectoryArray(:, 2); % altitude ASL [ft]
vel = trajectoryArray(:, 3); % inertial velocity [ft/s]
acc = trajectoryArray(:, 4); % inertial acceleration [ft/s^2]
mach = trajectoryArray(:, 5); % inertial mach number
mass = trajectoryArray(:, 6); % total rocket mass [lbm]
atmosPressure = trajectoryArray(:, 7); % atmospheric pressure [psi]
thrust = trajectoryArray(:, 8); % total thrust [lb-f]
pressThrust = trajectoryArray(:, 9); % pressure thrust [lb-f]
%atmosDensity = trajectoryArray(:, 10); % atmospheric density [kg/m^3]

% Find burn index
MECOIndex = find(thrust == 0, 1);

% Create matrix to hold time & altitude values, write to CSV
altDataMatrix = [time(1:end), alt(1:end)];
writematrix(altDataMatrix,'altData.csv');

% Create matrix to hold time & altitude values, write to CSV
velDataMatrix = [time(1:end), vel(1:end)];
writematrix(velDataMatrix,'velData.csv');

% Create matrix to hold time & atmospheric pressure values, write to CSV
atmosPressureDataMatrix = [time(1:end), atmosPressure(1:end)];
writematrix(atmosPressureDataMatrix,'atmosPressureData.csv');

% Create data plots
grapherPSP(1, "Altitude", time(1:MECOIndex), alt(1:MECOIndex), "Time [s]", "Altitude ASL [ft]",...
    'plotX2', time(MECOIndex:end), 'plotY2', alt(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(2, "Velocity", time(1:MECOIndex), vel(1:MECOIndex), "Time [s]", "Inertial Velocity [ft/s]",...
    'plotX2', time(MECOIndex:end), 'plotY2', vel(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(3, "Acceleration", time(1:MECOIndex), acc(1:MECOIndex), "Time [s]", "Acceleration [ft/s^2]",...
    'plotX2', time(MECOIndex:end), 'plotY2', acc(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(4, "Mach", time(1:MECOIndex), mach(1:MECOIndex), "Time [s]", "Mach",...
    'plotX2', time(MECOIndex:end), 'plotY2', mach(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(5, "Mass", time(1:MECOIndex), mass(1:MECOIndex), "Time [s]", "Mass [lbm]",...
    'plotX2', time(MECOIndex:end), 'plotY2', mass(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(6, "Ambient Pressure", time(1:MECOIndex), atmosPressure(1:MECOIndex), "Time [s]", "Atmospheric Pressure [psi]",...
    'plotX2', time(MECOIndex:end), 'plotY2', atmosPressure(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(7, "Total Thrust", time(1:MECOIndex), thrust(1:MECOIndex), "Time [s]", "Thrust [lb-f]",...
    'plotX2', time(MECOIndex:end), 'plotY2', thrust(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);

grapherPSP(8, "Total Thrust", time(1:MECOIndex), pressThrust(1:MECOIndex), "Time [s]", "Pressure Thrust [lb-f]",...
    'plotX2', time(MECOIndex:end), 'plotY2', pressThrust(MECOIndex:end), 'legendNames', ["Burn", "Coast"]);
end

%use export_fig?