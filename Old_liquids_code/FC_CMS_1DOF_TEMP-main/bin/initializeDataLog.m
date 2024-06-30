function [fileName] = initializeDataLog(prefix, writeDate, runDescript)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Data Write Script
Nick Mondora

Input variables:  prefix (RAW or GUD for raw and filtered data)
                  writeDate (date and time of file creation)
                  runDescript (description of the run)

Output variables: fileName (file name)
%}

%% CREATE VARIABLE NAMES
% these must be in the same order they are passed to dataWrite.m in
% frootLoops.m, C-Star is the only metric unit
varNames = ["Outer Diameter (in)", "Dry mass (lbm)", "Wet Mass (lbm)", "Mass of Propellant Tanks (lbm)", "Burn time (s)", ...
    "Expected Isp (s)", "Max Alt (ft)", "Rail Speed (ft/s)", "Max Mach", "Max Accel (g)", ...
    "Chamber Pressure (psi)", "O/F", "Total Mass Flow (lbm/s)", "Thrust (lbf)", ...
    "Height of rocket (ft)", "Aspect Ratio", "Fuel Tank Height (ft)", "Ox Tank Height (ft)", ...
    "Burnout Altitude (ft)", "Fuel Tank ID (in)", "Ox Tank ID (in)", "Expected C-Star (m/s)", ...
    "Total Throat Radius Ablation (in)", "Nozzle Exit Radius (in)", "Thrust to Weight Ratio", ...
    "Fuel Tank Volume (in^3)", "Ox Tank Volume (in^3)", "Helium Tank Volume (L)", "Ambient Pressure at Half Burn (psi)",...
    "Fuel Tank Pressure (psi)", "Ox Tank Pressure (psi)", ...
    "Ox Run Line Velocity (ft/s)", "Fuel Run Line Velocity (ft/s)", "Ox Run Line Size", "Fuel Run Line Size", "Exit Pressure (psi)",...
    "Combustion Temperature (K)", "Altitude At t=10s (ft)", "COPV Mass (lbm)", "Throat Radius (in)", "Vehicle Return Index","FFC (%)","FFC Diameter (in)","Num FFC Holes"];

%% WRITE TO THE CSV
fileName = string(prefix) + string(writeDate) + '.csv';

% insert run description
input = ["Run Description:" string(runDescript)];
input = table(input);
writetable(input, fileName, 'WriteVariableNames', 0, 'WriteMode', 'append');

% write the date of creation
input = ["Run date:" string(writeDate)];
input = table(input);
writetable(input, fileName, 'WriteVariableNames', 0, 'WriteMode','append');

% create the variable names
input = table(varNames);
writetable(input, fileName, 'WriteVariableNames', 0,'WriteMode','append');