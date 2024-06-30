%% Purdue Space Program - Liquids
% Subteam:         Propulsion
% Script:          thermalConductivity
% Description:     Script/function description
% Contributor(s):  Connor McGuire
% Input(s):
%          Variable:    Unit:    Info:
%          var1         -        relevant info for var1
%          var2         -        relevant info for var2
%          var3         -        relevant info for var3
% Output(s):
%          Variable:    Unit:    Info:
%          var4         -        relevant info for var4
%          var5         -        relevant info for var5
%          var6         -        relevant info for var6
% Example call:
%          [var4,var5,var6] = formatGuidelines(var1,var2,var3);

%% Function definition
function simpleIgnitionEnergy(scalePercent)
    fprintf("\nStarting Simple Ignition Energy\n")
    
    fdot = 4.94 * symunit().lbm / symunit().s; % Fuel mass flow rate in lbm/s
    fdot = unitConvert(fdot, symunit().kg / symunit().s); %Convert to kg/s
    
    odot = 5.93 * symunit().lbm / symunit().s; % Oxidizer mass flow rate in lbm/s 
    odot = unitConvert(odot, symunit().kg / symunit().s); %Convert to kg/s
    
    fdot = fdot * scalePercent;
    odot = odot * scalePercent;
    
    tempAdiabatic = 3730 * symunit().K; % The Adiabatic flame temperature of ethanol/Oxygen [4]
    tempAdiabatic = unitConvert(tempAdiabatic, symunit().K);
    
    tempIgnition = 365 * symunit().K;
    tempIgnition = unitConvert(tempIgnition, symunit().K);
    
    tempStarting = ["LOX", 100 * symunit().K; % Saying its starting at 100, for simpler calcs cause the schomate only goes down that far
                    "Ethanol", 293 * symunit().K];
    
    boilingStats = [ "LOX" , 90.19 * symunit().K, 213e3 * symunit().J / symunit().kg;
                   "Ethanol" , 351.5 * symunit().K, 919e3 * symunit().J / symunit().kg]; % [ NAME, Temp in K , enthalpy of evaporation (J/kg) ] 
    
    o2BoilingPower = odot * boilingStats(1,3);
    ethanolBoilingPower = fdot * boilingStats(2,3);
    
    deltaH_c_ox = SchomateOxygen([double(separateUnits(unitConvert(tempStarting(1,2), symunit().K))), double(separateUnits(tempIgnition))], "deltaH") * symunit().J / symunit().mol;
    deltaH_c_ethanol = SchomateEthanol([double(separateUnits(unitConvert(tempStarting(2,2), symunit().K))), double(separateUnits(tempIgnition))], "deltaH") * symunit().J / symunit().mol;
    
    o2MolarMass = 15.999 * symunit().g / symunit().mol;
    ethanolMolarMass = 46.07 * symunit().g / symunit().mol;
    
    deltaH_c_ethanol = deltaH_c_ethanol / ethanolMolarMass;
    deltaH_c_ox = deltaH_c_ox / o2MolarMass;
    
    o2HeatPower = odot * deltaH_c_ox;
    ethanolHeatPower = fdot * deltaH_c_ethanol;
    
    totalPower = o2HeatPower + o2BoilingPower + ethanolHeatPower + ethanolBoilingPower;
    fprintf("Total Power: ")
    apx(totalPower)
    
    fprintf("Finished Simple Ignition Energy \n")
end