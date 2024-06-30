% Written By Connor McGuire to evaluate the required energy to raise the
% propellent to the ignition temperature using the Schomate equation for
% pressure
fprintf("\n")

fdot = 4.94; % Fuel mass flow rate in lbm/s
fdot = unitConvert(fdot * symunit().lbm / symunit().s, symunit().kg / symunit().s); %Convert to kg/s

odot = 5.93; % Oxidizer mass flow rate in lbm/s 
odot = unitConvert(odot * symunit().lbm / symunit().s, symunit().kg / symunit().s); %Convert to kg/s

mdot = fdot + odot; %Total mass flow rate in kg/s
oxMassPercent = odot / mdot;

pressureChamber = 14.7 * symunit().psi;
pressureChamber = unitConvert(pressureChamber, symunit().Pa);

%%Critical Diameter Calculations
% D_c * Laminar combustion velocity / heat diffusivity = Shape Constant
% D_c = Shape Constant * heat diffusivity / laminar combustion Velocity
% heat diffusivity = Thermal conductivity / (volumetric density * heat Capacity)
% D_c = shape Constant * Thermal conductivity / (Laminar Combustion Velocity * volumetric density * heat Capacity)
propigationShape = "Sphere"; %"Plane" or "Sphere" 
if propigationShape == "Plane"
    propigationConstant = 4;
elseif propigationShape == "Sphere"
    propigationConstant = 12;
end    

speedLaminar = 0.45 * symunit().m / symunit().s; % The laminar flame speed of Ethanol. [3] This is for air as an oxidizer, but this gives us a lower bound
o2MolarMass = 32 * symunit().g / symunit().mol;
ethanolMolarMass = 46.07 * symunit().g / symunit().mol;
tempAdiabatic = 3730 * symunit().K; % The Adiabatic flame temperature of ethanol/Oxygen [4]

propMolarMass = oxMassPercent * o2MolarMass + ethanolMolarMass * (1 - oxMassPercent);
propMolRatio = (1 - oxMassPercent) * propMolarMass / ethanolMolarMass; %Ethanol to total
R_constant = 8.3145 * symunit().J / (symunit().mol * symunit().K);
densityGas = propMolarMass * pressureChamber / (R_constant * tempAdiabatic); %density = molar mass * pressure / R * temperature

temporaryValue = double(separateUnits(tempAdiabatic));
C_pOxygen = SchomateOxygen(double(separateUnits(tempAdiabatic)), "c_p") * symunit().J / (symunit().mol * symunit().K); %C_p in J/mol * k
C_pEthanol = SchomateOxygen(double(separateUnits(tempAdiabatic)), "c_p") * symunit().J / (symunit().mol * symunit().K); %C_p in J/mol * k
C_pOxygen = C_pOxygen / o2MolarMass; %Convert the specific heat to Energy / temperature * mass 
C_pEthanol = C_pEthanol / ethanolMolarMass; %Convert the specific heat to Energy / temperature * mass 
C_p = oxMassPercent * C_pOxygen + (1 - oxMassPercent) * C_pEthanol;
thermalConductivity = thermalConductivityFunc(tempAdiabatic, propMolarMass, densityGas, propMolRatio, C_pOxygen, pressureChamber);
heatDiffusivity = thermalConductivity / (densityGas * C_p);

% D_c = shape Constant * Thermal conductivity / (Laminar Combustion Velocity * volumetric density * heat Capacity)
D_c = propigationConstant * heatDiffusivity / speedLaminar;
D_c = unitConvert(D_c, symunit().m); %Convert the critical Diameter to meters
[outputVal, outputUnit] = separateUnits(unitConvert(D_c, symunit().m));
fprintf("Critical Diameter: %e %s", outputVal, outputUnit)

V_c = 4/3 * sym(pi) *(D_c / 2)^3;   % Calculate the critical volume
mass_c = densityGas * V_c;          % The mass of the critical volume
oxMass_c = oxMassPercent * mass_c;  % The mass of the oxygen in the critical volume
ethanolMass_c = (1 - oxMassPercent) * mass_c; % The mass of the ethanol in the critical volume
oxMol_c = oxMass_c / o2MolarMass;             % The number of mols of o2 in the critical volume
ethanolMol_c = ethanolMass_c / ethanolMolarMass; % The number of mols of ethanol in the critical volume

tempStarting = ["LOX", 100 * symunit().K; % Saying its starting at 100, for simpler calcs cause the schomate only goes down that far
                "Ethanol", 293 * symunit().K];

boilingStats = [ "LOX" , 90.19 * symunit().K, 213e3 * symunit().J / symunit().kg;
               "Ethanol" , 351.5 * symunit().K, 919e3 * symunit().J / symunit().kg]; % [ NAME, Temp in K , enthalpy of evaporation (J/kg) ] 

boilingEnergy = oxMass_c * boilingStats(1,3) + ethanolMass_c * boilingStats(2,3);
deltaH_c_ox = SchomateOxygen([double(separateUnits(unitConvert(tempStarting(1,2), symunit().K))), double(separateUnits(unitConvert(tempAdiabatic, symunit().K)))], "deltaH") * symunit().J / symunit().mol;
deltaH_c_ethanol = SchomateEthanol([double(separateUnits(unitConvert(tempStarting(2,2), symunit().K))), double(separateUnits(unitConvert(tempAdiabatic, symunit().K)))], "deltaH") * symunit().J / symunit().mol;
energyCritical = oxMol_c * deltaH_c_ox + ethanolMol_c * deltaH_c_ethanol;
[outputVal, outputUnit] = separateUnits(unitConvert(unitConvert(energyCritical, symunit().m^3 / symunit().pa),symunit().J));
fprintf("\nCritical Energy: %e %s", outputVal, outputUnit)

velInj = 4.2 * symunit().m / symunit().s; % Velocity of the injected fluid into the chamber in m/s (Assumes directly downwards)

lenChamber = 19.4 * symunit().inch; % the length of the combustion chamber in inches
lenChamber = unitConvert(lenChamber, symunit().m);

%tRes is the time that the propellent has to combust in the chamber. For simplicity's sake, we assume that the propellent is in free fall
%tRes = abs(sqrt(((vInj)^2 + 2 * 9.81 * (symunit().m / sy munit().s ^ 2) * lenChamber)) - vInj) / (9.81 * (symunit().m / symunit().s ^ 2));
%tRes = unitConvert(tRes, symunit().s);
%tres = tRes * symunit().s;

fprintf("\n")
%%Sources
%{
[1] Chapter 11 Propellant Ignition and Flame Propagation (AIAA) https://arc.aiaa.org/doi/pdf/10.2514/5.9781600866760.0405.0435
[2] Ethanol Equations (Thermal Fluids): http://www.thermalfluidscentral.org/encyclopedia/index.php/Thermophysical_Properties:_Ethanol
[3] Laminar Burning Velocity of ethanol : Gülder OL. Laminar burning velocities of methanol, ethanol and iso-octane-air mixtures. In: Nineteenth Symposium on combustion. Pittsburgh, PA: The Combustion Institute; 1982. p. 275–81.
[4] Adiabatic Flame Temperatures (Engineering Flame Temperatures) : https://www.engineeringtoolbox.com/adiabatic-flame-temperature-d_996.html
[5] Oxygen Schomate Variables (NIST) : https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1&Type=JANAFG&Plot=on
%}