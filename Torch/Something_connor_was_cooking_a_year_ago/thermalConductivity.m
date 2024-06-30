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
function k = thermalConductivity(temp, aveMolarMass, density, aveDiameter, heatCapacity)
temp = unitConvert(temp,symunit().K);
aveMolarMass = unitConvert(aveMolarMass, symunit().kg / (symunit().mol));
density = unitConvert(density, symunit().kg/(symunit().m ^ 3));
aveDiameter = unitConvert(aveDiameter, symunit().m);
heatCapacity = unitConvert(heatCapacity, symunit().J / (symunit().k * symunit().mol));

%% Code section title/overview
meanSpeed = meanParticleSpeed(temp, aveMolarMass); % The average speed of the gas molecules
meanFreePath = averageFreePath(temp, pressure, aveDiameter); %How far, on average, the molecules can travel without colliding with each other

k = 1 / 3 * heatCapacity * density * meanSpeed * meanFreePath;

end

%% Function definition
function aveSpeed = meanParticleSpeed(temp, aveMolarMass)
    R = 8.3144 * symunit().J / (symunit().K * symunit().mol); % universal Gas constant
    aveSpeed = sqrt((8 * R * temp) / (pi * aveMolarMass));
end

%% Function definition
function mfp = averageFreePath(temp, pressure, averageDiameter)
    kB = 1.380649e-23 *symunit().J / symunit().K; %Boltzmaan constant
    mfp = kB * temp / (sqrt(2) * pi() * (averageDiameter ^ 2) * pressure);
end
%{
    Sources:
    [1]: https://www.tec-science.com/thermodynamics/heat/thermal-conductivity-of-gases/
    [2]: https://en.wikipedia.org/wiki/Mean_free_path
    
%}
