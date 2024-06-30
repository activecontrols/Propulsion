%% Purdue Space Program - Liquids
% Subteam:         Propulsion
% Script:          thermalConductivity
% Description:     Script/function description
% Contributor(s):  Connor McGuire
% Input(s):
%          Variable:    Unit:    Info:
%          var1         -        relevant info for var1
% Output(s):
%          Variable:    Unit:    Info:
%          var4         -        relevant info for var4
% Example call:
%          [var4,var5,var6] = formatGuidelines(var1,var2,var3);

%% Function definition
function k = thermalConductivityFunc(temp, aveMolarMass, density, molRatio, heatCapacity, pressure) % molRatio is mols ethanol:total
aveDiameter = aveDiameterFunc(molRatio, true);
temp = unitConvert(temp,symunit().K);
aveMolarMass = unitConvert(aveMolarMass, symunit().kg / (symunit().mol));
avagadroNum = 6.02214e23 / symunit().mol;

averageDiameter = unitConvert(aveDiameter, symunit().m);
heatCapacity = unitConvert(heatCapacity, symunit().J / (symunit().K * symunit().mol));

%% Code section title/overview
meanSpeed = meanParticleSpeed(temp, aveMolarMass); % The average speed of the gas molecules
meanFreePath = averageFreePath(temp, pressure, averageDiameter); %How far, on average, the molecules can travel without colliding with each other

particleDensity = density * avagadroNum / aveMolarMass;
molarHeatCapacity = heatCapacity * aveMolarMass;
k = particleDensity * meanSpeed * meanFreePath * molarHeatCapacity / (3 * avagadroNum); %[6]
k = unitConvert(k, symunit().J / (symunit().m * symunit().s * symunit().K));
end

%% Function definition
function aveSpeed = meanParticleSpeed(temp, aveMolarMass)
    R = 8.3144 * symunit().J / (symunit().K * symunit().mol); % universal Gas constant
    aveSpeed = sqrt((8 * R * temp) / (pi * aveMolarMass));
end

%% Function definition
function diameter = aveDiameterFunc(molRatio, isEthanolOxygen)
    o2KenDiameter = 346e-12 * symunit().m;
    ethKenDiameter = 4.3e-10 * symunit().m;
    if isEthanolOxygen
        diameter = (ethKenDiameter * molRatio + o2KenDiameter * (1 - molRatio))/(molRatio + 1);
    else
        fprintf("\n\n Non- Ethanol Oxygen not implemented in thermalConductivityFunc");
    end
end

%{
    Sources:
    [1]: https://www.tec-science.com/thermodynamics/heat/thermal-conductivity-of-gases/
    [2]: https://en.wikipedia.org/wiki/Mean_free_path
    [3]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6838788/ //Ethanol Kinetic Diameter
    [4]: https://en.wikipedia.org/wiki/Kinetic_diameter#cite_note-Ismail14-2 //O2 Kinetic Diameter
    [5]: http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/menfre.html
    [6]: https://ocw.snu.ac.kr/sites/default/files/NOTE/543.pdf //page 22 Thermal Conductivity of a gas
%}
