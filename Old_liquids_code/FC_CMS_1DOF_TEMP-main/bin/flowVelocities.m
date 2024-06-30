%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Flow Velocities
Talal Zaim

Calculates flow velocities in different tube thicknesses
%}

%% Definitions
bzbLMDot = 3.0385; % BZB Lox mass flow rate, lb/s, taken from HF2
bzbFMDot = 1.1141; % BZB fuel mass flow rate, lb/s, taken from HF2
RocketThreeOF = 2.2;

LoxDensityBZB = 71.2303 / 12^3; %lb/in^3
LoxDensityR3 = 0.04364; % lb/in^3
MethDensity = 24.76 / 12^3; %lb/in^3
JetADensity = 0.0290463428; % lb/in^3

pipeArea.dash4 = pi / 4 * (0.25^2 - 0.035^2); % in^2

pipeArea.dash8_049 = pi / 4 * (0.5^2 - 0.049^2); % in^2
pipeArea.dash8_065 = pi / 4 * (0.5^2 - 0.065^2); % in^2
pipeArea.dash8_083 = pi / 4 * (0.5^2 - 0.083^2); % in^2

pipeArea.dash12_065 = pi / 4 * (0.75^2 - 0.065^2); % in^2
pipeArea.dash12_083 = pi / 4 * (0.75^2 - 0.083^2); % in^2
pipeArea.dash12_095 = pi / 4 * (0.75^2 - 0.095^2); % in^2
pipeArea.dash12_120 = pi / 4 * (0.75^2 - 0.120^2); % in^2

pipeArea.dash16_065 = pi / 4 * (1^2 - 0.065^2); % in^2
pipeArea.dash16_083 = pi / 4 * (1^2 - 0.083^2); % in^2
pipeArea.dash16_095 = pi / 4 * (1^2 - 0.095^2); % in^2
pipeArea.dash16_120 = pi / 4 * (1^2 - 0.120^2); % in^2
pipeArea.dash16_188 = pi / 4 * (1^2 - 0.188^2); % in^2

kFactor90Deg = 0.18;

names = fieldnames(pipeArea);

for i = 1: size(names,1)
    area = pipeArea.(names{i});
    BZBVelo(i, :) = [bzbLMDot / (LoxDensityBZB * area), bzbFMDot / (MethDensity * area)];
end

BZBVelo = BZBVelo ./ 12;

mdot = 1:1:20;
for index = 1 : length(mdot)
    mdotFu = mdot(index) / (RocketThreeOF + 1);
    mdotLox = mdot(index) - mdotFu;

    for i = 1: size(names,1)
        area = pipeArea.(names{i});
        RocketThreeVelo(i, :, index) = [mdotLox / (LoxDensityR3 * area), mdotFu / (JetADensity * area)];
        
    end
end

RocketThreeVelo = RocketThreeVelo ./ 12;
numMDots = size(RocketThreeVelo, 3);

figure();
plot(mdot, reshape(RocketThreeVelo(2, 1, :), [numMDots, 1, 1]), 'bs');
hold on;
plot(mdot, reshape(RocketThreeVelo(5, 1, :), [numMDots, 1, 1]), 'b^');
plot(mdot, reshape(RocketThreeVelo(9, 1, :), [numMDots, 1, 1]), 'bd');

plot(mdot, reshape(RocketThreeVelo(2, 2, :), [numMDots, 1, 1]), 'rs');
plot(mdot, reshape(RocketThreeVelo(5, 2, :), [numMDots, 1, 1]), 'r^');
plot(mdot, reshape(RocketThreeVelo(9, 2, :), [numMDots, 1, 1]), 'rd');

legend (["LOx -8", "LOx -12", "LOx -16", "Fu -8", "Fu -12", "Fu -16"]);
xlabel ("Total mdot [lb/s]");
ylabel ("Run Line Flow Velocity [ft/s]");
grid on;
title ("Run Lines Flow Velocities (OF = 2.5)");