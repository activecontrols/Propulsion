%% SIZING

% Simulation Loop
clc
clear;
addpath("lib\");
res = 20;
minTWR = 1;
maxTWR = 2.2;
TWRVec = linspace(minTWR, maxTWR, res);
PropMass = zeros(1,res);
FlightTime = zeros(1,res);
ThrustDev = zeros(1,res);

i = 1;

% Set second parameter to 0 for no flight profile graph generation 
for TWR = TWRVec
    [PropMass(i), FlightTime(i), ThrustDev(i)] = TOAD_3DoF_SIM(TWR, 0);
    i = i + 1;
end


%% Costs
weights = [1/5, 1/5, 3/5];
weights = weights / norm(weights);
cost = weights(1)*PropMass / max(PropMass) +...
       weights(2)*FlightTime / max(FlightTime) +...
       weights(3)*(ThrustDev / max(ThrustDev)).^2;

%% Plots
figure(3);
subplot(1,2,1);
plot(TWRVec, cost, 'g', 'LineWidth',1);
grid on
xlabel('TWR');
ylabel('Cost');
title('Weighted Cost Function');
xlim([minTWR maxTWR])

subplot(1,2,2);
hold on; grid on;
yyaxis left;
plot(TWRVec, PropMass, 'r', 'LineWidth',1);
plot(TWRVec, FlightTime, 'b-', 'LineWidth', 1);
ylabel('Prop Mass [kg], Flight Time [sec]');

yyaxis right;
plot(TWRVec, ThrustDev, 'm', 'LineWidth', 1);
xlabel('TWR');
ylabel('Thrust Deviation from 70% [N^2]');
title('Outputs');
xlim([minTWR maxTWR])
legend('Prop Mass', 'Flight Time','Cumulative Thrust Deviation from 70%',...
    'Location','southwest')

sgtitle('TOAD TWR Sizing Code');