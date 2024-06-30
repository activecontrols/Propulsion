%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Data Plotting
Jonah Fouts

Run after data_RUNME.m before clearing workspace
%}

% limits
altLimit = 150000;
altMin = 75000;

limit_boolean = (maxAlt < altLimit) & (maxAlt > altMin);

% created limited data
limit_AR = AR(limit_boolean);
limit_burnoutAlt = burnoutAlt(limit_boolean);
limit_burnTime = burnTime(limit_boolean);
limit_combustionTemperature = combustionTemperature(limit_boolean);
limit_exitPresHalfBurn = exitPresHalfBurn(limit_boolean);
limit_exitRadius = exitRadius(limit_boolean);
limit_expectedCStar = expectedCStar(limit_boolean);
limit_expectedIsp = expectedIsp(limit_boolean);
limit_heightFu = heightFu(limit_boolean);
limit_heightLOX = heightLOX(limit_boolean);
limit_heTankVolume = heTankVolume(limit_boolean);
limit_loxTankVolume = loxTankVolume(limit_boolean);
limit_maxAccel = maxAccel(limit_boolean);
limit_maxAlt = maxAlt(limit_boolean);
limit_maxMach = maxMach(limit_boolean);
limit_mDot = mDot(limit_boolean);
limit_mDry = mDry(limit_boolean);
limit_mTanks = mTanks(limit_boolean);
limit_mWet = mWet(limit_boolean);
limit_OF = OF(limit_boolean);
limit_outerDiameter = outerDiameter(limit_boolean);
limit_pressureChamber = pressureChamber(limit_boolean);
limit_pressureExit = pressureExit(limit_boolean);
limit_pressureFuelTank = pressureFuelTank(limit_boolean);
limit_pressureOxTank = pressureOxTank(limit_boolean);
limit_railSpeed = railSpeed(limit_boolean);
limit_rocketHeight = rocketHeight(limit_boolean);
rocketNum = 1:length(limit_maxAlt);
limit_rpTankVolume = rpTankVolume(limit_boolean);
limit_throatAblation = throatAblation(limit_boolean);
limit_thrust = thrust(limit_boolean);
limit_tubeIDFu = tubeIDFu(limit_boolean);
limit_tubeIDOx = tubeIDOx(limit_boolean);
limit_twr = twr(limit_boolean);
limit_velOx = velOx(limit_boolean);
limit_velFu = velFu(limit_boolean);
limit_runLineOx = runLineOx(limit_boolean);
limit_runLineFu = runLineFu(limit_boolean);
limit_radiusThroatInit = radiusThroatInit(limit_boolean);

% find max rocket
max_maxAlt = max(limit_maxAlt);
index = find(maxAlt == max_maxAlt);
max_burnTime = burnTime(index);
max_expectedIsp = expectedIsp(index);
max_mDry = mDry(index);
max_rocketHeight = rocketHeight(index);
max_mDot = mDot(index);
max_thrust = thrust(index);
max_pressureChamber = pressureChamber(index);
max_AR = AR(index);

% print max rocket
fprintf("Rocket that maximizes altitude: \n")
fprintf("Altitude: %.2f ft\n",max_maxAlt)
fprintf("Burn Time: %.2f s\n",max_burnTime)
fprintf("Isp: %.2f s\n",max_expectedIsp)
fprintf("Dry Mass: %.2f lbs\n",max_mDry)
fprintf("Rocket Height: %.2f ft\n",max_rocketHeight)
fprintf("M dot: %.2f lbs/s\n",max_mDot)
fprintf("Thrust: %.2f lbf\n",max_thrust)
fprintf("Chamber Pressure: %.2f psi\n",max_pressureChamber)
fprintf("Aspect Ratio: %.2f\n",max_AR)

f = figure(1);
h = plot(rocketNum, sort(limit_maxAlt), 'b.');
grid on
hold on
title("Maximum Altitude vs. Iteration")
xlabel("Rocket Number")
ylabel("Maximum Altitude (ft)")
PSPStyler (f, h, "Light");
saveas(f, "fig1.png");

%% PLOTS
f = figure(2);
h = plot(limit_pressureChamber,limit_maxAlt,'b.');
grid on
title("Maximum Altitude vs. Chamber Pressure")
xlabel("Chamber Pressure (psi)")
ylabel("Maximum Altitude")
PSPStyler (f, h, "Light");

f = figure(3);
h = plot(limit_OF, limit_maxAlt,'b.');
grid on
title('Max Altitude vs. OF')
xlabel('OF Ratio')
ylabel('Max Altitude (ft)')
PSPStyler (f, h, "Light");

f = figure(10);
h = plot(limit_OF, limit_expectedIsp,'b.');
grid on
title('ISP vs. OF')
xlabel('OF Ratio')
ylabel('Specific Impulse (s)')
PSPStyler (f, h, "Light");

f = figure(9);
h = plot(limit_OF, limit_expectedCStar,'b.');
grid on
title('C Star vs. OF')
xlabel('OF Ratio')
ylabel('C star (m/s)')
PSPStyler (f, h, "Light");

f = figure(8);
h = plot(limit_OF,limit_combustionTemperature,'b.');
grid on
title('Combustion Temperature vs. OF')
xlabel('OF Ratio')
ylabel('Combustion Temperature (K)')
PSPStyler (f, h, "Light");

f = figure(4)
h = plot(limit_burnTime,limit_maxAlt,'b.')
grid on
title("Maximum Altitude vs. Burn Time")
xlabel("Burn Time (s)")
ylabel("Maximum Altitude (ft)")
PSPStyler (f, h, "Light");

f = figure(5)
h = plot((limit_mDot ./ (1 + limit_OF)),limit_maxAlt,'b.')
grid on
title("Maximum Altitude vs. Fuel Mass Flow")
xlabel("Fuel Mass Flow Rate (lbs/s)")
ylabel("Maximum Altitude (ft)")
PSPStyler (f, h, "Light");

f = figure(6);
h = plot(limit_pressureExit,limit_maxAlt,'b.');
grid on
title("Maximium Altitude vs. Exit Pressure")
xlabel("Exit Pressure (psi)")
ylabel("Maximum Altitude (ft)")
PSPStyler (f, h, "Light");

f = figure(7);
h = plot(limit_AR,limit_maxAlt,'b.');
grid on
title("Maximum Altitude vs. Aspect Ratio")
xlabel("Aspect Ratio")
ylabel("Maximum Altitude (ft)")
PSPStyler (f, h, "Light");

%% Purdue Space Program - Liquids
%% Tango Zulu Package
% Function that alters an inputted plot into official PSP colors. Call this
% function after creation of the entire plot, grid, legend, and any
% additional curves.
%
% Inputs:   fig - figure [object], obtained by calling fig = figure(...);
%           plotIn - plot [object], obtained by calling plotIn = plot(...);
% Optional inputs:
%           colorMode - choose between "Light" and "Dark" modes when
%               graphing. Defaults to "Dark". [string]

function PSPStyler(fig, plotIn, colorMode)
    %% Constants Declaration
    lineWidth = 2;
    gold = '#DAAA00'; % curve 1
    dust = '#EBD99F'; % curve 2
    aged = '#8E6F3E'; % curve 3
    darkColor = '#252526';
    steel = '#555960'; % grid
    lightColor = '#F3F0E9'; % text
    railwayGray= '#9D9795';
    blue = '#4472C4';

    %% Figure Settings
    if nargin < 3
        colorMode = "Dark";
    end
    
    figAxes = gca;
    %dimmensions = numel(axis) / 2;
    plotIn.set('LineWidth',lineWidth)

    if strcmpi(colorMode, "light")
        figColor = lightColor;
        figAxesTitleColor = darkColor;
        figAxesColor = lightColor;
        figAxesXColor = darkColor;
        figAxesYColor = darkColor;
        figAxesZColor = darkColor;
        legendColor = lightColor;
        legendTextColor = darkColor;
        lineColors = {gold, railwayGray, darkColor};
    else
        figColor = darkColor;
        figAxesTitleColor = lightColor;
        figAxesColor = darkColor;
        figAxesXColor = lightColor;
        figAxesYColor = lightColor;
        figAxesZColor = lightColor;
        legendColor = darkColor;
        legendTextColor = lightColor;
        lineColors = {gold, dust, aged};
    end

    %% Figure Modification
    fig.Color = figColor;
    figAxes.Title.Color = figAxesTitleColor;
    figAxes.Color = figAxesColor;
    figAxes.XColor = figAxesXColor;
    figAxes.YColor = figAxesYColor;
    figAxes.ZColor = figAxesZColor;
    colororder(fig,lineColors);
    figAxes.GridColor = steel;
    figAxes.GridAlpha = 0.9;

    if ~isempty(figAxes.Legend)
        figAxes.Legend.Color = legendColor;
        figAxes.Legend.TextColor = legendTextColor;
    end
    
    set(fig, 'InvertHardCopy', 'off');
%     set(figAxes, 'FontName', 'Franklin Gothic Book');
%     set(gca, 'FontName', 'Franklin Gothic Book');
end

