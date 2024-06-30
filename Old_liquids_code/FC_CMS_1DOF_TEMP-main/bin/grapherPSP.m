%% Creates a PSP-Branded Figure
% Inputs:   figNum - the figure number [integer]
%           figTitle - the figure title [string]
%           plotX1 - the x dataset that is plotted [double]
%           plotY1 - the y dataset that is plotted [double]
%           xAxisLabel - label for x-axis [string]
%           yAxisLabel - label for y-axis [string]
%       Optional inputs:
%           limMatrix - 1x4 or 4x1 matrix defining the limits of the x- 
%               and y- axes in this order [xLower, xUpper, yLower, yUpper]
%           plotX2 - second x dataset that is plotted, corresponding to
%               plotY2 [double]
%           plotY2 - second y dataset that is plotted, corresponding to
%               plotX2 [double]
%           plotX3 - third x dataset that is plotted, corresponding to
%               plotY3 [double]
%           plotY3 - third y dataset that is plotted, corresponding to
%               plotX3 [double]
%           legendNames - nx1 or 1xn string matrix (with n = 1, 2, or 3)
%               listing the legend entries (in order by Y1, Y2, Y3) 
%               [string]

function grapherPSP(figNum, figTitle, plotX1, plotY1, xAxisLabel, yAxisLabel, varargin)
%% Constants Declaration
lineWidth = 2;
gold = '#DAAA00'; % curve 1
dust = '#EBD99F'; % curve 2
aged = '#8E6F3E'; % curve 3
background = '#252526';
steel = '#555960'; % grid
lightColor = '#F3F0E9'; % text

%% Variable Input Declaration
limMatrix = NaN;
plotX2 = NaN;
plotY2 = NaN;
plotX3 = NaN;
plotY3 = NaN;
legendNames = NaN;

if ~isempty(varargin)
    for index = 1:2:length(varargin)
        if varargin{index} == "limMatrix"
            limMatrix = varargin{index+1};

        elseif varargin{index} == "plotX2"
            plotX2 = varargin{index+1};

        elseif varargin{index} == "plotY2"
            plotY2 = varargin{index+1};

        elseif varargin{index} == "plotX3"
            plotX3 = varargin{index+1};
       
        elseif varargin{index} == "plotY3"
            plotY3 = varargin{index+1};

        elseif varargin{index} == "legendNames"
            legendNames = varargin{index+1};
        end
    end
end

%% Figure Creation
fig = figure(figNum);
hold on; grid on;
plot(plotX1, plotY1, 'Color', gold, 'LineWidth', lineWidth);
xlabel(xAxisLabel); 
ylabel(yAxisLabel);

figureTitle = title(figTitle);
figAxes = gca;
fig.Color = background;
figureTitle.Color = lightColor;
figAxes.Color = background;
figAxes.YColor = lightColor;
figAxes.XColor = lightColor;
figAxes.GridColor = steel;
figAxes.GridAlpha = 0.9;

%% Optional Feature Addition
if ~isnan (limMatrix)
    xlim = [limMatrix(1), limMatrix(2)];
    ylim = [limMatrix(3), limMatrix(4)];
end

if ~isnan(plotX2)
    plot(plotX2, plotY2, 'Color', dust, 'LineWidth', lineWidth);

    if ~isnan(plotX3)
        plot (plotX3, plotY3, '-', 'Color', aged, 'LineWidth', lineWidth);
    end

    l = legend(legendNames, 'location', 'best', 'TextColor', lightColor);
%     l.Title.Color = 'red';
end