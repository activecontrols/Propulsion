function [countGUD, dataRAW] = dataWrite(dataRAW, fileNameRAW, fileNameGUD, countGUD)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Data Write Script
Nick Mondora

Input Variables:  dataRAW (matrix containing all data values)
                  fileNameRAW (Name of CSV containing all data points)
                  fileNameGUD (Name of CSV containing limited data points)

Output Variables: countGUD (count of acceptable rockets)
                  dataRAW (matrix containing all data values)
%}

sizeRAW = size(dataRAW);
dataGUD = []; % Var to contain all accepted rockets

for i = 1:sizeRAW(1)
    alt = dataRAW(i, 7);
    railSpeed = dataRAW(i, 8);
    
    req1 = alt > 50e3; % max altitude [ft]
    req2 = railSpeed > 100; % rail speed [ft/s]
    
    if (req1 && req2)
        dataGUD = cat(1, dataGUD, dataRAW(i,:));
        countGUD = countGUD + 1;
    end
end

% if there were any acceptable rockets, save them
if height(dataGUD) >= 1
    writematrix(dataGUD, fileNameGUD, 'WriteMode', 'append');
end

% save all rockets despite acceptibility, separate CSV file
writematrix(dataRAW, fileNameRAW, 'WriteMode', 'append');
dataRAW = zeros(1, sizeRAW(2));