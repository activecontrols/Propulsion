function [dataReturn] = frootLoops(outputNum)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Loops Function
Nick Mondora, Cameron Williams, Talal Zaim, Jonah Fouts, Tyler Trostle,
Matt Currie, Nicholas Faber, Liam Schenk

Input file:   Input_Properties.xlsx

Output files: RAW-1DOF-mm-dd-yy-HH-MM.csv (Saves every run)
              GUD-1DOF-mm-dd-yy-HH-MM.csv (Saves rockets in bounds)

Input variables:  outputNum (number of outputted values)

Output variables: dataReturn (array of outputted values)
%}

%% INITIALIZATIONS AND GETTING STARTED
% ask for user running the script. will be added to CSV header
runDescript = input('Run Description: ', 's');

% initialize 2 CSV data logs: 1 for ALL raw data, 1 for rockets within
% acceptible ranges as defined in dataWrite.m
writeDate = datestr(now, 'mm-dd-yy-HH-MM-SS'); % determine date/time of script run
fileNameRAW = initializeDataLog('RAW-', writeDate, runDescript); % raw output file name
fileNameGUD = initializeDataLog('GUD-', writeDate, runDescript); % filtered output file name

% Load Excel file containing input parameters and assign to variables accordingly
inputProperties = readmatrix('Input_Properties.xlsx', 'sheet', 'StepSettings'); % contains iteration properties
tubeProperties = readmatrix('Input_Properties.xlsx', 'sheet', 'TubeSizes'); % contains tube properties
copvProperties = readmatrix('Input_Properties.xlsx', 'sheet', 'COPVSpecs'); % contains COPV properties

% chamber pressure range [psi]
pcLower = inputProperties(1, 2);
pcUpper = inputProperties(2, 2);
pcStep = inputProperties(3, 2);

% O/F range
lowerOF = inputProperties(1, 3);
upperOF = inputProperties(2, 3);
stepOF = inputProperties(3 ,3);

% burn time range [s]
lowerBurnTime = inputProperties(1, 4);
upperBurnTime = inputProperties(2, 4);
stepBurnTime = inputProperties(3, 4);

% m-dot range [lbs/s]
lowerMDotFu = inputProperties(1, 5);
upperMDotFu = inputProperties(2, 5);
stepMDotFu = inputProperties(3, 5);

% exit pressure [psi]
peLower = inputProperties(1, 6);
peUpper = inputProperties(2, 6);
peStep = inputProperties(3, 6);

% avionics mass [lbs]
mAvionics = inputProperties(1, 7);

% engine dry mass [lbs]
mEngine = inputProperties(1, 8);

% regen [boolean]
isRegen = inputProperties(1, 9);

% nozzle ablation rate [in/s]
ablationRate = inputProperties(1, 10);

% trajectory [boolean]
plotTraj = inputProperties(1, 11);

% pressure drops 
% tubeFrictionFactor = inputProperties(1, 12);
% injectorDrop = inputProperties(1, 13);
% regenDrop = inputProperties(1, 14);
% venturiDrop = inputProperties(1, 15);
lowerPlumbingDrop = inputProperties(1, 21);

%structures and stuff
% engineLength = inputProperties(1, 16);
% fincanLength = inputProperties(1, 17);
% midPlumLength = inputProperties(1, 18);
% recoveryLength = inputProperties(1, 19);
% avionicsLength = inputProperties(1, 20);

% tube sizes [in]
tubeODList = tubeProperties(1, 2 : end); % find list of different tube ODs
tubeWallThick = tubeProperties(2 : end, 2 : end); % get the thicknesses for all tubes

% COPV sizes
% copvID = copvProperties(:, 1);
pressureService = copvProperties(:, 2); % [psi]
copvLength = copvProperties(:, 3); % [in]
copvDiameter = copvProperties(:, 4); % [in]
mCopv = copvProperties(:, 5); % [lbs]
copvVolume = copvProperties(:, 6); % [L]

% atmos model [matrix]
atmosStep = 0.1; % atmos table altitude step [m]
atmosMax = 86000; % Max altitude that atmos table has data for [m]
[atmosModel(:, 1), atmosModel(:, 2), ~, atmosModel(:, 3)] = stdatmo(0 : atmosStep : atmosMax); % density, speed of sound curve

% calculate total number of iterations
countWritten = 0; % increments by 1 on every iteration. On iter. 100, triggers a csv write or mainOut
countCompleted = 0; % iterations completed
countGUD = 0; % number of acceptable rockets as specified in dataWrite.m
countEstimated = length(pcLower : pcStep : pcUpper) * length(lowerOF : stepOF : upperOF)...
    * length(lowerBurnTime : stepBurnTime : upperBurnTime) * length(lowerMDotFu : stepMDotFu : upperMDotFu)...
    * (size(tubeProperties, 2) - 1) * length(peLower : peStep : peUpper) + 1;

dataReturn = zeros(countEstimated, outputNum); % fixed size dataReturn matrix

% Generate progress bar
waitBar = waitbar(0, sprintf("0.00%% complete, INF s left, %.0f total runs",countEstimated), 'Name',' Simulation progress','CreateCancelBtn','setappdata(gcbf, ''canceling'', 1)');
setappdata(waitBar, 'canceling', 0);
tic; %start timer for time estimation

%% LOOPS
for pressureExit = peLower : peStep : peUpper
    for pressureChamber = pcLower : pcStep : pcUpper
        for OF = lowerOF : stepOF : upperOF
            %% PROPULSION
            nameString = strcat('PSP_BM1_estimates_pip_', num2str(int8(pressureChamber / pressureExit)), '_p_c_',...
                num2str(pressureChamber), '_O_F_', num2str(OF));

            inputName = append(nameString, '.inp');
            outputName = append(nameString, '.out');

            %run the CEA script to get expected (with efficiency loss) Isp [s]
            %, C* [m/s], expansionRatio and specificHeatRatio
            % [Isp, CStar, expansionRatio, specificHeatRatio, combustionTemperature, ~, ~, ~, ~] = PSP_1DOF_CEA_function_wrapper(pressureChamber,...
            %     pressureExit, OF, nameString, 0);
            efficiencyCStar = .9;
            efficiencyCF = .9;

            % Moves and deletes files
            % movefile(inputName, 'INP_OUT');
            % movefile(outputName, 'INP_OUT');
            % delete(append(pwd, '\PSP_CEA_function_wrapper\', inputName));

            for mDotFu = lowerMDotFu : stepMDotFu : upperMDotFu
                % run film cooling function
                mDot = (OF + 1) * mDotFu;
                %[fcPercent, fcDiam, fcHoles] = get_film_cooling_pcnt(OF, mDot, pressureChamber, true);
                fcDiam = 1;
                fcHoles = 1;
                fcPercent = .105;
                mDotFC = mDotFu * fcPercent;
                new_mDotFu = mDotFu + mDotFC;
                tank_OF = OF * mDotFu / new_mDotFu;

                % run CEA for propulsion evaluation
                [Isp, CStar, expansionRatio, specificHeatRatio, combustionTemperature, ~, ~, ~, ~] = PSP_1DOF_CEA_function_wrapper(pressureChamber,pressureExit, (fcPercent * tank_OF + (1 - fcPercent) * OF), nameString, 0);
                movefile(inputName, 'INP_OUT');
                movefile(outputName, 'INP_OUT');
                delete(append(pwd, '\PSP_CEA_function_wrapper\', inputName));
                expectedIsp = Isp / 9.81 * efficiencyCStar * efficiencyCF;
                expectedCStar = CStar * efficiencyCStar;

                for burnTime = lowerBurnTime : stepBurnTime : upperBurnTime
                    for tubeODIndex = 1 : length(tubeODList)
                        % Wait bar percent completed calculations & update
                        percentDone = countCompleted / countEstimated * 100;
                        percentLeft = 100 - percentDone;
                        if countCompleted ~= 0
                            avgIterTime = avgIterTime + (avgIterTime + toc) / countCompleted;
                        else
                            avgIterTime = 0;
                        end
                        timeLeft = seconds(percentLeft / (percentDone / avgIterTime));
                        timeLeft.Format = 'hh:mm:ss';

                        %% FLUID SYSTEMS
                        upScale = 1; % Multiply engine mass by this if using regen
                        mEngine = mEngine * (upScale ^ isRegen); % [lbm]

                        % tank thicknesses
                        tubeOD = tubeODList(tubeODIndex); % [in]
                        availThicknesses = tubeWallThick(:, tubeODIndex); % [in]

                        % run fluids script
                        [heightOx, heightFu, mTanks, mFluidSys, mFluids, volumeFuTank, volumeOxTank, ...
                            tubeIDFu, tubeIDOx, pressureFuTank, pressureOxTank, velOx, velFu, runLineOx, runLineFu, bestCopvIdx] = ...
                            tankSizing(burnTime, new_mDotFu, OF, pressureChamber, isRegen, availThicknesses, tubeOD, pressureService, copvVolume, mCopv, lowerPlumbingDrop);

                        if mFluidSys == intmax
                            % update wait bar
                            countCompleted = countCompleted + 1;

                            waitbar(percentDone / 100, waitBar, sprintf("%.2f%% complete, %s left, %.0f total runs", percentDone, timeLeft, countEstimated)); % update progress bar

                            if getappdata(waitBar,'canceling') % Check for clicked Cancel button
                                delete(waitBar);
                                error("Simulation cancelled by user after %.0f iterations.\nData may not have saved.", countCompleted);
                            end
                            continue
                        end

                        volumeHeTank = copvVolume(bestCopvIdx);

                        % Total height of combined tanks
                        tankHeight = heightOx + heightFu; % [in]
                        heightOx = heightOx / 12; % [ft]
                        heightFu = heightFu / 12; % [ft]

                        %% STRUCTURES
                        %run the structures script
                        [mStructures, rocketHeight] = structures(tubeOD, tankHeight, copvLength(bestCopvIdx));

                        % calculate aspect ratio
                        AR = rocketHeight / tubeOD;

                        % calculate rocket height [ft]
                        rocketHeight = rocketHeight / 12; % [in]

                        %% TRAJECTORY
                        printResults = false;

                        growthFactor = 1.3; % Account for additional mass
                        mDry = (mStructures + mAvionics + mFluidSys + mEngine) * growthFactor; % [lbm]
                        mWet = mDry + mFluids; % [lbm]

                        mDotTotal = mDotFC + mDotFu * (1 + OF); % [lbm/s]

                        thrust = expectedIsp * mDotTotal; % [lbf]

                        twr = thrust / mWet;

                        % Bad Rocket Checks
                        if (thrust > 4000 || AR > 35 || mDry > 300 || twr < 2.5 || runLineFu < 4 || runLineOx < 4 || copvDiameter(bestCopvIdx) > (tubeOD - 0.25) || tankHeight > 132)
                            % update wait bar
                            countCompleted = countCompleted + 1;

                            waitbar(percentDone / 100, waitBar, sprintf("%.2f%% complete, %s left, %.0f total runs", percentDone, timeLeft, countEstimated)); % update progress bar

                            if getappdata(waitBar,'canceling') % Check for clicked Cancel button
                                delete(waitBar);
                                error("Simulation cancelled by user after %.0f iterations.\nData may not have saved.", countCompleted);
                            end
                            continue
                        end

                        % call trajectory UDF
                        [maxAlt, railSpeed, maxMach, maxAccel, burnoutAlt, exitRadius, exitPresHalfBurn, altAtTenSeconds, radiusThroatInit, ~, ~, trajectoryArray] = ...
                            trajModel(mDry, mWet, tubeOD, thrust, burnTime, mDotTotal, printResults, atmosModel, ...
                            expectedCStar, pressureChamber, expansionRatio, pressureExit, ablationRate, specificHeatRatio, isRegen, plotTraj);

                        countCompleted = countCompleted + 1; % record that we've completed an additional iteration (waitBar)
                        countWritten = countWritten + 1; % ^same but this gets reset after every data write
                       
                        % save trajectory array to folder
                        if plotTraj
                            fileName = strcat('Trajectory Models\TrajArray_', num2str(countCompleted), '.mat');
                            save(fileName, 'trajectoryArray');
                        end
                        
                        throatAblation = ablationRate * burnTime;

                        % Log data
                        dataRAW(countWritten, :) = [tubeOD, mDry, mWet, mTanks, burnTime, expectedIsp, maxAlt, railSpeed, ...
                            maxMach, maxAccel, pressureChamber, OF, mDotTotal, thrust, rocketHeight, AR, ...
                            heightFu, heightOx, burnoutAlt, tubeIDFu, tubeIDOx, expectedCStar, throatAblation, exitRadius, twr, ...
                            volumeFuTank, volumeOxTank, volumeHeTank, exitPresHalfBurn, pressureFuTank, pressureOxTank, ...
                            velOx, velFu, runLineOx, runLineFu, pressureExit, combustionTemperature, altAtTenSeconds, mCopv(bestCopvIdx), radiusThroatInit, countCompleted, ...
                            fcPercent, fcDiam, fcHoles];

                        % Dump data:
                        if countWritten == 100
                            dataReturn = cat(1, dataReturn, dataRAW);
                            [countGUD, dataRAW] = dataWrite(dataRAW, fileNameRAW, fileNameGUD, countGUD);
                            countWritten = 0;
                        end

                        % update wait bar
                        waitbar(percentDone / 100, waitBar, sprintf("%.2f%% Complete, %s left, %.0f total runs", percentDone, timeLeft, countEstimated)); % update progress bar
                        if getappdata(waitBar, 'canceling') % Check for clicked Cancel button
                            delete(waitBar);
                            error("Simulation cancelled by user after %.0f iterations.\nData may not have saved.", countCompleted);
                        end
                    end
                end % end for innerDiameter
            end % end for mDotFu
        end % end for burnTime
    end % end for O_F
end % end for pressureChamber

%% FINISHING THE JOB
waitbar(countCompleted / countEstimated, waitBar, sprintf('Finishing up...')); % update progress bar

% if writeCount didn't reach 10, log any remaining data from mainOut
if countWritten ~= 0
    dataReturn = cat(1,dataReturn,dataRAW);
    [countGUD, ~] = dataWrite(dataRAW, fileNameRAW, fileNameGUD, countGUD);
end

% clean up the workspace
clearvars -except fileNameGUD fileNameRAW waitBar countCompleted countGUD dataReturn;
delete(waitBar);

fprintf("Total Iterations: %.0f\n", countCompleted); % total iterations
fprintf("Total Acceptable Iterations: %.0f\n", countGUD); % total iterations
fprintf("Raw Excel file: %s\n", fileNameRAW); % file name
fprintf("Reduced Excel file: %s\n", fileNameGUD); % file name
clear tEnd count_completed count_ACC waitBar fileName_RAW fileName_ACC;

% Restores Path Searches
restoredefaultpath;
end