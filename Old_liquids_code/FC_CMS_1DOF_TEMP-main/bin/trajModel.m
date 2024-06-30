function [maxAlt, railSpeed, maxMach, maxAccel, burnoutAlt, exitRadius, exitPresHalfBurn, altAtTenSeconds, radiusThroatInit, maxQ, maxQAlt, trajectoryArray] = ...
          trajModel(mDry, mWet, tubeOD, thrust, burnTime, mDotTotal, printResults, atmosModel, ...
          expectedCStar, pressureChamber, expansionRatio, pressureExit, ablationRate, specificHeatRatio, isRegen, plotTraj)
%{
Purdue Space Program - Liquids
Rocket 3 1DoF - Trajectory Model
Tyler Trostle, Talal Zaim, Cameron Williams

Input variables:  mDry [lbm]
                  mWet [lbm]
                  tubeOD [in]
                  thrust [lbf]
                  burnTime [s]
                  mDotTotal [lbm/s]
                  printResults (boolean)
                  atmosModel (atmospheric model matrix)
                  expectedCStar [m/s]
                  pressureChamber [psi]
                  expansionRatio
                  pressureExit [psi]
                  ablationRate [in/s]
                  specificHeatRatio
                  isRegen
Output variables: maxAlt [ft]
                  railSpeed [ft/s]
                  maxMach
                  maxAccel [g]
                  burnoutAlt [ft]
                  exitRadius [in]
                  exitPresHalfBurn [psi]
                  altAtTenSeconds [ft]
                  trajectoryArray [s, ft, ft/s, ft/s^2, Mach, lbm, psi, lbf, lbf]
                  plotTraj [boolean]
%}

% Define constants
g = 9.81; % grav accel m/s^2

% Define conversions
lbmToKilos = 0.453592;
kilosToLbm = 1 / lbmToKilos;
psiToPa = 6894.8;
inchesToMeters = 0.0254;
metersToFeet = 3.2808;
lbfToNewtons = 4.44822;

%% Define Rocket Properties
% obtain cd versus mach data
cdArray = readmatrix ("CD CMS.CSV");
cdArray = cdArray (:, [1,4,5]); % get power on and power off cd

% Dry mass
mDry = mDry * lbmToKilos; % lbs to kg
mWet = mWet * lbmToKilos; % lbs to kg
mass = mWet; % kg

% Reference area
refArea = pi * (tubeOD / 2 * inchesToMeters) ^ 2; % m^2

% Thrust
thrustInit = thrust * lbfToNewtons; % N
thrust = thrustInit; % N

% Calculate total impulse
impulse = thrust * burnTime; % N*s

% Calculate mass flow rate
mDotTotal = mDotTotal * lbmToKilos; % kg/s

% Nozzle geometry 
pressureChamber = pressureChamber * psiToPa; % PSI to pascals
pressureExit = pressureExit * psiToPa; % PSI to pascals
areaThroat = (expectedCStar * mDotTotal) / (pressureChamber); % [m^2]
radiusThroatInit = sqrt(areaThroat / pi); % [m]
radiusThroat = radiusThroatInit;
areaExit = areaThroat * expansionRatio; % [m^2]
exitRadius = sqrt(areaExit / pi); % [m]
ablationRate = ablationRate * inchesToMeters;% in/s to m/s

%{
if printResults
    % Print title
    fprintf('\n------ 1DOF SIMULATION ------\n\n')

    % Print input parameters
    fprintf('INPUT PARAMETERS\n')
    fprintf('Dry Mass  = %.2f lbm (with growth factor)\n', mDry)
    fprintf('Cd        = %.2f\n', Cd)
    fprintf('Area      = %.4f m^2\n', refArea)
    fprintf('Burn Time = %.2f s\n', burnoutTime)
    fprintf('Isp       = %.1f s\n', expectedIsp)
    fprintf('Thrust    = %.1f lbf\n\n', thrust)

    % Print calculated parameters
    fprintf('CALCULATED PARAMETERS\n')
    fprintf('Impulse         = %.1f lbf*s\n', impulse * c.N2LBF)
    fprintf('Propellant Mass = %.2f lbm\n', (mDotTotal * burnoutTime) * c.KG2LBM)
    fprintf('Total Mass      = %.2f lbm\n', (mWet) * c.KG2LBM)
    fprintf('Mass Flow Rate  = %.2f lbm/s\n\n', mDotTotal * c.KG2LBM)
end
%}

% Initial conditions
altFAR = 615.09; % altitude of FAR-MARS site (m)
altRail = altFAR + 60 * 12 * inchesToMeters; % altitude of the end of the rail (m)
alt = altFAR;
vel = 0; % m/s
mach = 0;
acc = 0;
burnoutAlt = intmin; % m
isOffRail = false;

% Time intervals
t = 0; % s
dt = 0.01; % time step for sim (s)

% Maximum values
maxVel = intmin; % m/s
maxMach = intmin;
maxAcc = 0; % m/s^2
maxVelAlt = intmin; % m
maxMachAlt = intmin; % m
maxQ = intmin; % pa
maxQAlt = intmin; % m
exitMachMin = intmax; % mach num
exitMachMax = intmin; % mach num

% Atmosphere Data
atmosStep = 0.1; % m

%[atmosData(:,2), atmosData(:,1), ~, atmosData(:,3)] = stdatmo(0:altStep:atmosAltMax);
atmosData(:,2) = atmosModel(:,1); %density (kg/m^3)
atmosData(:,1) = atmosModel(:,2); %speed of sound (m/s)
atmosData(:,3) = atmosModel(:,3); %pressure (Pa)


% Trajectory Time Array [time, position, velocity, acceleration, Mach,
% total mass, ambient pressure, thrust, pressure thrust]

trajectoryArray = [0, alt, vel, acc, mach, mass, 0, thrustInit, 0]; % intialize array

%% Run simulation
if printResults
    fprintf('Simulation running...')
end

while vel >= 0
    try % Check if rocket is not higher than 86 km
        speedOfSound = atmosData(round(alt / atmosStep), 1);
        densityAtmos = atmosData(round(alt / atmosStep), 2);
        pressureAtmos = atmosData(round(alt / atmosStep), 3);
    catch % return negative rocket if rocket goes to space
        maxAlt = intmin;
        railSpeed = intmin;
        maxMach = intmin;
        maxAccel = intmin;
        burnoutAlt = intmin;
        return;
    end
    
    if (t - burnTime / 2) < (dt * 2)
       exitPresHalfBurn = pressureAtmos;
    end

    if t <= burnTime
        % Calculate total mass
        mass = mass - mDotTotal * dt; %kg

        % Nozzle Geometry
        areaThroat = pi * radiusThroat ^ 2;
        %radiusExit = radiusExit + ablationRate * dt;
        %areaExit = pi * radiusExit ^ 2;
        expansionRatio = areaExit / areaThroat;

        exitMach = sqrt(((specificHeatRatio + 1) / (specificHeatRatio - 1)) ^ ((specificHeatRatio + 1) / 2) * expansionRatio ^ (specificHeatRatio - 1));

        if exitMach > exitMachMax
            exitMachMax = exitMach;
        elseif exitMach < exitMachMin
            exitMachMin = exitMach;
        end

        if isRegen == 0
            radiusThroat = radiusThroat + ablationRate * dt;
            thrust = thrustInit * (exitMach / exitMachMax);
        end

        pressureThrust = areaExit * (pressureExit - pressureAtmos);
        burnoutAlt = alt;
        if mach < 0.01
            Cd = 0.4;
        else
            Cd = cdArray(uint32(mach / 0.01), 3);
        end

        if round(t, 3) == 10
            altAtTenSeconds = alt; % m
        end
    else
        pressureThrust = 0;
        thrust = 0;
        mass = mDry;

        if mach < 0.01
            Cd = 0.4;
        else
            Cd = cdArray(uint32(mach / 0.01), 2);
        end
    end
    
    % Calculate dynamic pressure & drag
    if mach < 0.3
        pressureDynamic = 0.5 * densityAtmos * vel ^ 2; %incompressible dynamic pressure [Pa]
    else
        pressureDynamic = 0.5 * mach ^ 2 * 1.4 * pressureAtmos; %compressible dynamic presure [Pa]
    end
    
    drag = Cd * pressureDynamic * refArea; % [N]

    % Kinematics
    acc = ((thrust + pressureThrust - drag) / mass) - g; % m/s^2
    vel = vel + (acc * dt); % m/s
    alt = alt + (vel * dt);
    mach = vel / speedOfSound;
    
    % Checks for newest max alt, max vel, etc
    if (alt > altRail && isOffRail == false)
        railSpeed = vel * metersToFeet; % [ft/s]
        isOffRail = true;
    end
    
    if vel > maxVel
        maxVel = vel; % m/s
        maxVelAlt = alt; % m
    end
    
    if mach > maxMach
        maxMach = mach;
        maxMachAlt = alt; % m
    end

    if pressureDynamic > maxQ
        maxQ = pressureDynamic; % Pa
        maxQAlt = alt; % m
    end
    
    if abs(acc) > abs(maxAcc)
        maxAcc = acc; % m/s^2
    end

    if round(t, 3) == 10
        altAtTenSeconds = alt; % m
    end

    % Save values to Trajectory Array
    if plotTraj
      trajectoryArray = [trajectoryArray; t, alt, vel, acc, mach, mass, pressureAtmos, thrust, pressureThrust];
    end
    
    % Increment time and index
    t = t + dt; % s
end

%% Print results
maxAlt = alt * metersToFeet; % Max Alt in feet
maxVelAlt = maxVelAlt * metersToFeet; % Altitude of Max Vel (ft)
maxMachAlt = maxMachAlt * metersToFeet; % Max Mach Number
maxQAlt = maxQAlt * metersToFeet; % Max Q altitude [ft]
maxAccel = maxAcc / g; % Max acceleration in g's
exitPresHalfBurn = exitPresHalfBurn / psiToPa; %Pascals to PSI
exitRadius = exitRadius * metersToFeet * 12; % Meters to inches
radiusThroatInit = radiusThroatInit * metersToFeet * 12; % Meters to inches
burnoutAlt = burnoutAlt * metersToFeet; % Meters to Feet
altAtTenSeconds = altAtTenSeconds * metersToFeet; % altitude at ten seconds into flight [ft]

if plotTraj
    trajectoryArray(:, [2 3 4]) = trajectoryArray(:, [2,3,4]) .* metersToFeet; % convert spatial values to ft
    trajectoryArray(:, 6) = trajectoryArray(:, 6) ./ lbmToKilos; % convert mass to lbm
    trajectoryArray(:, 7) = trajectoryArray(:, 7) ./ psiToPa; % convert atmos pressure to psi
    trajectoryArray(:, 8) = trajectoryArray(:, 8) ./ lbfToNewtons; % convert Newtons to lb-f
    trajectoryArray(:, 9) = trajectoryArray(:, 9) ./ lbfToNewtons; % convert Newtons to lb-f

end

if printResults
    fprintf('done.\n')
    fprintf('Simulation time = %.3f s\n\n', toc)
    fprintf('RESULTS\n')
    fprintf('Max Altitude          = %.2f ft\n', maxAlt)
    fprintf('Time to Apogee        = %.1f s\n', t);
    fprintf('Max Speed             = %.1f m/s\n', maxVel)
    fprintf('Altitude at Max Speed = %.2f ft\n', maxVelAlt)
    fprintf('Max Mach number       = %.2f \n', maxMach)
    fprintf('Altitude at Max Mach  = %.2f ft\n', maxMachAlt)
    fprintf('Max Acceleration      = %.2f g''s \n', maxAcc/g)
    fprintf('Speed off of the rail = %.2f m/s \n', railSpeed)
    fprintf('Max Q                 = %.2f Pa \n', maxQ)
    fprintf('Altitude at Max Q     = %.2f m \n', maxQAlt)
    fprintf('Max Exit Velocity     = %.2f mach \n', exitMachMax)
    fprintf('Min Exit Velocity     = %.2f mach \n\n', exitMachMin)
end