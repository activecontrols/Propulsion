function [MaxAlt, rail_speed, Mach_num_max, MaxAccel] = ...
    Traj_1DoF_Model_Throttling_V3(m_dry, Cd, outerDiameter, thrust_lbf, t_b, isp_ex, mDotTotal, printResults, curve)

%{
Purdue Space Program - Liquids
Rocket 3 1DoF Trajectory Model (Throttling)
Tyler Trostle
Last Edit: 1/19/2022 (Tyler)

% THIS SCRIPT IS MEANT TO BE CALLED FROM OptimalThrustCurve.m SCRIPT

Inputs:
  1. Dry Mass: 'm_dry' [lb]
  2. Drag Coefficient: 'Cd' [assumed]
  3. Outer Diameter: 'outerDiameter' [in]
  4. Thrust: 'thrust_lbf' [lb]
  5. Burn Time: 't_b' [s]
  6. Expected Isp: 'isp_ex' [s]
  7. Total Mass Flow rate: ' mDotTotal [lb/s]
  8. Print Results: 'printResults' [boolean]
  9. minThrottle: minimum throttle percentage [0-1]
Outputs:
  1. Maximum Altitude: 'MaxAlt' [ft]
  2. Speed at end of rail: 'rail_speed' [m/s]
  3. Maximum Mach: 'Mach_num_max' [Mach]
  4. Maximum Accleration: 'MaxAccel' [g's]
%}
    

conversions;

%%EXAMPLE EXECUTION WITH BZB VALUES

% [MaxAlt, AltMaxV, Mach_num_max, MaxAccel] = Traj_1DoF_Model_Throttling_V3(100, .4, 6.5, 900, 10.3, 200, 4, 1, [1.3 1 1 1 1])

% Script Initialization
% c_order = get(gca, 'ColorOrder');
% figNum = 0;

% Define constants
g_0 = 9.8; % grav accel m/s^2

%% Define Rocket Properties

% Dry mass
m_dry = m_dry * c.LBM2KG; % convert to kg

% Reference area
S = pi * (outerDiameter / 2 * c.IN2M) ^ 2; % cross section area m^2

% Thrust
%thrust_lbf = 900; % lbf
thrust = thrust_lbf * c.LBF2N; % N

% Specific Impulse
%Isp = 200; % s

% Calculate total impulse (Max allowable = 9208 lbf*s)
impulse = thrust * t_b; % N*s

% Calculate propellant mass
mDotTotal = mDotTotal * c.LBM2KG;
m_prop = mDotTotal * t_b; % kg

% Calculate mass flow rate
m_dot_0 = mDotTotal; % kg/s

% Throttling Values
sum = 0;
n=10000;
for i = 1:n
    sum = sum + curve(ceil(size(curve,2)*i/n));
end
avgThrottle = sum/n;
t_b = impulse/(avgThrottle * thrust); %new burn time

if printResults
    % Print title
    fprintf('\n------ 1DOF SIMULATION ------\n\n')

    % Print input parameters
    fprintf('INPUT PARAMETERS\n')
    fprintf('Dry Mass  = %.2f lbm\n', m_dry * c.KG2LBM)
    fprintf('Cd        = %.2f\n', Cd)
    fprintf('Area      = %.4f m^2\n', S)
    fprintf('Burn Time = %.2f s\n', t_b)
    fprintf('Isp       = %.1f s\n', isp_ex)
    fprintf('Thrust    = %.1f lbf\n\n', thrust_lbf)

    % Print calculated parameters
    fprintf('CALCULATED PARAMETERS\n')
    fprintf('Impulse         = %.1f lbf*s\n', impulse * c.N2LBF)
    fprintf('Propellant Mass = %.2f lbm\n', m_prop * c.KG2LBM)
    fprintf('Total Mass      = %.2f lbm\n', (m_prop + m_dry) * c.KG2LBM)
    fprintf('Mass Flow Rate  = %.2f lbm/s\n\n', m_dot_0 * c.KG2LBM)
end

% Initial conditions
FAR_alt = 707; % altitude of FAR-MARS site (m)
rail_alt = FAR_alt + 60 * c.FT2M; % altitude of the end of the rail (m)
alt = FAR_alt;
v = 0; % m/s
Mach_num = 0;

% Time intervals
t = 0; % s
dt = 0.005; % time step for sim (s)

% Maximum values
v_max = 0; % m/s
Mach_num_max = 0;
acc_max = 0; % m/s^2
v_max_alt = 0;
Mach_num_max_alt = 0;

% Atmosphere Data
altStep = .1;
atmosAltMax = 84000; %max altitude for atmosphere data (above ~84000 meters values arent available)
[atmosData(:,2), atmosData(:,1), ~, atmosData(:,3)] = stdatmo(0:altStep:atmosAltMax);

%% Run simulation
if printResults
    fprintf('Simulation running...')
end
tic
index = 0;
while v >= 0
    index = index + 1;
    
    % Increment time
    t = t + dt; % s
    
    % Calculate total mass
    if t < t_b
        throttle = curve(ceil(size(curve,2)*t/t_b));
        thrustInst = thrust * throttle;
        thrustInstMat(index) = thrustInst;

        mDotTotal = m_dot_0 * throttle;
        m_dotMat(index) = mDotTotal;
    
        m_prop = m_prop - mDotTotal * dt; % kg
    else
        thrustInst = 0;
        thrustInstMat(index) = 0;

        mDotTotal = 0;
        m_dotMat(index) = 0;
        
        m_prop = 0; % kg
    end
    m = m_dry + m_prop; % kg
    
    % Determine speed of sound and density at current altitude
    if alt < atmosAltMax
        a = atmosData(round(alt / altStep), 1);
        rho = atmosData(round(alt / altStep), 2);
        pres = atmosData(round(alt / altStep), 3);
    else 
        a = 275.34;
        rho = 0;
        pres = 0;
    end
    %[~, a, ~, rho] = atmosisa(alt);
    %[~, a, ~, rho] = altcond1(alt);
    
    % Calculate drag force
    %drag = 0.5 * rho * v ^ 2 * Cd * S; % N
    if Mach_num < 0.3
        drag = 0.5 * rho * v * v * Cd * S; %incompressible dynamic pressure [Pa]
    else
        drag = 0.5 * Mach_num * Mach_num * 1.4 * pres * Cd * S; %compressible dynamic presure [Pa]
    end
    
    % Calculate instantaneous acceleration
    if t < t_b
        acc_i = ((thrustInst - drag) / m) - g_0; % m/s^2
    else
        acc_i = (-drag / m) - g_0; % m/s^2
    end
    
    acc = acc_i; % m/s^2
    
    % Determine if max acceleration has been achieved
    if acc > acc_max
        acc_max = acc; % m/s^2
    end
        
    % Integrate acceleration
    v = v + acc * dt; % m/s
    
    % Determine if max speed has been achieved
    if v > v_max
        v_max = v; % m/s
        v_max_alt = alt; % m
    end
    
    % Calculate Mach number
    Mach_num = v / a;
    
    % Determine if maximum Mach number has been achieved
    if Mach_num > Mach_num_max
        Mach_num_max = Mach_num;
        Mach_num_max_alt = alt; % m
    end
    
    % Integrate velocity
    alt = alt + v * dt;

    % Find speed at end of rail
    if (alt > (rail_alt - .5)) && (alt < (rail_alt + .5))
        rail_speed = v; %[m/s]
    end
end
run_time = toc;

%timeVec = 0:dt:t;

%% Print results
MaxAlt = alt * c.M2FT; %Max Alt in feet
AltMaxV = v_max_alt * c.M2FT; %Altitude of Max Vel (ft)
MachNumMaxAlt = Mach_num_max_alt * c.M2FT; %Max Mach Number
MaxAccel = acc_max / g_0; %Max acceleration in g's

if printResults
    fprintf('done.\n')
    fprintf('Simulation time = %.3f s\n\n', run_time)
    fprintf('RESULTS\n')
    fprintf('Max Altitude          = %.2f ft\n', MaxAlt)
    fprintf('Time to Apogee        = %.1f s\n', t);
    fprintf('Max Speed             = %.1f m/s\n', v_max)
    fprintf('Altitude at Max Speed = %.2f ft\n', AltMaxV)
    fprintf('Max Mach number       = %.2f \n', Mach_num_max)
    fprintf('Altitude at Max Mach  = %.2f ft\n', MachNumMaxAlt)
    fprintf('Max Acceleration      = %.2f g''s \n', acc_max/g_0)
    fprintf('Speed off of the rail = %.2f m/s \n\n', rail_speed)
end

%thrustInstMat * c.N2LBF
%m_dotMat * c.KG2LBM

