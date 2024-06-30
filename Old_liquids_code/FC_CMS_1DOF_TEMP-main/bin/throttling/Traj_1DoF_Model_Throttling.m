% Tyler Trostle and Jonah Fouts
% 11/06/2021

% NEED TO UPDATE WITH NEW TRAJ CODE

function [MaxAlt, AltMaxV, MachNumMaxAlt, MaxAccel] = ...
    Traj_1DoF_Model_Throttling(m_dry_lbm, Cd, R_OD, thrust_lbf, Isp, impulse_lbfs, printResults, minThrottle, minThrottleTime)

conversions;

%%EXAMPLE EXECUTION

%[MaxAlt, AltMaxV, MachNumMax, MaxAccel] = Traj_1DoF_Model_Throttling(64, .4, 6.5, 900, 200, 9208, 1, .8, .25);


%%INPUTS

%m_dry_lbm: dry mass in lbm
%Cd: coefficient of drag
%R_OD: rocket outer diameter (in)
%thrust_lbf: thrust of rocket (lbf)
%Isp: specific impulse (s)
%impulse_lbfs: total impulse (lbf * s)
%printResults: 0 (don't print results) or 1 (print results)

%minThrottle: minimum allowable throttle percentage (0-1)
%minThrottleTime: point in burn of lowest throttle (0-1)


%%OUTPUTS

%MaxAlt: max altitude (ft)
%AltMaxV: Altitude of Max Vel (ft)
%MachNumMax: Max Mach Number
%MaxAccel: Max acceleration in g's


%%BZB Values

%m_dry_lbm = 64
%Cd = 0.400;
%R_OD = 3.25*2
%thrust_lbf = 900
%Isp = 200;
%impulse_lbfs = 9208


% Script Initialization
% c_order = get(gca, 'ColorOrder');
% figNum = 0;


% Define constants
g_0 = 9.8; % grav accel m/s^2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Rocket Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dry mass
m_dry_lbm = m_dry_lbm * 1.3; % dry mass of rocket (lbm) - with mass growth factor of 1.3
m_dry = m_dry_lbm * c.LBM2KG; % convert to kg

% Drag coefficient
%Cd = 0.400;

% Reference area
S = pi * (R_OD / 2 * c.IN2M) ^ 2; % cross section area m^2

% Thrust
%thrust_lbf = 900; % lbf
thrust = thrust_lbf * c.LBF2N; % N
%minThrottle = .8;


% Specific Impulse
%Isp = 200; % s

% Total Impulse (Max allowable = 9208 lbf*s)
%impulse_lbfs = 9208; % lbf*s
impulse = impulse_lbfs * c.LBF2N; % N*s

% Calculate burn time
t_b = impulse / ((minThrottle + ((1 - minThrottle) / 2)) * thrust); % impulse over average thrust value for burn
minThrottleTime = minThrottleTime * t_b;

% Calculate propellant mass
m_prop_0 = impulse / (Isp * g_0); % kg
m_prop = m_prop_0;

% Calculate mass flow rate
m_dot = m_prop_0 / t_b; % kg/s

if printResults
    % Print title
    fprintf('\n------ 1DOF SIMULATION ------\n\n')

    % Print input parameters
    fprintf('INPUT PARAMETERS\n')
    fprintf('Dry Mass = %.2f lbm\n', m_dry_lbm)
    fprintf('Thrust   = %.1f lbf\n', thrust_lbf)
    fprintf('Isp      = %.1f s\n', Isp)
    fprintf('Impulse  = %.1f lbf*s\n', impulse_lbfs)
    fprintf('Cd       = %.2f\n', Cd)
    fprintf('S        = %.4f m^2\n\n', S)

    % Print calculated parameters
    fprintf('CALCULATED PARAMETERS\n')
    fprintf('Propellant Mass = %.2f lbm\n', m_prop_0 * c.KG2LBM)
    fprintf('Total Mass      = %.2f lbm\n', (m_prop_0 + m_dry) * c.KG2LBM)
    fprintf('Burn Time       = %.2f s\n', t_b)
    fprintf('Mass Flow Rate  = %.2f lbm/s\n\n', m_dot * c.KG2LBM)
end

% Initial conditions
alt = 707; % altitude of FAR-MARS site (m)
v = 0; % m/s
acc_old = 0; % m/s^2

% Time intervals
t = 0; % s
dt = 0.005; % time step for sim (s)

% Maximum values
v_max = 0; % m/s
Mach_num_max = 0;
acc_max = 0; % m/s^2

if printResults
    fprintf('Simulation running...')
end
tic
index = 0;
while v >= 0
    index = index + 1;
    
    % Increment time
    t = t + dt; % s
    
    thrustInst = thrust * thrustCurve2(t, t_b, minThrottle, minThrottleTime);
    thrustInstMat(index) = thrustInst;
    
    m_dot = thrustInst / (Isp * g_0);
    m_dotMat(index) = m_dot;
    
    % Calculate total mass
    if t < t_b
        m_prop = m_prop - m_dot * dt; % kg
    else
        m_prop = 0; % kg
    end
    m = m_dry + m_prop; % kg
    
    % Determine speed of sound and density at current altitude
    [~, a, ~, rho] = atmosisa(alt);
    %[~, a, ~, rho] = altcond1(alt);
    
    % Calculate drag force
    drag = 0.5 * rho * v^2 * Cd * S; % N
    
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
    v = v + acc*dt; % m/s
    
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
    alt = alt + v*dt;
    
end
run_time = toc;


%timeVec = 0:dt:t;

% Print results
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
    fprintf('Max Acceleration      = %.2f g''s \n\n', acc_max/g_0)
end

