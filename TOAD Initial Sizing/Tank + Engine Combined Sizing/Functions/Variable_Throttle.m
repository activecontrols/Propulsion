function [propellantMass, totalTime] = Variable_Throttle(start_TWR, mdot, min_throttle, maxThrust)
%% TOAD Flight
% Description: Flight Profile with Variable Throttle
% PD controlled profile follower targeting positions and velocities
% Author: Adam Grendys
% Last Edited: 3/5/2025

%% TOAD REFERENCES
% REQUIREMENTS DOC: https://docs.google.com/document/d/1jfazxSt6x4ROGItLOiyNnKDVktDiXMh2lE0mhNGMsWU/edit?usp=sharing
% SRR SLIDES: https://docs.google.com/presentation/d/151O5GhhcqatCP30IASsYGC5Nq8DB6PMOI8nIIVgjrB0/edit?usp=sharing
close all

%% INITIALIZATION
g_o = 32.174049; % (ft/sec^2)
rho_air = 0.0765; % (lb/ft^3)
Drag_Coef = .25; % Drag Coeffcient of TOAD (GUESS)
Cross_area = 5; % TOAD Cross sectional area (ft^2) (GUESS)
resolution = 1000;

tolerance = 0; % Acceptable Plus or Minus from referenced value (Deadband)
increment = 0; % Calculated change in throttle 
P_Gain = 0.01; % Gain for velocity error
D_Gain = 0.01; % Gain for acceleration error

%% CALCULATIONS
% Defining Profile Parameters
flightTime = 120; % Iterative Variable for Profile flight time (s)
time_hover = 5; % Hover time at apogee (s)
heightApogee = 165; % Apogee target height(ft)
UpSpeed = 15; % Ascent target speed (ft/s)
DownSpeed = -10; % Descent target speed (ft/s)
MaxLandingSpeed = 1.5; % Maximum allowable landing speed (ft/s) (1mph = ~1.5ft/s)
heightLanding = abs(DownSpeed); % Landing throttle taget start height (ft)
Status = 'Ascent'; % Flight Status String (Ascent, Hover, Descent, Landing)

indexHover = NaN;
indexDescent = NaN;
indexLanding = NaN;

% Defining Physics Parameters
mass_start = maxThrust / start_TWR;
mass_current = mass_start; % Iterative Variable for TOAD mass (lbm)
throttle = zeros(1, resolution); % Throttle percentage 40 to 100%
force = zeros(1, resolution); % (lbf)
position = zeros(1, resolution); % (ft)
velocity = zeros(1, resolution); % (ft/s)
acceleration = zeros(1, resolution); % (ft/s^2)
time = linspace(0, flightTime, resolution); % (s)

i = 1;
loopTrue = 1;
throttle(i) = 1 - increment;

while loopTrue
    i = i + 1;
    dt = time(i) - time(i - 1);

    % THROTTLE DETERMINATION
    switch Status
        case 'Ascent'
            %disp("ASCENT")
            increment = P_Gain * (UpSpeed - velocity(i - 1)) + D_Gain * (0 - acceleration(i - 1));

            if (velocity(i - 1) <= UpSpeed + tolerance) && (velocity(i - 1) >= UpSpeed - tolerance)
                throttle(i) = throttle(i - 1);
                disp("IN")
            else
                throttle(i) = throttle(i - 1) + increment;
                if increment > 0
                    %disp("UNDER")
                else
                    %disp("OVER")
                end
            end

            if position(i - 1) >= heightApogee - 8 * (velocity(i - 1) * dt - 1/2 * acceleration(i - 1) * dt^2)
                indexHover = i;
                Status = 'Hover';
            end
        case 'Hover'
            %disp("HOVER")
            increment = P_Gain * (0 - velocity(i - 1)) + D_Gain * (0 - acceleration(i - 1));

            if (velocity(i - 1) <= 0 + tolerance) && (velocity(i - 1) >= 0 - tolerance)
                throttle(i) = throttle(i - 1);
            else
                throttle(i) = throttle(i - 1) + increment;
                if increment > 0
                    %disp("UNDER")
                else
                    %disp("OVER")
                end
            end

            if time(i) - time(indexHover) >= time_hover
                indexDescent = i;
                Status = 'Descent';
            end
        case 'Descent'
            %disp("DESCENT")
            increment = P_Gain * (DownSpeed - velocity(i - 1)) + D_Gain * (0 - acceleration(i - 1));

            if (velocity(i - 1) <= DownSpeed + tolerance) && (velocity(i - 1) >= DownSpeed - tolerance)
                throttle(i) = throttle(i - 1);
            else
                throttle(i) = throttle(i - 1) + increment;
                if increment > 0
                    %disp("UNDER")
                else
                    %disp("OVER")
                end
            end

            if position(i - 1) <= heightLanding
                indexLanding = i;
                Status = 'Landing';
            end
        case 'Landing'
            %disp("LANDING")
            increment = P_Gain * (0 - velocity(i - 1)) + D_Gain * (0 - acceleration(i - 1));

            if (velocity(i - 1) <= 0 + tolerance) && (velocity(i - 1) >= 0 - tolerance)
                throttle(i) = throttle(i - 1) - dt/10;
            else
                throttle(i) = throttle(i - 1) + increment;
                if increment > 0
                    %disp("UNDER")
                else
                    %disp("OVER")
                end
            end

            if position(i - 1) <= 0 + tolerance && abs(velocity(i - 1)) < tolerance
                loopTrue = 0;
            end
    end

    if throttle(i) > 1
        throttle(i) = 1;
    elseif throttle(i) < min_throttle
        throttle(i) = min_throttle;
    end

    % PHYSICS LOOP
    mdot_current = mdot * throttle(i);
    thrust_current = maxThrust * throttle(i);
    mass_current = mass_current - mdot_current * dt;

    force(i) = (thrust_current - mass_current) - ... 
        1/2 * rho_air * velocity(i - 1)^2 * Drag_Coef * Cross_area;
    acceleration(i) = (force(i) / mass_current) * g_o;
    velocity = cumtrapz(time, acceleration);
    position = cumtrapz(time, velocity);

    % LOOP END CHECKS
    if position(i) <= 0 
        if abs(velocity(i)) > MaxLandingSpeed
            loopTrue = 0;
            fprintf("\nCRASH!")
        elseif abs(velocity(i)) < MaxLandingSpeed && Status == "Landing"
            loopTrue = 0;
            fprintf("\nWOAH GOOD!")
        end
    end

    if i > resolution - 1
        loopTrue = 0;
        fprintf("\nBOUNDS EXCEEDED")
    end
end

i = i - 1;
propellantMass = mass_start - mass_current;
totalTime = time(i);

%% FORMATTED OUTPUT
figure(1)
subplot(2,2,1)
sgtitle("TOAD Controlled Hop Profile")
plot(time(1:i), position(1:i))
xline(time(indexHover), 'r--')
xline(time(indexDescent), 'r--')
xline(time(indexLanding), 'r--')
ylabel("Position (ft)")
xlabel("Time (s)")
grid on

subplot(2,2,2)
plot(time(1:i), velocity(1:i))
xline(time(indexHover), 'r--')
xline(time(indexDescent), 'r--')
xline(time(indexLanding), 'r--')
xlabel("Time (s)")
ylabel("Velocity (ft/s)")
grid on

subplot(2,2,3)
plot(time(1:i), acceleration(1:i))
xline(time(indexHover), 'r--')
xline(time(indexDescent), 'r--')
xline(time(indexLanding), 'r--')
xlabel("Time (s)")
ylabel("Acceleration (ft/s^2)")
grid on

subplot(2,2,4)
plot(time(1:i), throttle(1:i) * 100)
xline(time(indexHover), 'r--')
xline(time(indexDescent), 'r--')
xline(time(indexLanding), 'r--')
xlabel("Time (s)")
ylabel("Throttle (%)")
grid on

fprintf("\nTotal Time: %.3f s", totalTime)
fprintf("\nPropellant Mass: %.3f lbm", propellantMass)

