%Stathis Kokkevis
%Toad preliminarythrottle profile calculation
%2/22/2025

clc
clear
close all;

%Define initial parameters

throttle = 0.1:0.01:1.0; %Throttle percentage
thrust = throttle * 550; %Thrust, in lbf
mdot = (1.484+1.236) * throttle; %Mass flow rate, constant O/F 1.2
TWR = 1.4; %Thrust to weight ratio of lander
fullyLoadedMass = thrust(end) / TWR; %Fully loaded mass of lander, in lbs

throttleDown = 0.2:0.01:(fullyLoadedMass / thrust(end)); % Max throttle to descend

timeAscent = 0:0.1:40; % Time interval for the ascent phase (seconds)
timeDescent = 40.1:0.1:100; % Time interval for the descent phase (seconds)


propellantMassFraction = 0.55; % Smaller than Masten Xodiac - Jan's calculation in slack
propellantMass = fullyLoadedMass * propellantMassFraction; % Propellant mass, in lbs


heightToad = zeros(1, (length(timeAscent) + length(timeDescent))); % Height of TOAD at every point in the mission
massToad = zeros(1, (length(timeAscent) + length(timeDescent)));



endAscentMass = 0; % Dummy variable to find maximum mass at end of ascent (least propellant mass used)

idealThrottleProfile = zeros(1, (length(timeAscent) + length(timeDescent))); %Ideal Throttle Profile for ascent


%Looping through all possible throttle, time of throttling up
%combinations to find best combination for fuel efficiency and target
%height
%Full means full throttle, down means throttled down
for i = 1:length(throttleDown)
    for j = 2:(length(timeAscent) - 1)
        
        %Define time interval of ascent as divided into full throttle
        %section and throttled down section
        ascentTimeIntervalFull = timeAscent(1:j);
        ascentTimeIntervalDown = timeAscent(j:end);

        % Calculate acceleration with mass flow rates of TADPOLE Engine (full
        % throttle)
        ascentMassFull = fullyLoadedMass - mdot(end) .* ascentTimeIntervalFull;
        ascentForceFull = thrust(end) - ascentMassFull;
        ascentAccelerationFull = ascentForceFull ./ascentMassFull;

        %Integrate to calculate velocities and resulting positions
        ascentVelocityFull = cumtrapz(ascentTimeIntervalFull, ascentAccelerationFull);
        ascentPositionFull = cumtrapz(ascentTimeIntervalFull, ascentVelocityFull);

        %Throttled down acceleration calculations
        ascentMassDown = ascentMassFull(end) - mdot(i) .* timeAscent(j:end);
        ascentForceDown = thrust(i) - ascentMassDown;
        ascentAccelerationDown = ascentForceDown ./ ascentMassDown;
        
        %Throttled down position/velocity
        ascentVelocityDown = ascentVelocityFull(end) + cumtrapz(ascentTimeIntervalDown, ascentAccelerationDown);
        ascentPositionDown = ascentPositionFull(end) + cumtrapz(ascentTimeIntervalDown, ascentVelocityDown);

        
        %Finding max height
        [maxHeight, idxMax] = max(ascentPositionDown);

        

        %Conditions for max height met? (50-55 m)
        if((164.042 < maxHeight) && (maxHeight < 180.4462))
            velocityAtMax = ascentVelocityDown(idxMax);
            
            %Conditions for velocity at apogee met?
            if (-1.5 < velocityAtMax) && (velocityAtMax < 1.5)
                
                %Highest mass calculated out of "accepted throttle profiles"
                if((ascentMassDown(end) > endAscentMass))
                    
                    %Assigning values for descent stage
                    endAscentMass = ascentMassDown(end);
                    endVelocity = ascentVelocityDown(end);
                    endPos = ascentPositionDown(end);
                    throttleEnd = throttleDown(i);
                    
                    %Ideal throttle profile found for ascent
                    idealThrottleProfile(1:j) = throttle(end);
                    idealThrottleProfile(j:length(timeAscent)) = throttleDown(i);
                    
                    %Height & mass of toad found throughout ascent
                    heightToad(1:j) = ascentPositionFull;
                    heightToad(j:length(timeAscent)) = ascentPositionDown;

                    massToad(1:j) = ascentMassFull;
                    massToad(j:length(timeAscent)) = ascentMassDown;

                end
            end
        end
    end
end



%Throttle level that will allow in ascent;
throttleUp = (thrust(1)/endAscentMass):0.01:1.0;

finalMass = 0;

cutoffTime = -1;

for i = 1:length(throttleUp)
    for j = 2:(length(timeDescent) - 1)

        %Define time interval of ascent as divided into full throttle
        %section and throttled down section
        descentTimeIntervalDown = timeDescent(1:j);
        descentTimeIntervalFull = timeDescent(j:end);

        % Calculate acceleration with mass flow rates of TADPOLE Engine (full
        % throttle)
        descentMassDown = endAscentMass - mdot(1) .* descentTimeIntervalDown;
        descentForceDown = thrust(1) - descentMassDown;
        descentAccelerationDown = descentForceDown ./descentMassDown;

        %Integrate to calculate velocities and resulting positions
        descentVelocityDown = endVelocity + cumtrapz(descentTimeIntervalDown, descentAccelerationDown);
        descentPositionDown = endPos + cumtrapz(descentTimeIntervalDown, descentVelocityDown);

        %Throttled down acceleration calculations
        descentMassFull = descentMassDown(end) - mdot(i) .* timeDescent(j:end);
        descentForceFull = thrust(i) - descentMassFull;
        descentAccelerationFull = descentForceFull ./ descentMassFull;

        %Throttled down position/velocity
        descentVelocityFull = descentVelocityDown(end) + cumtrapz(descentTimeIntervalFull, descentAccelerationFull);
        descentPositionFull = descentPositionDown(end) + cumtrapz(descentTimeIntervalFull, descentVelocityFull);


        %Finding max height
        [minHeight, idxMin] = min(descentPositionFull);
        
        
        

        %Conditions for min height met? (0-0.5 m)
        if((0 < minHeight) && (minHeight < 1.5))
            velocityAtMin = descentVelocityFull(idxMin);
            
            %Conditions for velocity at apogee met?
            if (-0.5 < velocityAtMin) && (velocityAtMin < 0.5)
                
                %Lowest mass calculated out of "accepted throttle profiles"
                if((descentMassFull(idxMin) > finalMass))
                    
                    

                    finalMass = descentMassFull(idxMin);
                    % %Assigning values for descent stage
                    % endMass = descentMassDown(end);
                    % endVelocity = descentVelocityDown(end);
                    % endPos = descentPositionDown(end);
                    % throttleEnd = throttleDown(i);

                    %Ideal throttle profile found for descent
                    idealThrottleProfile(length(timeAscent):length(timeAscent)+j) = throttle(1);
                    idealThrottleProfile(length(timeAscent)+j:end) = throttleUp(i);



                    %Height of toad found throughout descent
                    heightToad(length(timeAscent)+1:length(timeAscent)+j) = descentPositionDown;
                    heightToad(length(timeAscent)+j:end) = descentPositionFull;

                    massToad(length(timeAscent)+1:length(timeAscent)+j) = descentMassDown;
                    massToad(length(timeAscent)+j:end) = descentMassFull;

                    tolerance = 1e-6;
                    idx = find(abs(massToad - finalMass) < tolerance, 1);

                    cutoffTime = idx/10;
                    
                    
                    

                end
            end
        end
    end
end



if cutoffTime ~= -1
    time = 0:0.1:cutoffTime;
    
    
    figure;
    plot(time, idealThrottleProfile(1:idx+1), "-b", LineWidth=2)
    grid on;
    xlabel("Time (s)")
    ylabel("Throttle (%)")
    title('Time(s) vs. Throttle (%)')
    
    figure;
    plot(time, heightToad(1:idx+1), "-r", LineWidth=2);
    grid on;
    xlabel("Time (s)")
    ylabel("Altitude (feet)")
    title('Time (s) vs. Altitude (feet)')
    hold off
    
    propellantMassFinal = propellantMass - (fullyLoadedMass - finalMass);
    
    fprintf("Remaining Propellant Mass: %f lbs\n", propellantMassFinal);
else
    fprintf("You didn't close it buddy, try again (change input parameters)")
end

