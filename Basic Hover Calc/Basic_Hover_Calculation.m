%% Baasic Hover Calculations
%--------------------------------------------------------------------------
% Definition: Used to see the impact of different variables on the fuel
% usage during a fixed trajectory flight to 52 meters and back
%--------------------------------------------------------------------------
% Sections are as followed ->
% Accend: Acceleration, Constant Thrust, Deceleration
% Hover: Constant Thrust
% Decend: Deceleration, Constant Thrust, Acceleration
%--------------------------------------------------------------------------
% Written by Ryan Lukow (rlukow5@gmail.com)
% Created: 03/23/2023
% Last Updated: 04/19/2023
%--------------------------------------------------------------------------

clear;
clc;

%% Input Values

%Standard Constants
gravity = 9.80665 ; %m/s^2
ShutOff = false; %State of engine
ResultsMatrix = zeros(1);

Spreadsheet = readmatrix('Basic Hover Inputs.xlsx');

%Input from excel doc
DryMass = Spreadsheet(4,4); %Kg                
WetMass = Spreadsheet(3,4); %kg                              
HeightLimit = Spreadsheet(6,2); %meters                      
ascendSpeed = Spreadsheet(5,2); %m/s                          
TimeStep = Spreadsheet(7,2); %Sec                          
MaxThrust = Spreadsheet(2,4); %N
IndividualPlots = Spreadsheet(13,2); %T/F
IterateOverSafetyFactor = Spreadsheet(11,2); %T/F
IterateOverDryWeight = Spreadsheet(12,2); %T/F

CurrentDryWeight = DryMass;
x = 0;
VelocityError = 0;
%% Dry Weight Loop
while(VelocityError == 0)
x = x + 1;

CurrentMaxFuel = WetMass;
z = 0;
FactorOfSafety = 1;
%% Fuel vs Factor Of Saftey Loop
while(FactorOfSafety >= 1)
z = z + 1;

    clear WetWeight
    clear MassFlow

    %Determine remaining variables
    DryWeight = CurrentDryWeight * gravity; %N
    WetWeight(1) = CurrentMaxFuel * gravity; %N
    MassFlow(1) = 0;
    descendSpeed = -1 * ascendSpeed; %m/s   
    VelocityError = 0;
    accelerationTime = 8;
    accelerationFlightValue = ascendSpeed / accelerationTime;
    deltaX = 0.5 * accelerationFlightValue * accelerationTime^2;
    ConstantThrustTime = (HeightLimit - (2 * deltaX)) / ascendSpeed;
    
    %Center of mass calculations
    DensityLox = 1140; %Kg/m^3
    DensityIpa = 785; %Kg/m^3
    TankVolume = 0.037067; %m^3
    radiusCap = 0.1524; %m
    heightCyl = 0.3048; %m
    TankVolumeCap = ((2 / 3) * pi * radiusCap^3); %m^3
    TankVolumeCylinder = (pi * radiusCap^2 * heightCyl); %m^3
    distanceBetweenTanks = 0.1524; %m = 6 inches
    distanceToIpa = 2 * radiusCap + heightCyl + distanceBetweenTanks;
    
    %Initial values for Loop
    Height(1) = 0;
    Velocity(1) = 0;
    Acceleration(1) = 0;
    Force(1) = 0;
    Time(1) = 0;
    Index = 1;
    TotalWeight(1) = WetWeight(1) + DryWeight;
    TotalMass(1) = TotalWeight(1) / gravity;
    
    
    %% Trajectory Loop
    while(Height(Index) > 0 || Index < 100)
        
        %Change time
        Index = Index + 1;
        Time(Index) = Time(Index - 1) + TimeStep;
    
    
        %% Set throttle values
        %Initial accelerating ascent
        if Time(Index) <= accelerationTime
             throttle(Index) = (((ascendSpeed / accelerationTime) * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        %Constant velocity ascent
        elseif Time(Index) <= ConstantThrustTime + accelerationTime
            throttle(Index) = ((0 * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        %deceleration ascent
        elseif Time(Index) <= ConstantThrustTime + 2 * accelerationTime
            throttle(Index) = (((descendSpeed / accelerationTime) * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
       
        %Hover
        elseif Time(Index) <= ConstantThrustTime + 2 * accelerationTime + 10
            throttle(Index) = ((0 * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        
        %deceleration descent
        elseif Time(Index) <= ConstantThrustTime + 3 * accelerationTime + 10
            throttle(Index) = (((descendSpeed / accelerationTime) * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        %Constant velocity descent
        elseif Time(Index) <= 2 * ConstantThrustTime + 3 * accelerationTime + 10
            throttle(Index) = ((0 * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        %acceleration descent for landing
        elseif Time(Index) <= 2 * ConstantThrustTime + 4 * accelerationTime + 10
            throttle(Index) = (((ascendSpeed / accelerationTime) * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        %if not any above, throttle = 0
        else
            throttle(Index) = ((-0.2 * TotalMass(Index - 1)) + (TotalWeight(Index - 1))) / MaxThrust;
        end
    
        if(throttle(Index) > 1 || throttle(Index) < 0 )
            VelocityError = 1;
        end

        %% Update values
        %Massflow, ChamberPressure, & ExitPressure calculations
        %Lines of best fit provided based on Kamon's throttle code.
        MassFlow(Index) = 0.453592 * (2.45 * (throttle(Index)) + 0.29);
        ChamberPressure(Index) = 246.55 * (throttle(Index)) + 28.49;
        ExitPressure(Index) = 15.95 * (throttle(Index)) + 2.06;
    
        %Update values
        WetWeight(Index) = WetWeight(Index - 1) - abs((MassFlow(Index) * TimeStep * gravity));
        TotalWeight(Index) = WetWeight(Index) + DryWeight;
        TotalMass(Index) = TotalWeight(Index) / gravity;
        Force(Index) = (MaxThrust * throttle(Index)) - TotalWeight(Index);
        Acceleration(Index) = Force(Index) / TotalMass(Index);
        Velocity(Index) = Velocity(Index - 1) + Acceleration(Index) * TimeStep;
        Height(Index) = Height(Index - 1) + Velocity(Index) * TimeStep;
    
        %If the rockets landed, then 0 values
        if(Height(Index) <= 0 && Index > 2)
            Acceleration(Index) = 0;
            Velocity(Index) = 0;
            Height(Index) = 0;
        end
    
        %If the rockets runs out of fuel record shut off point
        if(WetWeight(Index) <= 0 && ShutOff == false)
            ShutOff = true;
            TimeOfShutoff = Time(Index);
            HeightOfShutoff = Height(Index);
            VelocityOfShutoff = Velocity(Index);
            AccelerationOfShutoff = Acceleration(Index);
            ForceOfShutoff = Force(Index);
        end
       
        %% Center of mass calculation
        WeightLox(Index) = (WetWeight(Index) / 2.3) * 1.3;
        WeightIpa(Index) = WetWeight(Index) - WeightLox(Index);
        MassLox(Index) = WeightLox(Index) / gravity;
        MassIpa(Index) = WeightIpa(Index) / gravity;
        VolumeLox(Index) = MassLox(Index) / DensityLox; %m^3
        VolumeIpa(Index) = MassIpa(Index) / DensityIpa; %m^3
    
        if(VolumeLox(Index) >= TankVolumeCylinder + TankVolumeCap) %in top cap
            V = VolumeLox(Index) - (TankVolumeCylinder + TankVolumeCap);
            a = 1;
            b = 0;
            c = 3 * radiusCap^2;
            d = ((2 * radiusCap^3) - (3 * V/ pi));
            HeightLoxRoots = roots([a b c d]);
            LoxRoot = 0;
            for i = 1:HeightLoxRoots(end)
                if(isreal(HeightLoxRoots(i)) && HeightLoxRoots(i) > 0)
                    LoxRoot = HeightLoxRoots(i);
                end
            end
            HeightLox(Index) = (heightCyl + radiusCap) + LoxRoot;
        elseif((VolumeLox(Index) >= TankVolumeCap)) %in cylinder
            V = VolumeLox(Index) - TankVolumeCap;
            HeightLox(Index) = (radiusCap) + (V / (pi * radiusCap^2));
        else %in bottom cap
            V = VolumeLox(Index);
            a = 1;
            b = 0;
            c = 3 * radiusCap^2;
            d = (3 * V/ pi);
            HeightLoxRoots = roots([a b c d]);
            LoxRoot = 0;
            for i = 1:HeightLoxRoots(end)
                if(isreal(HeightLoxRoots(i)) && HeightLoxRoots(i) > 0)
                    LoxRoot = HeightLoxRoots(i);
                end
            end
            HeightLox(Index) = LoxRoot;
    
        end
        
        if(VolumeIpa(Index) >= TankVolumeCylinder + TankVolumeCap) %in top cap
            V = VolumeIpa(Index) - (TankVolumeCylinder + TankVolumeCap);
            a = 1;
            b = 0;
            c = 3 * radiusCap^2;
            d = ((2 * radiusCap^3) - (3 * V/ pi));
            HeightIpaRoots = roots([a b c d]);
            IpaRoot = 0;
            for i = 1:HeightIpaRoots(end)
                if(isreal(HeightIpaRoots(i)) && HeightIpaRoots(i) > 0)
                    IpaRoot = HeightIpaRoots(i);
                end
            end
            HeightIpa(Index) = (heightCyl + radiusCap) + IpaRoot;
        elseif((VolumeIpa(Index) >= TankVolumeCap)) %in cylinder
            V = VolumeIpa(Index) - TankVolumeCap;
            HeightIpa(Index) = (radiusCap) + (V / (pi * radiusCap^2));
        else %in bottom cap
            V = VolumeIpa(Index);
            a = 1;
            b = 0;
            c = 3 * radiusCap^2;
            d = (3 * V/ pi);
            HeightIpaRoots = roots([a b c d]);
            IpaRoot = 0;
            for i = 1:HeightIpaRoots(end)
                if(isreal(HeightIpaRoots(i)) && HeightIpaRoots(i) > 0)
                    IpaRoot = HeightIpaRoots(i);
                end
            end
            HeightIpa(Index) = IpaRoot;
        end
    
        %Center of mass height
        CenterIpa(Index) = (HeightIpa(Index) / 2);
        CenterLox(Index) = (HeightLox(Index) / 2);
        
        CenterOfMass(Index) = (((distanceToIpa + CenterIpa(Index)) * MassIpa(Index) + (CenterLox(Index) * MassLox(Index))) / (MassIpa(Index) + MassLox(Index)));
    
    end
    
    %Set default values
    throttle(1) = throttle(2);
    MassFlow(1) = MassFlow(2);
    ExitPressure(1) = ExitPressure(2);
    ChamberPressure(1) = ChamberPressure(2);

    %record total fuel usage
    FuelUsage = abs(WetWeight(end) - WetWeight(1)) / gravity;
    FuelUsagePercent = (FuelUsage / CurrentMaxFuel) * 100;
    CenterChange = CenterOfMass(2) - CenterOfMass(end);
    
    FactorOfSafety = 100 / FuelUsagePercent;
    
        if(IterateOverSafetyFactor == 1)
        SafetyFactorMatrix(z) = FactorOfSafety;
        MaxFuelMatrix(z) = CurrentMaxFuel * 2.20462;
        CurrentMaxFuel = CurrentMaxFuel - 1;
        else
            FactorOfSafety = 0;
        end
end
   
    if(IterateOverDryWeight == 1)
        SafetyFactorMatrix2(x) = 100 / FuelUsagePercent;
        DryWeightMatrix(x) = CurrentDryWeight * 2.20462;
        CurrentDryWeight = CurrentDryWeight + 4.53592;  
    else
        VelocityError = 1;
    end

    if((CurrentDryWeight + WetMass) * gravity >= MaxThrust * 0.8)
        VelocityError = 1;
    end

end

if(IterateOverDryWeight == 0)
    VelocityError = 0;
end

%% Plots


if(IndividualPlots == 1)

    %Subplot of force, velocity, acceleration, and height
    figure(1)
    subplot(2,2,1)
    grid on
    sgtitle("")
    plot(Time, Height)
    grid on
    if(ShutOff == true)
        hold on
        plot(TimeOfShutoff,HeightOfShutoff,'r*')
        hold off
    end
    title("Height vs Time")
    xlabel("Time [s]")
    ylabel("Height [m]")
    
    subplot(2,2,2)
    plot(Time, Force)
    grid on
    title("Force vs Time")
    if(ShutOff == true)
        hold on
        plot(TimeOfShutoff,ForceOfShutoff,'r*')
        hold off
    end
    xlim([TimeStep (length(Time) * TimeStep) ])
    xlabel("Time [s]")
    ylabel("Force [N]")
    
    subplot(2,2,3)
    plot(Time, Acceleration)
    grid on
    title("Acceleration vs Time")
    if(ShutOff == true)
        hold on
        plot(TimeOfShutoff,AccelerationOfShutoff,'r*')
        hold off
    end
    xlim([TimeStep (length(Time) * TimeStep) ])
    xlabel("Time [s]")
    ylabel("Acceleration [m/s^2]")
    
    subplot(2,2,4)
    plot(Time, Velocity)
    grid on
    if(ShutOff == true)
        hold on
        plot(TimeOfShutoff,VelocityOfShutoff,'r*')
        hold off
    end
    title("Velocity vs Time")
    xlabel("Time [s]")
    ylabel("Velocity [m/s]")

    %Chamber Pressure over time
    figure(3);
    plot(Time,ChamberPressure);
    grid on
    title("Chamber Pressure of the Lander Throughout Flight")
    xlabel("Time [s]")
    ylabel("Chamber Pressure [psi]")

    %Exit Pressure over time
    figure(4);
    plot(Time,ExitPressure);
    grid on
    title("Exit Pressure of the Lander Throughout Flight")
    xlabel("Time [s]")
    ylabel("Exit Pressure [psi]")
    
    %Mass Flow over time
    figure(2);
    plot(Time, MassFlow);
    grid on
    title("Mass flow rate of the Lander Throughout Flight")
    xlabel("Time [s]")
    ylabel("Mass flow rate [kg/s]")

    %Throttle over time
    figure(5);
    plot(Time,throttle, "-");
    grid on
    title(" Throttle percentage of the Lander Throughout Flight")
    xlabel("Time [s]")
    ylabel("Percentage of Max thrust [%]")
    ylim([0 1])

    %Height, Velocity, and Acceleration over time
    figure(6)
    plot(Time, Height, 'b')
    grid on
    hold on
    plot(Time, Velocity, 'm')
    plot(Time, Acceleration, 'g')
    if(ShutOff == true)
        xline(TimeOfShutoff, 'r--')
        legend("Height [m]", "Velocity [m/s]", "Acceleration [m/s^2]", "Time of Engine Shutoff")
    else
        legend("Height [m]", "Velocity [m/s]", "Acceleration [m/s^2]")
    end
    hold off
    xlabel("Time [s]")

end

if(IterateOverSafetyFactor == 1)
    figure(7)
    plot(MaxFuelMatrix, SafetyFactorMatrix)
    grid on
    xlabel("Initial Fuel Weight [lbs]")
    ylabel("Factor of Saftey [%]")
    title("Initial Fuel Weight vs Factor of Saftey With a Dry Mass of " + (DryMass / 0.45359237) + " lbs")
    xline(123, 'r--')
    legend("", "Maximum Fuel Weight for Tanks", 'Location','northwest')
end

if(IterateOverDryWeight == 1)
    figure(8)
    plot(DryWeightMatrix, SafetyFactorMatrix2)
    grid on
    xlabel("Initial Dry Weight [lbs]")
    ylabel("Factor of Saftey [%]")
    title("Initial Dry Weight vs Factor of Saftey With a Wet Mass of " + round(WetMass * 2.20462) + " lbs")
end
%% Results

if(WetWeight(end) < 0)
    fprintf("Error: Fuel Ran Out\n")
elseif(VelocityError == 1)
    fprintf("Error: Ascend speed too large\n")
else
    fprintf("Total Fuel Usage: %f/%f Kg\n" ,FuelUsage, WetMass)
    fprintf("Total Fuel Usage: %.3f%%\n" ,FuelUsagePercent)
    fprintf("Time of Flight: %.2f seconds\n", Time(end))
    fprintf("Velocity at landing: %.2f m/s\n", Velocity(end - 1))
    fprintf("Center of mass change: %.4f m\n", CenterChange)
end

%If we use 95% of our fuel, the center of mass of the fuel changes by
%around 8.7 inches.

%If we use 24% of our fuel in a typical flight, the center of mass changes
%by only around 1.75 inches



