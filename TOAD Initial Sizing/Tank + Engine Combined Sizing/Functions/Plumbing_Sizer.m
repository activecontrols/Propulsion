clear
clc

% Display introductory message
disp("Hi! This is a calculator to calculate the pressure drop across the entire upper plumbing system of the TOAD test vehicle");

% Constants
Friction_factor = 0.21751;
Pipe_Inner_Dia = 0.0127;  % Pipe diameter in meters (0.5 inches)
g = 9.81; % Gravitational acceleration (m/s^2)
inch_to_meter = 0.0254; % Conversion factor  1 inch = 0.0254 m

% Fitting data from Table 8.4
fittings = {
"90° Elbow - Flanged regular", 0.3;
"90° Elbow - Flanged long radius", 0.2;
"90° Elbow - Threaded regular", 1.5;
"90° Elbow - Threaded long radius", 0.7;
"90° Elbow - Miter", 1.30;
"90° Elbow - Miter with vanes", 0.20;
"45° Elbow - Threaded regular", 0.4;
"45° Elbow - Flanged long radius", 0.2;
"Tee, Straight through flow - Threaded", 0.9;
"Tee, Straight through flow - Flanged", 0.2;
"Tee, Branching flow - Threaded", 2.0;
"Tee, Branching flow - Flanged", 1.0;
"Globe valve - Open", 10;
"Angle valve - Open", 5;
"Gate valve - Open", 0.20;
"Gate valve - 75% open", 1.10;
"Gate valve - 50% open", 3.6;
"Gate valve - 25% open", 28.8;
"Ball valve - Open", 0.5;
"Ball valve - 1/3 closed", 5.5;
"Ball valve - 2/3 closed", 200;
"Water meter", 7;
"Coupling", 0.08
};

% Display the list of fittings
disp("Available fittings and their K-values:");

for i = 1:length(fittings)
    fprintf("%d: %s (K = %.2f)\n", i, fittings{i,1}, fittings{i,2});
end

% User Inputs
Flow_Rate_lb_s = input("Enter mass flow rate through system (lb/s): ");
Flow_Rate_kg_s = Flow_Rate_lb_s * 0.4536;  % Convert lb/s to kg/s

% Fluid selection
Fluid_Choice = input("\nEnter choice of fluid for line system (1 = IPA, 0 = LOX): ");

if Fluid_Choice == 1
    Density = 786; % IPA density in kg/m^3
    Density_Imperial = Density / 27680;
else
    Density = 1141;  % LOX density in kg/m^3
    Density_Imperial = Density / 27680;
end

% Cross-sectional area and velocity calculations
Area_Cross_Section = pi * (Pipe_Inner_Dia / 2)^2;  % A = πr²
Volumetric_Flow_Rate = Flow_Rate_kg_s / Density;  % Q = ṁ / ρ
Pipe_Velocity = Volumetric_Flow_Rate / Area_Cross_Section;  % V = Q / A

% Initialize arrays for plotting
pressure_drop_points = [0]; % Start at zero pressure drop
component_labels = {"Start"};

% Number of pipes
Number_Pipes = input("Enter number of pipes in system: ");
Total_Pipe_Loss = 0;

if Number_Pipes > 0
    for i = 1:Number_Pipes

        fprintf("\nFor Pipe %d:\n", i);
        lengthOfPipe_in = input("Enter length of pipe (inches): ");
        lengthOfPipe_m = lengthOfPipe_in * inch_to_meter;

        % Pressure drop calculation for this pipe segment (converted to Psi)
        Pressure_Drop = ((Friction_factor * lengthOfPipe_m * (Pipe_Velocity^2) * Density) / (2 * Pipe_Inner_Dia)) / 6895;
        fprintf("Pressure drop for pipe segment %d: %.2f Psi\n", i, Pressure_Drop);

        Total_Pipe_Loss = Total_Pipe_Loss + Pressure_Drop;
        pressure_drop_points(end + 1) = Total_Pipe_Loss;
        component_labels{end + 1} = sprintf("Pipe %d", i);
    end
end

% Fittings consideration
Number_Fittings = input("\nEnter the number of fittings and bends to consider: ");
Total_Fitting_Loss = 0;

if Number_Fittings > 0

    for i = 1:Number_Fittings
        fprintf("\nFor Fitting / Bend %d:\n", i);
        Fitting_Index = input("Enter the index number of the fitting from the list: ");

        if Fitting_Index >= 1 && Fitting_Index <= length(fittings)
            K_Value = fittings{Fitting_Index, 2};
            fprintf("Selected: %s (K = %.2f)\n", fittings{Fitting_Index, 1}, K_Value);
        else
            disp("Invalid selection. Using default K = 0.");
            K_Value = 0;
        end

        % Calculate pressure loss for the fitting
        Head_Loss = (K_Value * Pipe_Velocity^2) / (2 * g);
        Fitting_Dp = (Head_Loss * g * Density) / 6895;
        fprintf("Pressure loss for fitting %d: %.2f Psi\n", i, Fitting_Dp);

        Total_Fitting_Loss = Total_Fitting_Loss + Fitting_Dp;
        pressure_drop_points(end + 1) = Total_Pipe_Loss + Total_Fitting_Loss;
        component_labels{end + 1} = sprintf("Fitting %d", i);
    end
end

% Regulators consideration
Number_Regulators = input("\nEnter the number of regulators in the system: ");
Total_Regulator_Loss = 0;

if Number_Regulators > 0
    for i = 1:Number_Regulators

        fprintf("\nFor Regulator %d:\n", i);
        Cv = input("Enter Cv for Throttle valve: \n") * 3.85 ;
        Regulator_Drop = (Flow_Rate_lb_s ^2) / ((Cv ^2) * Density_Imperial);
        fprintf("\nPressure drop for regulator %d: %.2f psi", i, Regulator_Drop);

        Total_Regulator_Loss = Total_Regulator_Loss + Regulator_Drop;
        pressure_drop_points(end + 1) = Total_Pipe_Loss + Total_Fitting_Loss + Total_Regulator_Loss;
        component_labels{end + 1} = sprintf("Regulator %d", i);
    end
end

% Final Output
Total_Pressure_Loss = Total_Pipe_Loss + Total_Fitting_Loss + Total_Regulator_Loss;
fprintf("\nTotal pressure drop due to all components: %.2f Psi\n", Total_Pressure_Loss);

% Plot the pressure drop along the system
figure;
plot(1:length(pressure_drop_points), pressure_drop_points, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xticks(1:length(component_labels));
xticklabels(component_labels);
xtickangle(45);
ylabel("Cumulative Pressure Drop (Psi)");
xlabel("Component");
title("Pressure Drop Along the System");
grid on;
