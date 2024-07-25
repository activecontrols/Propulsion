function [] = optimizeMomentumRatios(radialDensity, annularDensity, radialVelocity, massFlowRateAnnular, massFlowRateRadial, throttleLevel)
    
    %{
       radialDensity: the density of the radial fluid, in kg/m^3
       annularDensity: the density of the annular fluid, in kg/m^3
       radialVelocity: the velocity of the radial fluid, in m/s
       annularVelocity: the velocity of the annular fluid, in m/s
    %}
    %defaultRadialDiameter = 0.05; %m
    %defaultAnnularGapDistance = 0.05; %m

    
    %radialDiameter = 0.01:0.01:0.10; % the diameter of the annular gap, in m
    %annularGapDistance = 0.01:0.01:0.10; % the distance of the annular gap, in m

    %LMRRadial = 1/4 .* (radialDensity .* pi .* radialDiameter.^2 .* radialVelocity^2) ./ (annularDensity .* radialDiameter .* defaultAnnularGapDistance .* annularVelocity^2); % Local Momentum ratio, varying the radial diameter of the pintle injector
    %LMRAnnular = 1/4 .* (radialDensity .* pi .* defaultRadialDiameter^2 .* radialVelocity^2) ./ (annularDensity .* defaultRadialDiameter .* annularGapDistance .* annularVelocity^2); % Local Momentum ratio, varying the annular gap distance of the pintle injector
    

    TMR = 1;
    outerRadius = 0.785/2; % m

    annularArea = (TMR .* massFlowRateAnnular.^2) ./ (annularDensity .* massFlowRateRadial .* radialVelocity); % m^2

    radialArea = (massFlowRateRadial .^ 2 * annularDensity .* annularArea)./(massFlowRateAnnular.^2 .*radialDensity .* TMR); % m^2

    slotHeight = radialArea ./ (0.025*0.0254*40); % m
    
    actuationLength = slotHeight - slotHeight(length(slotHeight)); % m

    throttleLevel = throttleLevel * 100; % (%)

    

    %Plots

    subplot(3,1,1)
    plot(throttleLevel, radialArea, '-or')
    grid on
    hold on
    xlabel("Throttle Level (%)")
    ylabel("Radial Area (m^2)")
    title('Throttle Level (%) vs. Radial Area (m^2)')
    hold off

    subplot(3,1,2)
    plot(throttleLevel, annularArea, '-ob')
    grid on
    hold on
    xlabel("Throttle Level (%)")
    ylabel("Annular Area (m^2)")
    title('Throttle Level (%) vs. Annular Area (m^2)')
    hold off
    
    subplot(3,1,3)
    plot(throttleLevel, actuationLength, '-og')
    grid on
    hold on
    xlabel("Throttle Level (%)")
    ylabel("Actuation Length (m)")
    title('Throttle Level (%) vs. Actuation Height (m)')
    hold off


