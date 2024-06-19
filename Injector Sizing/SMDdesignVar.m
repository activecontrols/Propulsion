function [D32] = SMDdesignVar(lOpen, rhoOx, Ufuel, Uox, surfaceTensionFuel, pintleTipAngles)

    %{

    lOpen: the opening distance of the pintle, from sleeve to tip (m)
    rhoOx: density of oxidizer (kg/m^3)
    Ufuel: speed of fuel (m/s)
    Uox: speed of oxidizer (m/s)
    surfaceTensionFuel: surface tension of fuel (N/m)

    %}
    
    
    xi = (90 - pintleTipAngles) / 90; % Normalized angle
    
    q = 3.455 - 0.225 * xi; % Empirical relation for q
    
    Weber = (rhoOx * lOpen * ((Uox - Ufuel)^2)) / surfaceTensionFuel; % Weber Number of fuel-oxidizer impingement

    D32 = lOpen .* ((1 ./ xi) .* exp(4 - q .* (Weber ^ (0.1)))); % Sauter Diameter


    
    %{
minValue = inf;
    minIndex = 0;
    
    for x = 1:numel(pintleTipAngles)
        if D32(x) < minValue
            minValue = D32(x);
            minIndex = x;
        end
    end
    
    pintleTipAngle = pintleTipAngles(minIndex);
    
    fprintf("Minimum Sauter Diameter: %0.4f m at %d degrees\n", minValue, pintleTipAngle);
end
    %}