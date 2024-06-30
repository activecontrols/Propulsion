%% Function definition
function mfp = averageFreePath(temp, pressure, averageDiameter)
    kB = 1.380649e-23 *symunit().J / symunit().K; %Boltzmaan constant
    mfp = kB * temp / (sqrt(2) * pi() * (averageDiameter ^ 2) * pressure);
end
