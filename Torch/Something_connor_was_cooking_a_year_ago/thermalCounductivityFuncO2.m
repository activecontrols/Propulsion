% Written By Connor McGuire
%Generate a thermal conductivity value from the input temperature in kelvin
function k = thermalCounductivityFuncO2(t)
    t = t - 297.15;
    k = -2.2143*10^-8 * t^2 + 9.2414*10^-5 * t + 4.8000*10^-4;
end

% Source
% Thermal conductivity of hydrogen and oxygen as a function of temperature : https://www.researchgate.net/figure/Thermal-conductivity-of-hydrogen-and-oxygen-as-a-function-of-temperature_fig12_334170896