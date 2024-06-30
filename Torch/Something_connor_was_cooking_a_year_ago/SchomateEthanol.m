function solution = SchomateEthanol(temps, result) 
    %Temps in K (Scalar, not unit-ed by symunit) should be a single variabl
                            %A             B         C           D         E          F      G       H
    schomateCoefficients = [-8.87256 , 282.389 , -178.85 ,   46.3528 , 0.48364 , -241.239 , 0 , -241.239]; %[1]
    if result == "c_p"
        solution = c_pcalc(temps, schomateCoefficients);
    elseif result == "deltaH"
        solution = deltaH(temps, schomateCoefficients);
    end
end

function C_p = c_pcalc(temp, coefficents)
    C_p = coefficents(1) + coefficents(2) * temp + coefficents(3) * temp ^ 2 + coefficents(4) * temp ^ 3 + coefficents(5) / (temp ^ 2); 
end

function dH = deltaH(Temps, coefficents)
    startTemp = Temps(1);
    endTemp = Temps(2);
   
          
    startingH = schomateInt(startTemp, coefficents);
    endH = schomateInt(endTemp, coefficents);
    
    dH = endH - startingH;
end

function h = schomateInt(temp, coefficients)
    temp = temp / 1000;
    h = coefficients(1) * temp + 1/2 * coefficients(2) * temp ^ 2 + 1/3 * coefficients(3) * temp ^ 3 + 1/4 * coefficients(4) * temp ^ 4 - coefficients(5) / temp + coefficients(6) - coefficients(7);
end

%Sources
%{
    [1] Ethanol Schomate Coefficients (<> ADD NAME <>): https://strathprints.strath.ac.uk/6704/6/strathprints006704.pdf
%}