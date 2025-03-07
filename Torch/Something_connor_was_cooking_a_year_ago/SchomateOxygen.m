function solution = SchomateOxygen(temps, result) %Temps in K (Scalar, not united by symunit)
    schomateCoefficientsIndexes = [ 100 , 700 ; 
                                    700 , 2000; 
                                    2000 , 6000];

                            %A             B           C           D           E           F           G       H
    schomateCoefficients = [ 31.32234 , -20.23531 , 57.86644  , -36.50624 , -0.007374 , -8.903471 , 246.7945 , 0.0
                             30.03235 , 8.772972  , -3.988133 , 0.788313  , -0.741599 , -11.32468 , 236.1663 , 0.0
                             20.91111 , 10.72071  , -2.020498 , 0.146449  , 9.245722  , 5.337651  , 237.6185 , 0.0];
    
    if result == "c_p"
        solution = c_pcalc(temps, schomateCoefficientsIndexes, schomateCoefficients);
    elseif result == "deltaH"
        solution = deltaH(temps, schomateCoefficientsIndexes, schomateCoefficients);
    end
end



function C_p = c_pcalc(temp, tempIndexes, coefficents)
    selectedRow = -1;
    for row = (1:height(tempIndexes))
        %fprintf("\n%d vs %d: %d in %d - %d\n",selectedRow, row, temp, tempIndexes(row, 1), tempIndexes(row, 2))
        if (temp <= tempIndexes(row, 2)) && (temp >= tempIndexes(row, 1))
            selectedRow = row;
        end
    end

    if selectedRow == -1
        fprintf("\nError C_p temp Out of range\n")
    end
    temp = temp / 1000;
    coefficents = coefficents(selectedRow,:);   
    C_p = coefficents(1) + coefficents(2) * temp + coefficents(3) * temp ^ 2 + coefficents(4) * temp ^ 3 + coefficents(5) / (temp ^ 2); 
end

function dH = deltaH(Temps,tempIndexes, coefficents)
    startTemp = Temps(1);
    endTemp = Temps(2);

    startingRow = rowSelect(startTemp, tempIndexes);
    if startingRow == -1
        fprintf("\nError deltaH startingtemp Out of range\n")
    end

    endingRow = rowSelect(endTemp, tempIndexes);
    if endingRow == -1
        fprintf("\nError deltaH startingtemp Out of range\n")
    end

    %fprintf("Rows %d-%d\n", startingRow, endingRow)
   
          
    startingH = schomateInt(startTemp, coefficents(startingRow,:));
    endH = schomateInt(endTemp, coefficents(endingRow,:));
    
    dH = endH - startingH;
end


function returnRow = rowSelect(value, table)
    returnRow = -1;
    for row = (1:height(table))
        %fprintf("Start Temp %d vs %d: %d in %d - %d\n",returnRow, row, value, table(row, 1), table(row, 2))
        v2 = table(row, 2);
        v1 = table(row, 1);
        if (value <= table(row, 2)) && (value >= table(row, 1))
            returnRow = row;
        end
    end
end

function h = schomateInt(temp, coefficients)
    temp = temp / 1000;
    h = 0;
    for row = (1:height(coefficients))
        h = h + coefficients(1) * temp + 1/2 * coefficients(2) * temp ^ 2 + 1/3 * coefficients(3) * temp ^ 3 + 1/4 * coefficients(4) * temp ^ 4 - coefficients(5) / temp + coefficients(6) - coefficients(7);
    end

end
