function throttle = thrustCurve2(time, burnTime, minThrottle, minTime)


if time < (minTime * burnTime)
    
    throttle = 1 - (1 - minThrottle) * (time / (minTime*burnTime));
    
elseif time <= burnTime
    
    throttle = minThrottle + (1 - minThrottle) * ((time - burnTime*minTime) / (burnTime - burnTime*minTime));
    
else
    throttle = 0;
    
end

