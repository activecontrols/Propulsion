function throttle = thrustCurve1(time, burnTime, minThrottle)

%linear throttle down to min throttle across whole burn
throttle = 1 - (1 - minThrottle) * (time / burnTime);