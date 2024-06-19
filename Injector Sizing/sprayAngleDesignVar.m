function [sprayAngles] = sprayAngleDesignVar(pintleTipAngles, mDotFuel, Ufuel, mDotOx, Uox)

% calculate total mometum ratio
TMR = (mDotFuel .* Ufuel .* cosd(pintleTipAngles)) ./ ... 
       (mDotOx .* Uox + mDotFuel .* Ufuel .* sind(pintleTipAngles));

% spray angle formula in degrees
sprayAngles = acosd(1 ./ (1 + TMR));