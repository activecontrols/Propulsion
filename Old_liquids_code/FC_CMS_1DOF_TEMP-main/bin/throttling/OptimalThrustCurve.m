% Tyler Trostle
% Last Edit: 1/06/2022 (Tyler)

% Creates array of all possible step function thrust curves
minThrottle = .3;
maxThrottle = 1.3;

%for offset = 0:.1:.3

% minThrottle = .7 + offset;
% maxThrottle = 1 + offset;
    
    
steps = 6;
n = 4;
curves = n ^ steps;
throttleMat = ones(curves, steps);
for i = 0:(curves-1)
    curve = str2num(dec2base(i,n));
    for j = steps:-1:1
        scale = floor(curve / (10^(j-1))) / (n-1);
        curve = curve - floor(curve / (10^(j-1))) * 10^(j-1);
        throttleVal = minThrottle + (maxThrottle-minThrottle) * scale;
        throttleMat(i+1,steps-j+1) = throttleVal;
    end
end

% Vehicle parameters for Traj code
m_dry = 330;
Cd = .4;
outerDiameter = 8.5;
thrust_lbf = 2500;
t_b = 20;
isp_ex = 215;
mdot = 12; %mdot total

% Uncomment and run to get optimal thrust curve

alt = zeros(1,curves);
parfor i = 1:curves
    alt(i) = Traj_1DoF_Model_Throttling_V3(m_dry, Cd, outerDiameter, thrust_lbf, t_b, isp_ex, mdot, 0, throttleMat(i,:));
    
    fprintf("%d / %d\n", i, curves);
end


maxAlt = 0;
for i = 1:curves
    alt_current = alt(i);
    if alt_current > maxAlt
        maxAlt = alt_current;
        maxCurve = i;
    end
end

%Graphs optimal thrust curve
%(with n = steps = 5, maxCurve = 2501)
curveNum = maxCurve;
figure(1);
xMax = 100;
x = 1:.1:xMax;
y = throttleMat(curveNum,ceil(steps*x/xMax));
plot(x/xMax,y)
axis([0 1 0 1.1*maxThrottle])
title("Optimal Thrust Curve");
xlabel("Time Elapsed [0-1]");
ylabel("Throttle Setting [0-1]");


fprintf("\nTrajectory with NO throttling");
altNoThrottle = Traj_1DoF_Model_Throttling_V3(m_dry, Cd, outerDiameter, thrust_lbf, t_b, isp_ex, mdot, 1, ones(1,steps));
fprintf("\nTrajectory with throttling");
altThrottle = Traj_1DoF_Model_Throttling_V3(m_dry, Cd, outerDiameter, thrust_lbf, t_b, isp_ex, mdot, 1, throttleMat(maxCurve,:));


fprintf("\nMinThrottle: %0.2f    MaxThrottle: %0.2f\n", minThrottle, maxThrottle);
fprintf("Alt without throttling: %6.0f ft\nAlt with throttling:    %6.0f ft\n", altNoThrottle, altThrottle);
fprintf("Altitude gained: %6.0f ft (%2.2f%%)\n", altThrottle - altNoThrottle, (altThrottle/altNoThrottle-1)*100);
%end
