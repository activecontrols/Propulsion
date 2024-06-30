% Tyler Trostle
% 11/06/2021

conversions;

%%initial rocket 3 parameters
mDry = 150;
Cd = 0.4;
tubeOD = 8.125;
burnoutTime = 20;
isp = 215;

%%calculated parameters
mProp = (impulse_lbfs * c.LBF2N / (isp * 9.8)) * c.KG2LBM; %prop mass lbm
mWet = mDry + mProp; %rocket mass (lbm)

%%iterated parameters
%sets min thrust to twice the rocket mass rounded up to the nearest multiple of 25
minThrust = 2 * 25 * ceil(mWet / 25); 
thrust = minThrust:25:3000; %thrust array

maxAltMat = zeros(1,size(thrust,2)); %initialize maximum Altitude array

for i = 1:size(thrust,2) %iterate over thrust

    %call trajectory model and save max altitude to array
    [maxAltMat(i), ~, ~, ~] = Traj_1DoF_Model(mDry, Cd, tubeOD, thrust(i), burnoutTime, isp, 0);
    
end

maxAlt = max(maxAltMat); %highest altitude achieved (ft)
thrustMaxAlt = thrust(maxAltMat == max(maxAltMat)); %thrust of highest altitude (lbf)
TWRMaxAlt = thrustMaxAlt / mWet; %initial TWR of highest altitude

fprintf("\nMaximum Altitude of %0.0f ft with Thrust = %0.0f lbf and initial TWR = %0.2f\n\n", maxAlt, thrustMaxAlt, TWRMaxAlt);

figure(1);
plot(thrust, maxAltMat,'LineWidth',4);
grid on
title("Maximum Altitude (ft) vs. Thrust (lbf) with Impulse = " + impulse_lbfs + " lbf*s");
xlabel("Thrust (lbf)");
ylabel("Maximum Altitude (ft)");

figure(2);
plot(thrust / mWet, maxAltMat,'LineWidth',4);
grid on
title("Maximum Altitude (ft) vs. Initial TWR with Impulse = " + impulse_lbfs + " lbf*s");
xlabel("Thrust to Weight Ratio at t = 0");
ylabel("Maximum Altitude (ft)");