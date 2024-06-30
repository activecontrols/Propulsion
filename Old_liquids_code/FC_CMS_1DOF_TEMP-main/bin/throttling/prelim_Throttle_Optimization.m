% Tyler Trostle
% 11/06/2021

conversions;

%%initial rocket 3 parameters
m_dry_lbm = 150;
Cd = 0.4;
R_OD = 8;
Isp = 215;
impulse_lbfs = 30000; %9208 lbf*s or 24000 lbf*s ?

%%calculated parameters
m_prop_0 = (impulse_lbfs * c.LBF2N / (Isp * 9.8)) * c.KG2LBM; %prop mass lbm
init_mass = m_dry_lbm + m_prop_0; %rocket mass (lbm)

%%iterated parameters
%sets min thrust to twice the rocket mass rounded up to the nearest multiple of 25
%minThrust = 2 * 25 * ceil(init_mass / 25); 
thrust_lbf = 2500;
minThrottle = .7:.1:1;
minThrottleTime = 0:.05:1;


maxAltMat = zeros(size(minThrottleTime,2),size(minThrottle,2)); %initialize maximum Altitude array

for j = 1:size(minThrottle,2)
    for i = 1:size(minThrottleTime,2)

        %call trajectory model and save max altitude to array
        [maxAltMat(i,j), ~, ~, ~] = Traj_1DoF_Model_Throttling(m_dry_lbm, Cd, R_OD, thrust_lbf, Isp, impulse_lbfs, 0, minThrottle(j), minThrottleTime(i));
    end
end

%maxAlt = max(maxAltMat); %highest altitude achieved (ft)
%throttleMaxAlt = minThrottle(maxAltMat == max(maxAltMat)); %thrust of highest altitude (lbf)
%throttleTimeMaxAlt = minThrottleTime(maxAltMat == max(maxAltMat)); %initial TWR of highest altitude

%fprintf("\nMaximum Altitude of %0.0f ft with Max Throttle = %0.1f lbf and Throttle Time = %0.2f\n\n", maxAlt, throttleMaxAlt, throttleTimeMaxAlt);

figure(1);

plot(minThrottleTime, maxAltMat(:,1), 'LineWidth',2);
hold on
plot(minThrottleTime, maxAltMat(:,2), 'LineWidth',2);
plot(minThrottleTime, maxAltMat(:,3), 'LineWidth',2);
plot(minThrottleTime, maxAltMat(:,4), 'LineWidth',2);

legend("Min Throttle: .7", "Min Throttle: .8", "Min Throttle: .9", "Min Throttle: 1")
grid on
title("Maximum Altitude (ft) vs. Throttling Amount and Time with Impulse = " + impulse_lbfs + " lbf*s");
xlabel("Point of Minimum Thrust (fraction of burn time)");
ylabel("Maximum Altitude (ft)");




