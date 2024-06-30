% Tyler Trostle
% Last Edit: 1/05/2022 (Tyler)

conversions;

%%initial rocket 3 parameters
m_dry_lbm = 3000; %[lbm]
Cd = 0.4; 
R_OD = 8.5; %[in]
thrust = 2500; %[lbf]
t_b = 15; %[s]
Isp = 215; %[s]

%%calculated parameters
m_prop_0 = (thrust * c.LBF2N * t_b / (Isp * 9.8)) * c.KG2LBM; %prop mass lbm
init_mass = m_dry_lbm + m_prop_0; %rocket mass (lbm)

%%iterated parameters
minThrottle = .6:.01:1;
%minThrottleTime = 0:.05:1;


maxAltMat = zeros(size(minThrottle,2),1); %initialize maximum Altitude array

for j = 1:size(minThrottle,2)

    %call trajectory model and save max altitude to array
    [maxAltMat(j), ~, ~, ~] = Traj_1DoF_Model_Throttling_V2(m_dry_lbm, Cd, R_OD, thrust, t_b, Isp, 0, minThrottle(j));
end

%maxAlt = max(maxAltMat); %highest altitude achieved (ft)
%throttleMaxAlt = minThrottle(maxAltMat == max(maxAltMat)); %thrust of highest altitude (lbf)
%throttleTimeMaxAlt = minThrottleTime(maxAltMat == max(maxAltMat)); %initial TWR of highest altitude

%fprintf("\nMaximum Altitude of %0.0f ft with Max Throttle = %0.1f lbf and Throttle Time = %0.2f\n\n", maxAlt, throttleMaxAlt, throttleTimeMaxAlt);

figure(1);

plot(minThrottle, maxAltMat, 'LineWidth',2);

grid on
title("Maximum Altitude (ft) vs. Throttling Amount");
xlabel("Minimum Throttle Ammount");
ylabel("Maximum Altitude (ft)");


figure(2)
plot(1:100, thrustCurve1(1:100, 100, .8))
axis([0 100 0 1]);
xlabel("Elapsed Time (% of burn)");
ylabel("Throttle Amount");
title("Thrust vs Time");


