
burnTime = 100;
minThrottle = .7;
minThrottleTime = .8;

for t = 1:1:burnTime
    throttle1(t) = thrustCurve1(t, burnTime, minThrottle);
    %throttle2(t) = thrustCurve2(t, burnTime, minThrottle, minThrottleTime);
end

time = 1:1:100;

figure(1)
plot(time, throttle1)
hold on
%plot(time, throttle2)
axis([0 burnTime 0 1]);
xlabel("Elapsed Time (% of burn)");
ylabel("Throttle Amount");
title("Thrust vs Time");