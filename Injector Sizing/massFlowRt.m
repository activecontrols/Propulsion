function [mDot] = massFlowRt(area, velocity, density)

% equation for mass flow rate 
mDot = density * area * velocity;