function [PropMass, FlightTime, ThrustDev] = TOAD_3DoF_SIM(TWR, plotmode)
addpath("lib\")
%Variable set-up
syms r_1 r_2 v_1 v_2 theta theta_dot m m_dot thrust alpha_ang
%Constants
syms g init_m radius lengthTOAD max_mdot max_thrust
%Time variable 
syms t 

%Declare constants 
consts = transpose([g init_m radius lengthTOAD max_mdot max_thrust]);
r = [r_1; r_2];
v = [v_1; v_2];
dictionary = [r,     v;
              theta, theta_dot;
              m,    m_dot];

%Set up assumptions for symbolic variables
assume(imag(dictionary) == 0);
assume(imag(consts) == 0);
assume(consts > 0);
assume(thrust > 0);
assume(alpha_ang > 0);
assume(m > 0);

%State Vector
x = [r; v; theta; theta_dot; m];

%Input vector
u = [thrust; alpha_ang];

%BF Thrust
BF = thrust*[sin(alpha_ang); cos(alpha_ang)];

%Rotation matrix from Inertial frame to Body frame
I2B = [cos(theta), -sin(theta);
       sin(theta),  cos(theta)];

%Rotation matrix from Body frame to Inertial frame
B2I = [cos(-theta), -sin(-theta);
       sin(-theta),  cos(-theta)];

%Inertial frame net force
IF = B2I*BF - [0; m*g];

%Derivative of position
r_dot = v;

%Derivative of velocity
v_dot = IF/m;

%Derivative of mass
m_dot = -thrust / max_thrust * max_mdot;

%Torque
tau = lengthTOAD*thrust*sin(alpha_ang);

%Moment of inertia
I = (1/3)*m*radius^2+(1/4)*m*2*lengthTOAD^2;

%Derivative of angular velocity
theta_ddot = tau/I;

%Derivative of state vector
x_dot = [r_dot; v_dot; theta_dot; theta_ddot; m_dot];

%Create numerical functions
matlabFunction(x_dot, 'File', './lib/odefcn2.m', 'Vars', {t, x, u, consts});
%%
%Define constants and initial conditions

g = 9.81;                           %m/s^2
thrust_max = 2446.52;               %Newtons
mass_i = thrust_max / (TWR * g);    %kg
rad = 0.0254;                       %m
len = 4;                            %m
mdot_max = 1.234;                   %kg/s

constants = [g; mass_i; rad; len; mdot_max; thrust_max];
x_0 = [0; 0; 0; 0; 0; 0; mass_i];

%Calculate jacobians for X and U
J_x = jacobian(x_dot, x);
J_u = jacobian(x_dot, u);

%System critical points
x_crit = [0;0;0;0;0;0;mass_i];
u_crit = [constants(2)*constants(1);0];

%Generate A and B matrices using jacobians, constants, and critical points
A = double(subs(J_x, [x; u; consts], [x_crit; u_crit; constants]));
B = double(subs(J_u, [x; u; consts], [x_crit; u_crit; constants]));

%C and D matrices
C = eye(size(x, 1));
D = 0;

%Set up LQR controller
a_weights = [5; 8; 1; 5; 2; 2; 0];
b_weights = [1; 1];
maxA = 0.5*g;

a_weights = a_weights / norm(a_weights);
b_weights = b_weights / norm(b_weights);

max_x = [10, 150, 1, 4, pi/6, pi/6, 1000];
max_u = [(maxA / g + 1) * mass_i* 0.5 * g , pi/24];

Q = eye(size(x,1)) .* a_weights ./ max_x.^2;
R = eye(size(u,1)) .* b_weights ./ max_u.^2;

% Gain Matrix
[K, ~, ~] = lqr(A, B, Q, R);

%%
% Simulation
clear ref_generator;
clear inputfcn;
tspan = [0 45];

% Solve using ODE45 and nonlinear dynamics
[tsim, xsim] = ode45(@(tsim, xsim) odefcn2(tsim, xsim, inputfcn(K, xsim, tsim), constants), tspan, x_0);

% Solve for inputs
u = zeros(size(xsim,1 ), 2);
for i = 1:1:size(xsim, 1)
    x_i = xsim(i,1:7).';
    t_i = tsim(i);
    u(i,1:2) = inputfcn(K, x_i, t_i);
end
%% OUTPUTS
EndIndex = find(xsim(100:end,2) < 0.1 & abs(xsim(100:end,4)) < 10);
MaxAccel = 0.5*g;
if isempty(EndIndex)
    PropMass = 0;
    FlightTime = 0;
    ThrustDev = 0;
else
    EndIndex = EndIndex(1) + 100;
    PropMass = mass_i - xsim(EndIndex, 7);
    FlightTime = tsim(EndIndex);
    Accel = (xsim(2:end,2) - xsim(1:end-1,2)) ./ (tsim(2:end) - tsim(1:end-1));

    ThrustDev = u(:,1) - 0.7 * thrust_max;
    ThrustDev = sum(ThrustDev.^2);
end



%%
% Plots
if plotmode == 1

    Cscale = 1.5;
    
    figure(1);
    subplot(2,2,1)
    Color1 = [0.9290 0.6940 0.1250] + (TWR - 1.6) / Cscale;
    Color1 = max([0 0 0], Color1);
    Color1 = min([1 1 1], Color1);
    
    Color2 = [0.8500 0.3250 0.0980] + (TWR - 1.6) / Cscale;
    Color2 = max([0 0 0], Color2);
    Color2 = min([1 1 1], Color2);
    
    plot(tsim, xsim(:,1), 'Color',Color1, 'LineWidth',1);
    hold on; grid on;
    plot(tsim, xsim(:,3), 'Color',Color2, 'LineWidth',1);
    title('Lateral Velocity and Position');
    legend('Positon', 'Velocity');
    xlabel('Time [s]');
    ylabel('Velocity and Position [m/s] [m]');
    
    subplot(2,2,2)
    Color1 = [0 0.4470 0.7410] + (TWR - 1.6) / Cscale;
    Color1 = max([0 0 0], Color1);
    Color1 = min([1 1 1], Color1);
    
    Color2 = [0.4660 0.6740 0.1880] + (TWR - 1.6) / Cscale;
    Color2 = max([0 0 0], Color2);
    Color2 = min([1 1 1], Color2);
    
    plot(tsim, xsim(:,2), 'Color',Color1, 'LineWidth',1);
    hold on; grid on;
    plot(tsim, xsim(:,4), 'Color',Color2, 'LineWidth',1);
    title('Vertical Velocity and Position');
    legend('Positon', 'Velocity');
    xlabel('Time [s]');
    ylabel('Velocity and Position [m/s] [m]');
    
    subplot(2,2,3)
    
    Color1 = [0.3010 0.7450 0.9330] + (TWR - 1.6) / Cscale;
    Color1 = max([0 0 0], Color1);
    Color1 = min([1 1 1], Color1);
    
    Color2 = [0.4940 0.1840 0.5560] + (TWR - 1.6) / Cscale;
    Color2 = max([0 0 0], Color2);
    Color2 = min([1 1 1], Color2);
    
    plot(tsim, xsim(:,5)*180/pi, 'Color',Color1, 'LineWidth',1);
    hold on; grid on;
    plot(tsim, xsim(:,6)*180/pi, 'Color',Color2, 'LineWidth',1);
    title('Angular Velocity and Position');
    legend('Positon', 'Velocity');
    xlabel('Time [s]');
    ylabel('Angular Velocity and Position [deg/s] [deg]');
    
    subplot(2,2,4)
    
    Color1 = [0.6350 0.0780 0.1840] + (TWR - 1.6) / Cscale;
    Color1 = max([0 0 0], Color1);
    Color1 = min([1 1 1], Color1);
    
    plot(tsim, u(:,1) / thrust_max * 100, 'Color', Color1, 'LineWidth',1);
    hold on; grid on;
    title('TOAD Throttle');
    xlabel('Time [s]');
    ylabel('Throttle [%]');
    ylim([0 100]);
    yline(40, 'r--');
    
    sgtitle('TOAD Hop 3DoF Simulation');

end