%% Pendulum Model with Controller (with noise/disturbances added %%

% This function takes in the desired angle vs. time and returns the
% "measured" angle, rate, acceleration, and motor commands to attempt to
% follow the desired trajectory
function [y,yd, ydd, u_cmd, t] = pend(ydes, dt, tfinal,control_on) 

m = 0.75; %kg
g = 9.81; %m/s^2
l = 0.33; %m
C = -0.064; %1/s
D = -0.011; %unitless
T = 0.2; %1/(kg*m^2)

t = 0:dt:tfinal; % Create time vector

N = length(t); % Total number of time steps/length of simulation vector

% Initialize output vectors
y = zeros(N,1); 
yd = zeros(N,1);
ydd = zeros(N,1);
u = zeros(N,1);
u_cmd = zeros(N,1); % Simply a converted input vector to convert to Voltage (to look more like real data)
error = zeros(N,1); % Error between desired and measured trajectory
error_dot = zeros(N,1);
error_integ = zeros(N,1);

if control_on == 1 % 1 = PID control active, 0 = no controller

    % PID Controller Gains
    Kp = 150.0;
    Kd = 10;
    Ki = 30.;


    % Run the simulation (through each time step), start at 2nd time step
    % (initialize with 0s for all initial states
    for i = 2:N
        % Update motor input from PID controller
        error(i) = ydes(i-1) - y(i-1); % Error components need for PID controller
        error_dot(i) = (error(i) - error(i-1))/dt;
        error_integ(i) = error_integ(i-1) + error(i)*dt; 
        u(i) = Kp*error(i) + Kd*error_dot(i) + Ki*error_integ(i);

        % Update pendulum dynamics/states
        ydd(i) = -g/l*2.0*sin(y(i-1)) + T/m*u(i)*abs(u(i)) + C*yd(i-1) + D*yd(i-1)*abs(yd(i-1)) + (rand()-0.5)*30;
        yd(i) = yd(i-1) + ydd(i)*dt + (rand()-0.5)*0.1;
        y(i) = y(i-1) + yd(i)*dt + ydd(i)*dt^2*0.5;
        
    end
else % No control, ydes is the motor command (not desired angle)
    u = ydes;
    for i = 2:N
        % Update pendulum dynamics/states
        ydd(i) = -g/l*2.0*sin(y(i-1)) + T/m*u(i)*abs(u(i)) + C*yd(i-1) + D*yd(i-1)*abs(yd(i-1)) + (rand()-0.5)*30;
        yd(i) = yd(i-1) + ydd(i)*dt + (rand()-0.5)*0.1;
        y(i) = y(i-1) + yd(i)*dt + ydd(i)*dt^2*0.5;
    end
end

u_cmd = u;









