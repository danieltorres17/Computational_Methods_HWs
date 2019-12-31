% Daniel Torres
% Homework 8 
% CSCI 3656
% Fall 2019

clear vars; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 1 - Solve first order ODE using Forward Euler and RK4

k = 1:15; % integer vector to determine step sizes
h = 2.^-k; % step sizes
t1 = [0 8]; % time interval 
y0 = 0; % initial condition
t_n = (t1(2)-t1(1))./h; % number of steps vector

yp = @(t,y) -y+sin(t); % ODE 
y1 = @(t) 0.5*(exp(-t)+sin(t)-cos(t)); % analytical solution to yp

y1_true = y1(t1(2)); % true solution value at 8 seconds

% Part A - Verify Solution
% Check the write up for the verification

% Part B - Solve ODE using Forward Euler
rel_err_fwdE = zeros(1,length(h));
for i = 1:length(h)
    % approximate Forward Euler solution
    [t_fwdE,y_fwdE] = fwd_euler(yp,y0,h(i),t_n(i));
    % calculate relative error
    rel_err_fwdE(i) = abs((y1_true-y_fwdE(end)))/y1_true;
end
figure(1);
loglog(h,rel_err_fwdE,'-o','LineWidth',2); % plot Forward Euler results
hold on; 
grid on;
axis square;

% The asymptotic regime for this method covers the range of h's we chose.
% The observed convergence rate is O(h). I plotted a line to show that this
% has a first order convergence rate. 

% Part C - Solve ODE using Runge-Kutta 4
rel_err_rk4 = zeros(1,length(h));
for i = 1:length(h)
     % approximate RK4 solution
    [t_rk4,y_rk4] = rk4(yp,y0,h(i),t_n(i));
    % calculate relative error
    rel_err_rk4(i) = abs((y1_true-y_rk4(end)))/y1_true;
end
loglog(h,rel_err_rk4,'-o','LineWidth',2); % plot RK4 results

% The asymptotic regime for this method is from 10^-3 to about 10^(-.5).
% The observed convergence rate is O(h^4). I plotted a line to show that this
% has a fourth order convergence rate. 

% Estimate convergence rate
loglog(h,h,'LineWidth',2); % for Forward Euler convergence rate comparison
loglog(h,h.^4,'LineWidth',2); % for RK4 convergence rate comparison 

title('log(h) vs log(Relative Error)');
xlabel('h');
ylabel('Error');
legend('Fwd Euler','RK4','1st Order','4th Order','Location','SouthEast');
legend boxoff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 - Solve 3rd order ODE using ode45

u0 = [-3 -2 2]; % initial conditions
t2 = [0 5]; % time span [s]
y2 = @(t) -sin(2*t)+t^2-3; % analytical solution to ODE
y2_true = y2(t2(2)); % true solution value at 5 seconds

% Part A - Verify Solution
% Check notes included in the write up for the verification

% Part B - Rewrite Problem as First Order System
% Check notes included in the write up 

% Part C - Using ode45 to approximate solution y(t)
[t2,y2] = ode45(@(t,y) func(t,y),t2,u0); % solve using ode45 
figure(2);
plot(t2,y2(:,1),'LineWidth',2); % plot results
grid on; axis square;
title('ODE Solution Using ode45');
xlabel('t');
ylabel('v(t)');

% Part D - Tolerance Parameter Experiment
k2 = 1:10; % integer vector to determine tolerances
tol = 10.^-k2; % tolerance vector
rel_err2 = zeros(1,length(tol)); % memory allocation

for i = 1:length(tol)
    % approximate ode45 solution
    [t3,y3] = ode45(@(t,y) func(t,y),t2,u0,odeset('RelTol',tol(i)));
    % isolate y(t) approx solution vector
    y_val = y3(:,1);
    % calculate relative error 
    rel_err2(i) = abs((y2_true-y_val(end)))/y2_true;
end

% create table with tolerance and relative error values
T = table(tol',rel_err2'); 
T.Properties.VariableNames = {'Tolerance','Norm_of_rel_error'};
disp(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ODE Solver Functions 

function [t,y] = fwd_euler(f,y0,h,n)
    % This method was implemented using the formula given by Professor 
    % Constantine in class

    % INPUTS: 
    % f: function to evaluate
    % y0: initial condition
    % h: step size
    % n: number of steps 
    
    % OUTPUTS:
    % t: time vector
    % y: approximate solution vector
    
    t = zeros(1,n);  % time vector memory allocation
    y = zeros(1,n);  % approximated solution vector memory allocation
    y(1) = y0; % place initial condition as first value in solution vector
    
    for i = 1:n
        t(i+1) = t(i)+h;
        y(i+1) = y(i)+h*f(t(i),y(i));
    end
    
end

function [t,y] = rk4(f,y0,h,n)
    % This method was implemented using the formula found here:
    % https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#

    % INPUTS: 
    % f: function to evaluate
    % y0: initial condition
    % h: step size
    % n: number of steps 
    
    % OUTPUTS:
    % t: time vector
    % y: approximate solution vector

    t = zeros(1,n); % time vector memory allocation
    y = zeros(1,n); % approximated solution vector memory allocation
    y(1) = y0; % place initial condition as first value in solution vector
    
    for i = 1:n
        t(i+1) = t(i)+h;
        k1 = h*f(t(i),y(i));
        k2 = h*f(t(i)+(h/2),y(i)+(k1/2));
        k3 = h*f(t(i)+(h/2),y(i)+(k2/2));
        k4 = h*f(t(i)+h,y(i)+k3);
        y(i+1) = y(i)+(1/6)*(k1+2*k2+2*k3+k4); 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for ode45 solver

function [dx] = func(t,y)
    % check the notes included for how I broke down the system from a 
    % third-order to a first-order system

    dx(1) = y(2); % u2(t)
    dx(2) = y(3); % u3(t)
    dx(3) = 4*t^2+8*t-10-y(3)-4*y(2)-4*y(1);
    
    dx = dx';

end


