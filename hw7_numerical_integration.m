% Daniel Torres
% HW 7
% CSCI 3656
% Fall 2019

clear all; close all; clc;

interval = [1 1.5]; % interval over which to integrate
a = interval(1); % first point in interval
b = interval(2); % second point in interval
f = @(x) cos(3*pi*x); % function to numerically integrate
int_f = @(x) (1/(3*pi))*sin(3*pi*x); % actual integral 

% P1 Part A: Computing definite integral of f(x) on interval [1,1.5]

true_value = int_f(b)-int_f(a); % actual value of integral in domain

% P1 Part B: Estimating integral 

k = 1:10;
n = 2.^k+1; % number of intervals

y_trap = zeros(1,length(n)); % allocate space
y_simp = zeros(1,length(n)); % allocate space
y_cc = zeros(1,length(n)); % allocate space

for i=1:length(n)
    y_trap(i) = trap(f,a,b,n(i)); % integrate using trapezoidal method
    y_simp(i) = simpson(f,a,b,n(i)); % integrate using simpsons method
    [x_vals,weights] = fclencurt(n(i),a,b); % calculate weights 
    y_cc(i) = sum(f(x_vals).*weights); % integrate using Clenshaw
end

err_trap = abs(y_trap-true_value); % calculate trapezoidal method error
err_simp = abs(y_simp-true_value); % calculate Simpson method error
err_cc = abs(y_cc-true_value); % calculate Clenshaw-Curtis method error

figure(1);
loglog(n,err_trap,'LineWidth',2); % plot trapezoidal error
hold on 
loglog(n,err_simp,'LineWidth',2); % plot Simpsons error
loglog(n,err_cc,'LineWidth',2); % plot Clenshaw-Curtis error
loglog(n,n.^-1); % plotted to estimate rate of convergence
loglog(n,n.^-2); % plotted to estimate rate of convergence
xlim([1 10^4]);
axis square
grid on
title('Part 1: Error vs n');
xlabel('n');
ylabel('Error');
legend('Trapezoid','Simpson','CC','Linear','Quadratic','Location','Southwest');
legend boxoff

% P1 Part C: Identifying values of n in asymptotic region

% Using the lines I plotted in the last section it can be seen that the
% trapezoid method converges in a linear rate, Simpson's method converges
% at a quadratic rate and I wasn't able to determine the rate of
% convergence for the Clenshaw-Curtis method. 

% I was able to read the values of n that constitude the asymptotic regime
% from the plot. The trapezoidal method and Simpsons method both had the
% same values of n which were 3-1000. The Clenshaw-Curtis method only had
% n-values from 4-11. 

% Problem 2

y_trap2 = zeros(1,length(n)); % allocate space
y_simp2 = zeros(1,length(n)); % allocate space
y_cc2 = zeros(1,length(n)); % allocate space

for i=1:length(n)
    y_trap2(i) = trap(@sign,a,b,n(i)); % integrate using trapezoidal method
    y_simp2(i) = simpson(@sign,a,b,n(i)); % integrate using simpson method
    [x_vals2,weights2] = fclencurt(n(i),a,b); % calculate weights
    y_cc2(i) = sum(sign(x_vals2).*weights2); % integrate using clenshaw method
end

true_value2 = 0;

err_trap2 = abs(y_trap2-true_value2); % calculate trapezoidal method error
err_simp2 = abs(y_simp2-true_value2); % calculate Simpsons method error
err_cc2 = abs(y_cc2-true_value2); % calculate Clenshaw-Curtis method error

figure(2);
loglog(n,err_trap2,'LineWidth',2); % plot trapezoidal error
hold on
loglog(n,err_simp2,'LineWidth',2); % plot Simpsons error
loglog(n,err_cc2,'LineWidth',2); % plot Clenshaw-Curtis error
loglog(n,n.^-1); % plotted to estimate rate of convergence
loglog(n,n.^-2); % plotted to estimate rate of convergence
xlim([1 10^4]);
axis square
grid on
title('Part 2: Error vs n');
xlabel('n');
ylabel('Error');
legend('Trapezoid','Simpson','CC','Location','Southwest');
legend boxoff

% The observed rates of convergence for the sign function differ from the
% original function we used in the first section because now Simpson's
% method seems to converge at a linear rate. Not sure why the trapezoidal
% method and the Clenshaw-Curtis method don't look like they work correctly
% here but they did in the section above. 

% Numerical Integration Methods

function [y] = trap(f,a,b,n) % trapezoidal method implementation

    % INPUTS:
    % f: function to evaluate
    % a: first point in interval
    % b: last point in interval
    % n: number of points to use for integration
    
    % OUTPUT:
    % y: estimated integral
 
    h = (b-a)/n; % calculate distance of intervals
    x = a:h:b; % split up interval 
    y = h*(sum(f(x))+0.5*(f(a)+f(b))); % add up values and multiply by h

end

function [y] = simpson(f,a,b,n) % Simpson's method implementation

    % used this link for implementation:
    % https://www.mathworks.com/matlabcentral/fileexchange/28726-simpson-s-rule-integration
    % INPUTS:
    % f: function to evaluate
    % a: first point in interval
    % b: last point in interval
    % n: number of points to use for integration
    
    % OUTPUT:
    % y: estimated integral

    h = (b-a)/n;  % calculate distance of intervals
    x = a:h:b; % split up interval 
    y = h/3*(f(x(1))+2*sum(f(x(3:2:end-2)))+4*sum(f(x(2:2:end)))+f(x(end)));

end

function [y] = sign(x) % sign function implementation from homework
    
    x = x-0.2; 
    a = 0;
    
    if x > 0
        a = 1;
    elseif x == 0
        a = 0;
    elseif x < 0
        a = -1;
    end
    
    a = a+1;
    y = a;
end

