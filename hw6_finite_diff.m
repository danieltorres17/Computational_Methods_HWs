% Daniel Torres - 105600180
% HW 6 
% CSCI 3656
% Fall 2019

clear all; close all; clc;

k = 5:24; % exponent vector
h = 2.^-k; % h values for finite diff estimation
x0 = 0.2; % value at which to approximate derivative at

y = @(x) sin(4.8*pi.*x); % f(x)
y_prime = @(x) 4.8*pi*cos(4.8*pi*x); % f'(x)

%% Part 1: Implement forward, central and backward difference 

yp_forward = forward_(y,x0,h); % estimate derivative using forward diff
yp_central = central_(y,x0,h); % estimate derivative using central diff
yp_backward = backward_(y,x0,h); % estimate derivative using backward diff
yp_fin_diff = fin_diff(y,x0,h); % estimate derivative using formula from part 2

true_value = y_prime(x0); % true derivative value using f'(x)


% For the following calculations of the convergence rates I used points
% from the linear sections of each of the plots. 

% forward diff estimated rate of convergence
% Used points 8 and 20 from h vector and yp_forward vector
forward_conv = (yp_forward(8)-yp_forward(20))/(h(20)-h(8));
fprintf('Forward difference rate of convergence: %f.\n',forward_conv);

% central diff estimated rate of convergence
% Used points 1 and 17 h vector and yp_central vector
central_conv = (yp_central(1)-yp_central(17))/(h(1)-h(17));
fprintf('Central difference rate of convergence: %f.\n',central_conv);

% backward diff estimated rate of convergence
% Used points 1 and 20 h vector and yp_backward vector
backward_conv = (yp_backward(1)-yp_backward(20))/(h(1)-h(20));
fprintf('Backward difference rate of convergence: %f.\n',backward_conv);

% part 2 diff estimated rate of convergence
% Used points 1 and 11 h vector and yp_fin_diff vector
part2_conv = (yp_fin_diff(11)-yp_fin_diff(1))/(h(11)-h(1));
fprintf('Part 2 difference rate of convergence: %f.\n',part2_conv);

% I wasn't able to figure out why my plot for part 1 didn't match more
% closely to the plot provided by Professor Constantine. My plot looks
% similar but it seems to be missing some points for the forward difference
% estimation. 

% plot results for part 1
figure(1);
loglog(h,abs(yp_forward-true_value),'LineWidth',2); % plot forward diff
hold on
loglog(h,abs(yp_central-true_value),'LineWidth',2); % plot central diff
hold on
loglog(h,abs(yp_backward-true_value),'LineWidth',2); % plot backward diff
title('Finite Difference Method Comparisons');
xlabel('h');
ylabel('Error');
legend('Forward','Central','Backward','Location','SouthEast');
legend boxoff
xlim([10^-8 10^0]);
grid on
axis square 

% plot results for part 2
figure(2);
loglog(h,abs(yp_fin_diff-true_value),'LineWidth',2); % plot part 2 diff
title('Part 2: Finite Difference Method');
xlabel('h');
ylabel('Error');
legend('Part 2 Formula','Location','SouthEast');
legend boxoff
xlim([10^-8 10^0]);
grid on
axis square 

%% Derivative Functions

function [yp] = forward_(fun,x0,h) % forward difference function
    % Inputs:
    % fun: anonynomous function to evaluate
    % x: point at which to estimate derivative
    % h: step size
    
    % Outputs:
    % yp: estimated derivative at x0
    
    yp = (fun(x0+h)-fun(x0))./h;

end

function [yp] = backward_(fun,x0,h) % backward difference function
    % Inputs:
    % fun: anonynomous function to evaluate
    % x: point at which to estimate derivative
    % h: step size
    
    % Outputs:
    % yp: estimated derivative at x0

    yp = (fun(x0)-fun(x0-h))./h;
    
end

function [yp] = central_(fun,x0,h) % central difference function 
    % Inputs:
    % fun: anonynomous function to evaluate
    % x: point at which to estimate derivative
    % h: step size
    
    % Outputs:
    % yp: estimated derivative at x0

    yp = (fun(x0+h)-fun(x0-h))./(2*h);

end

function [yp] = fin_diff(fun,x0,h) % part 2 finite difference formula
    % Inputs:
    % fun: anonynomous function to evaluate
    % x: point at which to estimate derivative
    % h: step size
    
    % Outputs:
    % yp: estimated derivative at x0

    yp = (1./(6*h)).*(2*fun(x0+h)+3*fun(x0)-6*fun(x0-h)+fun(x0-(2.*h)));

end
