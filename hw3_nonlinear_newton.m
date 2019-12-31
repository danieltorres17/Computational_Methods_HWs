%%
% Daniel Torres
% Homework 3 Problem 3
% CSCI 3656 
% Fall 2019

close all; clear all; clc;

%%
% define functions with two inputs
f = @(u,v) [u^2+v^2-1; u^3-v^3+u];
df = @(u,v) [2*u, 2*v; 3*u^2+1, -3*v^2]; % Jacobian

% define functions with one vector input
f_vec = @(x) [x(1)^2+x(2)^2-1; x(1)^3-x(2)^3+x(1)];
df_vec = @(x) [2*x(1), 2*x(2); 3*x(1)^2+1, -3*x(2)^2]; % Jacobian

%% 
% Part 1: Plot for f1 = 0 and f2 = 0

% Used this to implement fsolve: https://www.mathworks.com/help/optim/ug/fsolve.html
x1 = [-1;-1];
rts1 = fsolve(f_vec,x1); % check f1 and f2 intersection
x2 = [1;1];
rts2 = fsolve(f_vec,x2); % check f1 and f2 intersection

u1 = -1:0.01:-.05; % first interval to evaluate f1 over
u2 = 0.05:0.01:1; % second interval to evaluate f1 over
v1 = -1:0.01:1; % interval to evaluate f2 over

f1 = @(x) abs((x.^3+x).^(1/3)); % f1 function
f2 = @(x) sqrt((1-x.^2)); % f2 function 

a1 = f1(u1); % evaluate f1 in first interval 
a2 = f1(u2); % evaluate f1 in second interval
b = f2(v1); % evaluate f2

plot(u1,-a1,'k','LineWidth',2); % plot f1
hold on
plot(u2,a2,'k','LineWidth',2); % plot f1
hold on
plot(v1,b,'r','LineWidth',2); % plot f2 (top half)
hold on
plot(v1,-b,'r','LineWidth',2); % plot f2 (bottom half)
hold on
plot(rts1(1),rts1(2),'o','LineWidth',3); % plot roots estimated with fzero
hold on
plot(rts2(1),rts2(2),'o','LineWidth',3); % plot roots estimated with fzero
axis square
grid on
xlim([-1.1 1.1]); % set x axis limits
ylim([-1.1 1.1]); % set y axis limits
title('u vs v');
xlabel('u');
ylabel('v');
legend('f1','Location','SouthEast');
legend boxoff

%% 
% Part 3: run Newton with several initial guesses and make table

x0 = [1;2]; % initial guesses for Newton's method
r = my_2d_newton(x0,f_vec,df_vec); % call Newton's method function

% try several starting points
n = 15; % number of random points to evaluate
init = zeros(n,2); % preallocate space for initial guess vector
roots = zeros(n,2); % preallocate space for estimated root vector
for i=1:n
   % create random starting point
   x0 = randn(1,2); % get two random normalized guesses
   r2 = my_2d_newton(x0',f_vec,df_vec); % call Newton method function
   init(i,:) = x0'; % put initial guesses into a vector
   roots(i,:) = r2; % put estimated roots into a vector
end

% All the estimated root pairs that were calculated with Newton's method are
% exactly what I calculated with fsolve and both roots are either positive or
% negative which makes sense since this system has two solutions.

% create table with four individual columns
T = table(init(:,1),init(:,2),roots(:,1),roots(:,2));
T.Properties.VariableNames = {'Init_1', 'Init_2', 'Newton_1' 'Newton_2'};
disp(T);

%%
% Part 4: Finding a point where Newton's method fails

x1 = [0;0];
r1 = my_2d_newton(x1,f_vec,df_vec);

% I chose to use [0;0] as the initial guesses because this is the only
% place that I found that Newton's method didn't work. I tried using other
% starting points very far out like [100;100] but that also converged. 
%The reason it doesn't converge at [0;0] is because the derivative of f2 = 0; 

%%
% This function was implemented with the help of Dr. Constantine in the
% Matlab review session he held on Thursday September 26th
function [r] = my_2d_newton(x0,f,df)
% INPUTS
% x0: initial guess, column vector
% f: system of functions
% df: Jacobian

% OUTPUTS
% r = estimated roots

    TOL = 1e-9; % error tolerance
    err = 1e9;
    iter_count = 0; % iteration counter initialization
    
    while err > TOL 
        % evaluate Jacobian and right hand side
        J = df(x0);
        R = -f(x0); 

        % solve linear system for Newton step
        p = J \ R; % DONT DO THIS -> p = inv(J)*R;

        % Newton step
        x1 = x0+p;
        
        % calc error
        err = norm(x0-x1)/norm(x0);

        % reset variables
        x0 = x1;
        iter_count = iter_count+1;
        
        % max iterations is 100
        if iter_count > 100 
            break;
        end
    end
    
    r = x0;
end