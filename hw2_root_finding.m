% Daniel Torres
% Homework 2 Problems 1,2,3
% CSCI 3656 - Fall 2019

clear all; close all; clc;

format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Implement bisection, Newton and secant methods for root
% finding

% check below for method implementations

% For the bisection method I used the textbook Numerical Analysis 2nd
% Edition by Sauer, pages 27 and 28

% For the Newton's method implementation I used both the textbook Numerical 
% Analysis 2nd Edition by Sauer and this website that I found:
% http://hplgit.github.io/Programming-for-Computations/pub/p4c/._p4c-bootstrap-Python027.html

% For the secant method implementation I used both the textbook Numerical 
% Analysis 2nd Edition by Sauer and this website that I found:
% http://hplgit.github.io/Programming-for-Computations/pub/p4c/._p4c-bootstrap-Python028.html#4th:NonlinAlgEq:Secant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem 2 - Compute roots of f(x) = 4x^2-3x-3 using the quadratic
% function 

% function f(x) coefficients
a = 4; 
b = -3;
c = -3;

r = quadratic(a,b,c); % true value of roots
root_1 = r(1);
root_2 = r(2);
fprintf('Root 1: %f\n',root_1);
fprintf('Root 2: %f\n',root_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem 3 - Plot the error of each method

f = @(x) 4*x^2-3*x-3; 
df = @(x) 8*x-3;

tol = 1e-5; % tolerance 

% bisection method results and plot
[r_b,B_values,B_iters] = bisection_(f,-1,0,tol);
fprintf('Bisection root found: %f\n',r_b);
bisec_errors = abs(B_values-root_2);
semilogy(bisec_errors,'LineWidth',2);

% Newton's method results and plot
[r_N,N_values,N_iters] = Newton_(f,df,-1,tol);
fprintf('Newton root found: %f\n',r_N);
newton_errors = abs(N_values-root_2);
hold on;
semilogy(newton_errors,'LineWidth',2);

% secant method results and plot
[r_s,s_values,s_iters] = secant_(f,-1,0,tol);
fprintf('Secant root found: %f\n',r_s);
secant_errors = abs(s_values-root_2);
hold on;
semilogy(secant_errors,'LineWidth',2);

% adding title, axes labels, and legend to plot
title('Error as function of iterations');
xlabel('Iterations');
ylabel('Error');
legend('Bisection','Newton','Secant');
legend('boxoff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,values,iters] = bisection_(f,a,b,tol) % Bisection method implementation
    %%%%%%%%%%%%%%%%%%%%
    % Inputs:
    % f: function handle
    % a: starting bound
    % b: ending bound
    % tol: tolerance 
    
    % Outputs:
    % r: estimated root
    % values: vector of function values at each iteration
    % iters: number of iterations
    %%%%%%%%%%%%%%%%%%%%
    
    % this algorithm was implemented using the code from Numerical Analysis
    % 2nd Edition by Sauer, pages 27 and 28 
    
    iter = 0; % initialize iteration count variable
    i = 0; % initialize counter for values vector
  
    if(sign(f(a)*f(b)) > 0) % check to see if root exists between a and b
       error('No root can be found'); 
    end
    f_a = f(a); % evaluate function at a
    
    while((b-a)/2 > tol) % check for convergence
       c = (a+b)/2;
       f_c = f(c);
       err(i+1) = c; % add value to err vector
       i = i+1;  % increase value vector counter size
       if(f_c == 0)
          break; 
       end
       if(sign(f_a*f_c) < 0) % adjust bounds
           b = c;
       else % adjust bounds
           a = c;
       end
       iter = iter + 1; % increase iteration count
    end
    
    r = (a+b)/2; % return estimated root
    values = err; % return values at every iteration
    iters = iter; % return number of iterations
end

function [r,values,iters] = Newton_(f,df,x,tol) % Newton's method implementation
    %%%%%%%%%%%%%%%%%%%%
    % Inputs:
    % f: function handle
    % df: derivative of function
    % x: initial guess
    % tol: tolerance 
    
    % Outputs:
    % r: estimated root
    % values: vector of function values at each iteration
    % iters: number of iterations
    %%%%%%%%%%%%%%%%%%%%
    
    % partly implemented by using this website:
    % http://hplgit.github.io/Programming-for-Computations/pub/p4c/._p4c-bootstrap-Python027.html
    
    iter = 0; % initialize iteration count variable
    i = 0; % initialize counter for values vector
    
    while(abs(f(x)) > tol) % check for convergence
       x = x - f(x)/df(x); % iteration 
       err(i+1) = x; % add value to err vector
       i = i + 1; % increase value vector counter size
       iter = iter + 1; % increase iteration count
    end
    r = x; % return estimated root
    values = err; % return values at every iteration
    iters = iter; % return number of iterations
end

function [r,values,iters] = secant_(f,x0,x1,tol) % secant method implementation
    %%%%%%%%%%%%%%%%%%%%
    % Inputs:
    % f: function handle
    % x0: first initial guess
    % x1: second initial guess
    % tol: tolerance 
    
    % Outputs:
    % r: estimated root
    % values: vector of function values at each iteration
    % iters: number of iterations
    %%%%%%%%%%%%%%%%%%%%
    
    % partly implemented by using this website:
    % http://hplgit.github.io/Programming-for-Computations/pub/p4c/._p4c-bootstrap-Python028.html#4th:NonlinAlgEq:Secant
    
    iter = 0; % initialize iteration count variable
    i = 0; % initialize counter for values vector
    while(abs(f(x1)) > tol) % check for convergence
       x = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0)); % iteration
       err(i+1) = x; % add value to err vector
       i = i + 1; % increase value vector counter size
       x0 = x1; 
       x1 = x;
       iter = iter + 1; % increase iteration count
    end
    r = x; % return estimated root
    values = err; % return values at every iteration
    iters = iter; % return number of iterations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r] = quadratic(a,b,c) % quadratic function implementation
    %%%%%%%%%%%%%%%%%%%%
    % Inputs:
    % a: coefficient 
    % b: coefficient
    % c: coefficient
    
    % Outputs:
    % r: vector of actual roots
    %%%%%%%%%%%%%%%%%%%%
    
    % implemented by using Numerical Analysis 2nd Edition by Sauer page 17
    
    r1 = (-b + sqrt(b^2-4*a*c))/(2*a);
    r2 = (-b - sqrt(b^2-4*a*c))/(2*a);
    
    r = [r1 r2];
end