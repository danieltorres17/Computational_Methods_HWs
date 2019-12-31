%
% Daniel Torres
% Homework 5
% CSCI 3656 
% Fall 2019

close all; clear all; clc;

% Part 1: Generating Training Data
n1 = 7; % number of partitions in vector
theta = 1; % theta parameter for sigmoid function
x1 = linspace(-5,5,n1); % x-points vector
y1 = sigmoid(x1,theta); % get y-data from sigmoid function

% create table and display the x and y value pairs calculated above
T1 = table(x1',y1');
T1.Properties.VariableNames = {'x1', 'y1'};
disp(T1);

% Part 2: Training the Model

V = fliplr(vander(x1)); % create Vandermonde matrix using x1 vector
a = V \ y1'; % solve for polynomial coefficients 

% Part 3: Generating Testing Data

n3 = 101; % number of partitions in vector
x3 = linspace(-5,5,n3); % x-points vector
y3 = sigmoid(x3,theta); % get y-data from sigmoid function
y3_mean = mean(y3); % calculate y3 mean
y3_stdev = std(y3); % calculate y3 standard deviation

% Part 4: Computing the Testing Error

% Links below were used to implement certain functions:
% https://www.mathworks.com/help/matlab/ref/polyval.html
% https://www.mathworks.com/help/matlab/ref/fliplr.html
% https://www.mathworks.com/help/matlab/ref/vander.html

y_p = polyval(flip(a),x3); % evaluate 101 point x vector using coefficients
err = max(abs(y3-y_p)); % max error of difference between true value and estimated 
fprintf('Max error calculated for difference between true vs polynomial estimation when theta = 1: %f.\n',err);

% Part 5: Repeating steps 1-4 with theta = 10

theta2 = 10; % theta parameter for sigmoid function
y5 = sigmoid(x1,theta2); % get y-data from sigmoid function

% create table and display the x and y value pairs calculated above
T2 = table(x1',y5');
T2.Properties.VariableNames = {'x2', 'y2'};
disp(T2);

V5 = fliplr(vander(x1)); % create Vandermonde matrix using x1 vector
a5 = V \ y5'; % solve for polynomial coefficients

y6 = sigmoid(x3,theta2); % get y-data from sigmoid function
y6_mean = mean(y6); % calculate y3 mean
y6_stdev = std(y6_mean); % calculate y3 standard deviation

y_p5 = polyval(flip(a5),x3); % evaluate 101 point x vector using coefficients
err5 = max(abs(y6-y_p5)); % max error of difference between true value and estimated 
fprintf('Max error calculated for difference between true vs polynomial estimation when theta = 10: %f.\n',err5);

disp(['The error was much bigger when theta = 10 than when theta = 1. This is ', ...
    'due to the fact that the bigger oscillations in the first and last intervals.', ...
    ' I''m assuming that the higher the theta value the higher the error.\n']);

figure(1);
plot(x1,y1,'LineWidth',2); % part 1 plot
hold on
plot(x3,y3,'LineWidth',2); % part 2 plot
hold on
plot(x3,y_p,'LineWidth',2); % part 3 plot
hold on
plot(x1,y5,'LineWidth',2); % part 4 plot
hold on
plot(x3,y_p5,'LineWidth',2); % part 5 plot
grid on
axis square
title('True vs Polynomial Interpolation');
xlabel('x');
ylabel('f_{\theta}');
legend('\theta = 1, n = 7','\theta = 1, n = 101','\theta = 1, poly eval',...
    '\theta = 10, n = 7','\theta = 10, poly eval','Location','SouthEast');
legend boxoff

% Sigmoid Function
function [y] = sigmoid(x,theta)
% function returns the sigmoid value of an input vector with a given theta
%%%%%%%%%%%
%   INPUT:
%   x: x-value vector
%   theta: specified theta value
%
%   OUTPUT:
%   y: sigmoid values of vector x
%%%%%%%%%%%

    y = 1./(1+exp(-theta.*x));
    
end




