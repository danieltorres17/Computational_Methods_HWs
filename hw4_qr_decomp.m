%%
% Daniel Torres
% Homework 3a
% CSCI 3656 
% Fall 2019

close all; clear all; clc;

%% Part 1: n_k x n_k Vandermonde and tridiagonal matrices

% Create a Vandermonde matrix 
n = 11;
x = linspace(-1, 1, n);
Vandermonde = fliplr(vander(x)); % orients columns appropriately

% Create an n x n tridiagonal matrix
tri_diag_mat = full(gallery('tridiag', n, 1, -2, 1));

% create n_k vector 
iters = 6;
n_k = zeros(1,iters);
for k =1:iters
    n_k(k) = 2^k+1; 
end

% Create n_k x n_k Vandermonde matrix and n_k x n_k tridiagonal matrix
V_cond_numbers = zeros(1,iters); % preallocate space
tri_cond_numbers = zeros(1,iters); % preallocate space
for i = 1:length(n_k)
    
    % create Vandermonde matrices
    x = linspace(-1,1,n_k(i));
    V = fliplr(vander(x));
    V_cond_numbers(i) = cond(V); % find condition number and store it
    
    % create tridiagonal matrices
    tri_diag_mat = full(gallery('tridiag', n_k(i), 1, -2, 1));
    tri_cond_numbers(i) = cond(tri_diag_mat); % find condiion number 
    
end

figure(1);
loglog(n_k,V_cond_numbers,'-o','LineWidth',2);
hold on
loglog(n_k,tri_cond_numbers,'-o','LineWidth',2);
axis square
grid on
title('Condition Number');
xlabel('Matrix Dimension');
ylabel('Condition Number');
legend('Vandermonde','Tridiagonal','Location','Northwest');
legend boxoff

% It can be seen that the condition numbers for the Vandermonde matrices increase
% much faster than for the tridiagonal matrices. 

%% Part 2: 33 x 33 Vandermonde matrix

% Create 33 x 33 Vandermonde matrix
n_2 = 33;
x_2 = linspace(-1, 1, n_2);
A = fliplr(vander(x_2)); % orients columns appropriately
b = rand(33,1); % create random b vector
x_true = A \ b; % solve for true solution using backslash

%% Part 3: QR Decomposition, my_qr_solver function

x_est = my_qr_solver(A,b); % estimated solution with my qr solver function
rel_err = norm(x_true-x_est)./norm(x_true); % compute relative error
fprintf('The norm of the relative error was calculated to be: %f.\n\n',rel_err);

%% Part 4: Timing Study

time_iters = 25; % time study iterations
slash_time = zeros(1,time_iters); % preallocate space
qr_time = zeros(1,time_iters); % preallocate space
n = 30; % size of A matrix and b vector 
x = linspace(-1, 1, n);
A = fliplr(vander(x)); % orients columns appropriately

% for this time study I created a random vector b for every iteration and
% solve that using the same A matrix. 

for i = 1:time_iters 
   
   b = randn(n,1); % create random number vector b 
   tic; % start clock for backslash operation
   x_sol = A \ b; % solve system using backslash
   slash_time(i) = toc; % stop clock and store time in backslash vector
   
   tic; % start clock for my_qr_function 
   x_est = my_qr_solver(A,b); % solve system using my_qr_solver function
   qr_time(i) = toc; % stop clock and store time in qr vector
end

avg_slash_time = mean(slash_time); % find average duration time for backslash
avg_qr_time = mean(qr_time); % find average duration time for qr function

fprintf('The average duration for the backslash was: %f seconds.\n',avg_slash_time);
fprintf('The average duration for qr decomp was: %f seconds.\n',avg_qr_time);
% It can be seen that the backslash which actually uses LU decomp is faster
% than qr decomposition.

function [x] = my_qr_solver(A,b) % qr decomp implementation
%   INPUTS:
%   A: matrix
%   b: vector
%
%   OUTPUT:
%   x: estimated solution using QR decomposition

% used this link to find out how linsolve works: 
% https://www.mathworks.com/help/matlab/ref/linsolve.html

    [Q,R] = qr(A); % compute QR decomposition of A
    y = Q'*b; % solve for y
    x = linsolve(R,y); % solve for x using Matlab backsolver
end



