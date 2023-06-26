function [c_t,A,b,d,p] = ERK435Tableu
%RK45TABLEU This function defines an embedded Runge-Kutta tableu as the
%vector c-tilde, the matrix A, and the vectors b and d. This tableu
%was given by Prof. Thor Soerevik in an obligatory assignment for the
%course MAT260 - Numerical Solution of Differential Equations, during the
%spring semester of 2015. 
%
%   Input: 
%   
%   This function takes no input. 
%
%   Output: 
%
%   c_t - The nodes of the higher order method.
%
%   A - The stage weights matrix.
% 
%   b - The scheme weigths of the lower order method.
% 
%   d - The difference between the scheme weights of the higher order method
%   and the scheme weights of the lower order method embedded within. 
%
%   p - The order of the lower order method. 

c_t = [0; 1/3; 1/3; 1/2; 1]; % The nodes of the higher order method.

% The stage weights matrix:
A = zeros(5,5); % A matrix with zeros where other values are not specified.
A(2:5,1) = [1/3; 1/6; 1/8; 1/2]; % Non-zero elements of first column.
A(3,2) = 1/6; % Non-zero element of second column. 
A(4:5,3) = [3/8; -3/2]; % Non-zero elements of third column.
A(5,4) = 2; % Non-zero element of fourth column. 

b = A(5,:); % The scheme weights of the lower order method (b = [b 0]). 

b_t = [1/6 0 0 2/3 1/6]; % The cheme weights of the higher order method

d = b_t - b; % The weights used to approximate the local error. 

p = 3; % The order of the lower order method. 
end

