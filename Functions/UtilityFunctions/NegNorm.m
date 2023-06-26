function negNorm = NegNorm(x,type)
%NEGNORM This function takes a weighted norm of the array x, where all
%positive elements have zero weight and all other elements have unit
%weight. 
%
%   Input: 
%
%   x - The array for which the weighted norm is calculated. 
%
%   type - The type of norm to take on the array x; must be one of the
%   types accepted by the MATLAB(R) function NORM (If x is a vector: any
%   positive scalar, Inf, or -Inf. If x is a matrix: 1, 2, Inf, or 'fro').
%
%   Output: 
%
%   negNorm - The weighted norm of x, where all positive elements have zero
%   weight and all other elemens have unit weight. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 


% Defining any positive elements as 0, i.e., giving them 0
% weight in the calculation of this norm:
x(x>0) = 0;

% Calculating the maximum absolute row sum of the negative values:
negNorm = norm(x,type);

end