function [c_t,A,b,d,p] = ERK436Tableu
%ERK34 Summary of this function goes here
%   Detailed explanation goes here

% ARK4(3)6L[2]SA–ERK in Additive Runge–Kutta schemes for
% convection–diffusion–reaction equations, Kennedy et. al. 2003. 

c_t = [0; 1/2; 83/250; 31/50; 17/20; 1]; % The c-tilde vector.

A = zeros(6,6); % A matrix with zeros where other values are not specified.

% The non-zero elements in the first column: 
A(2:6,1) = [1/2; 13861/62500; -116923316275/2393684061468; ...
    -451086348788/2902428689909;  647845179188/3216320057751];

% The non-zero elements in the second column: 
A(3:6,2) = [6889/62500; -2731218467317/15368042101831;  ...
    -2682348792572/7519795681897; 73281519250/8382639484533]; 

% The non-zero elements in the third column: 
A(4:6,3) = [9408046702089/11113171139209; ...
    12662868775082/11960479115383; 552539513391/3454668386233];

% The non-zero elements in the fourth column:
A(5:6,4) = [3355817975965/11060851509271; 3354512671639/8306763924573];

% The non-zero elements in the fifth column:
A(6,5) = 4040/17871;

b = [2889/524892 0 15625/83664 69875/102672 -2260/8211 1/4];

% The b-tilde vector:
b_t = [4586570599/29645900160 0 178811875/945068544 814220225/1159782912 ...
    -3700637/11593932 61727/225920];

d = b_t - b; 

p = 3; % The order of the lower order method.

end

