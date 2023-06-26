% Calculates the growth of osmo- and phagotrophs. 
function growth = Growth(un,mu,Y,I,expectedIngestion,N_z)
%GROWTH is a function that calculates the growth in B, A, D, H, C, Z, and F
%for each time step. 
%
%   Input:
%
%   un - The solutions from the previous time step. 
%
%   mu - A structure containing the growth rates fo B,
%   A, and D, as well as the maximum growth rates for B, A, D, H, C, Z, and
%   F.
%
%   Y - A structure containing the yields for H, C, Z, F, and B. 
%
%   I - A structure containing the ingestion rates for
%   H, C, Z, and F, as well as the maximum ingestion rates for H, C, Z, and
%   F.
%
%   expectedIngestion - The expected value of the ingestion rate for F in
%   the water column if fish is subject to DVM, the depth specific
%   ingestion rate for F (I.F) if not. 
%
%   N_z - The number of depths in zspan. 
%
%   Output: 
%
%   growth - The growth in B, A, D, H, C, Z, and F in the current time
%   step. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Defining a storage matrix for the growth: 
growth = zeros(N_z,7); % Growth in osmo- and phagotrophs.

% The growth in the osmotrophs:
growth(:,1) =  mu.B.*un(:,1); % Growth in B.
growth(:,2) = mu.A.*un(:,2); % Growth in A.
growth(:,3) = mu.D.*un(:,3); % Growth in D.

% The growth in the phagotrophs:
growth(:,4) = Y.H.*I.H.*un(:,4); % Growth in H.
growth(:,5) = Y.C.*I.C.*un(:,5); % Growth in C.
growth(:,6) = Y.Z.*I.Z.*un(:,6); % Growth in Z.
growth(:,7) = Y.F.*expectedIngestion.*un(:,7); % Growth in F.
end
