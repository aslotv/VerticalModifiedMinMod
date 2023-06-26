% Calculates the net growth in all state variables.
function s = netGrowth(un,growth,mu,Y,I,f,delta,k,P,R,rho,DOC_HCZ,DOC_F)
% This function calculates the right hand side of equations S1 to S15, in
% the Supporting Information. 
%
%   Input: 
%   un - The solutions u at t = t_n. 
%
%   mu - A structure containing the maximum growth rates for B, A, D, H, C,
%   Z, and F, as well as the growth rates for B, A, and
%   D.
%   
%   growth - The growth in osmo- and phagotrophs.
% 
%   Y - The yields. 
% 
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as the ingestion rates, and the
%   denominator from the expression defining them.
%
%   f - The photosynthetic carbon overflow, the fraction of total loss that
%   enters detritus, and the fraction of total loss that is respired.
%
%   delta - Specific mortality rates. 
% 
%   k - Dissolution, leakage and fragmentation rates, as well as empirical
%   coefficients related to chlorophyll light absorption.
% 
%   P - Photosynthetic quotient, inorganic phosphate concentration at deep
%   boundary (z_max), and particulate losses.
%   
%   R - Respiratory quotient, and respiration losses. 
% 
%   rho - The stoichiometric ratios and conversion factors.
%   
%   DOC_HCZ - DOC losses to L, from H, C, and Z. 
% 
%   DOC_F - DOC loss to L, from F. 
%
%   Output:
%   
%   s - The net growth in all state variables.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

%% Pre calculation set up

% Storage vectors and matrices:
s = zeros(size(un)); % Net growth of all variables.

%% Net osmotroph growth

%%% Bacteria (B):
s(:,1) = growth(:,1) - I.H.*un(:,4); 

%%% Autotrophic flagellate (A):
s(:,2) =  growth(:,2)- I.CA.*un(:,5); 

%%% Diatom (D):
s(:,3) = (mu.D - delta.D).*un(:,3) - I.ZD.*un(:,6); 

%% Net phagotroph growth

%%% Heterotrophic flagellate (H):
s(:,4) = growth(:,4) - I.CH.*un(:,5);

%%% Ciliate (C):
s(:,5) = growth(:,5) - I.ZC.*un(:,6);

%%% Mesozooplankton (Z):
s(:,6) = growth(:,6) - I.F.*un(:,7) - delta.Z.*un(:,6); 

%%% Fish (F):
s(:,7) = growth(:,7) - delta.F.*un(:,7); 
%% Net detritus growth

%%% Slow sinking detritus (Det_s):
s(:,8) = P.HC + k.frag.*un(:,9) - k.frag.*un(:,8);

%%% Fast sinking detritus (Det_f):
s(:,9) = delta.D.*un(:,3) + P.Z - k.frag.*un(:,9) ...
         + delta.F.*un(:,7) + P.F;

%%% Non-sinking detritus (Det_n):
s(:,10) = k.frag.*un(:,8) - k.l.*un(:,10);
%% Net growth in P, L, S, S_opal, and O

%%% Inorganic phosphate (P):
s(:,11) = -sum(s(:,1:10),2);

%%% Labile dissolved organic carbon (L):
s(:,12) = f.coc.*rho.CP.*sum(growth(:,2:3),2) ...
        + k.l.*rho.CP.*un(:,10) + DOC_HCZ ...
        - growth(:,1)/Y.BL + DOC_F;

%%% Silicate (S):
s(:,13) = k.opal.*un(:,14) - rho.DS.*growth(:,3);

%%% Particulate organic silicate (S_opal):
s(:,14) = rho.DS.*delta.D.*un(:,3) ...
        + rho.DS.*I.ZD.*un(:,6) - k.opal.*un(:,14);

%%% Dissolved oxygen (O):
s(:,15) = (1 + f.coc).*sum(growth(:,2:3),2).*rho.CP.*P.Q ...
    - (R.B + R.HCZ + R.F)./R.Q;
  
end
