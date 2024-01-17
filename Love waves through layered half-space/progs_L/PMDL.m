
% Date : 16/01/2024
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : Main_code.m
%  Sub functions   : None

% Description::
% This function calculates the stiffness matrices for elastic half-space
% and merge with the layer stiffness matrices.

% Input :
%       Layer stiffness matrices (A, G, M)
%       cs_hs - shear wave velocity of half-space
%       ro_hs - density of the half-space
%       PML - number of PMDL element

% Output :
%       Stiffness matrices of the model (A, G, M)
%      Lpml - thickness of half-space layers

%%

function [A, G, M, Lpml] = PMDL(A, G, M, cs_hs, ro_hs, PML)

mu   = ro_hs .* cs_hs.^2; % Lamé's first parameter
Dxx  = mu;
GP   = [0 2];             % Mid-point integration co-ordinate (0) and weight (2)

% thickness of each PMDL layers
% for details please refer to the corresponding paragraph in the paper.

% Lpml = 1*2.2.^(0:PML-1);   % for m scale profiles 
Lpml = 1 * 2.^(0:PML-1); %   % for km level crustal

sz   = length(A)+PML;
rc   = length(A):1:sz;
rc1  = rc + 1;

A(sz,sz) = 0; G(sz,sz) = 0; M(sz,sz) = 0;

s    = GP(1); wt = GP(2);
N    = [(1-s) (1+s)]./2;
Bs   = [-1 1]./2;

for k = 1:PML
    
    L  = Lpml(k); L  = L/2; B  = Bs/L;
    pp = rc(k) : rc1(k);
    %% Calculate elements of equation (15) in paper
    A(pp,pp) = A(pp,pp)   + (N.' * mu * N) * (L * wt);
    G(pp,pp) = G(pp,pp)   + (B.' * Dxx * B) * (L * wt);
    M(pp,pp) = M(pp,pp)   + (N.' * eye(1) * N) * (L * wt * ro_hs);
    
end

end

%%