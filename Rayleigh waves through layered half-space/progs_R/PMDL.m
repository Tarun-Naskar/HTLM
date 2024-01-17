
% Date : 16/01/2024
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : stiffness.m
%  Sub functions   : None

% Description::
% This function calculates the stiffness matrices for elastic half-space
% and merge with the layer stiffness matrices.

% Input :
%       Layer stiffness matrices (Kxx, Kxz, Kzz, M)
%       Halpf-space parameters (cs_hs, cp_hs, ro_hs)

% Output :
%       Stiffness matrices of the model

%%

function [K0, K2, M, Lpml] = PMDL(A, B, C, M, cs_hs, cp_hs, ro_hs, PML)

mu   = ro_hs .* cs_hs.^2;            % Lamé's first parameter
Lam  =  ro_hs .* cp_hs.^2 - 2 * mu;  % Lamé's second parameter

% Material sub-matrices. Please refer to equation(5) of the paper

Kxx  = [Lam+2*mu 0 ; 0 mu];
Kxz  = [0 Lam ;  mu 0];
Kzz  = [mu 0 ; 0 Lam+2*mu];

GP   = [0 2];                 % Mid-point integration co-ordinate (0) and weight (2)

% thickness of each PMDL layers
% for details please refer to the corresponding paragraph in the paper.
Lpml = 1*2.^(0:PML-1);         % for m scale profiles 
% Lpml = 1.5*2.2.^(0:PML-1);   % for km level crustal


sz   = length(C) + 2*PML;

% Calculating some index values to perform globalization
rc   = length(C)-1:2:sz; 
rc1  = rc + 3;

% Allocating space for the stiffness matrices
A(sz,sz)=0; B(sz,sz)=0; C(sz,sz) = 0; M(sz,sz)=0;

p    = GP(1); wt = GP(2);
N    = [(1-p) (1+p)]./2;
N    = kron(N,eye(2));
Bs   = [-1 1]./2;

for k = 1 : PML
    
    L  = Lpml(k);
    L  = L/2;
    N_  = Bs/L;
    N_  = kron(N_,eye(2));
    pp = rc(k):rc1(k);
    %% Calculate elements of equation (11) in paper
    A(pp,pp) = A(pp,pp)   + (N.' * Kxx * N) * (L * wt);
    B(pp,pp) = B(pp,pp)   - (N.' * Kxz * N_) * (L * wt);
    B(pp,pp) = B(pp,pp)   + (N_.' * Kxz.' * N) * L * wt;
    C(pp,pp) = C(pp,pp)   + (N_.' * Kzz * N_) * (L * wt);
    M(pp,pp)   = M(pp,pp) + (N.' * eye(2) * N) * (L * wt * ro_hs);
    
end

A = A(1:end-2,1:end-2); B = B(1:end-2,1:end-2);
C = C(1:end-2,1:end-2); M = M(1:end-2,1:end-2);

%% %% Arranging the horizontal DOF first then vertica
Z  = zeros(size(A,1)/2);
x  = 1 : 2 : size(A,1)-1;
z  = 2 : 2 : size(A,1);

% z- odd points (Horizontal); y-even points (Vertical); Refer to equation
% (18) in the paper

K2 = ([A(x,x) Z; -B(z,x) A(z,z)]); 
K0 = ([C(x,x) B(x,z);Z C(z,z)]);
M  = ([M(x,x) Z;Z M(z,z)]);

end