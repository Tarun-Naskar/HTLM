% Date : 13/07/2023
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : Main_code.m
%  Sub functions   : None

% Description::
% This function calculates the eigenvalue and eigenvectors

% Input :
%       Stiffness matrices

% Output :
%       Eigenvalues (kz) and eigenvectors (evi)


function[kz, evi] = Eigen(w, M, K0, K2, fmode)


kz  = 1i.*zeros( length(M), length(w));
%     evi = zeros( length(M), length(M), length(w));
evi = [];
%% Eigen value problem

for ii = 1 : length(w)
    
    f        = 2 * pi * w(ii);
    [kzi]    = eig(K0-f^2*M, -K2);
    %         [ev,kzi]  = eig(K0-f^2*M, -K2); kzi = diag(kzi); evi(:,:,ii) = ev;
    kzi      = sqrt((kzi));
    kz(:,ii) = kzi;
    
end

% Calculating the eigenvector at given frequency fmode
if fmode > 1        
    for jj = 1:length(fmode)
        f         = 2 * pi * fmode(jj);
        [ev, ~]  = eig(K0-f^2*M, -K2);
        evi       = ev;
    end
    
end
end
