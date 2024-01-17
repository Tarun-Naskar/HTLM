

% Date : 13/12/2023
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : Main_code.m
%  Sub functions   : None

% Description::
% This function calculates the dispersion curves

% Input :
%       w  - frequency vector
%       kz - filtered eigenvalues (wavenumbers)
%       cp - P-eave velocity vector

% Output :
%        v - sorted phase velocity
%        v1- unsorted phase velocity
%%
function [v, v1] = disp_curve(w, kz, cp)

kz(kz==0) = NaN;
f = repmat(w,size(kz,1),1);
f = f.*2*pi;          % frequency in rad/s
v = f./kz;            % Phase velocity

v(v>max(cp)-1) = NaN; % neglecting waves with velocity higher than maximum P-wave velocity

v1 = v;
v  = sort(v);         % Sorting the phase velocity from minumum to maximum
v  = v';

end
%%