
% Date : 16/01/2024
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : Main_code.m
%  Sub functions   : None

% Description::
% This function calculates the dispersion curves

% Input :
%       w  - frequency vector
%       kz - filtered eigenvalues (wavenumbers)
%       cs - S-wave velocity vector of finite layers
%       cs_hs - S-wave velocity of half-space

% Output :
%        v - sorted phase velocity
%        v1- unsorted phase velocity
%%
function [v, v1] = disp_curve(w, kz, cs, cs_hs)

kz(kz==0) = NaN;
f         = repmat(w, size(kz,1), 1);
f         = f .* 2 * pi;   % frequency in rad/s
v         = f./kz;         % Phase velocity
 
v(v>max([cs cs_hs])-1) = NaN;      % neglecting waves with velocity higher than maximum S-wave velocity

vmin      = 0.5 * min([cs cs_hs]); % neglecting waves with velocity lower than half of minimum S-wave velocity
v(v<vmin) = NaN;
v1        = v;
v         = sort(v);              % Sorting the phase velocity from minumum to maximum
v         = v';
v=real(v);

end
%%