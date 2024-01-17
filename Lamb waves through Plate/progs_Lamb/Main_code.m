

% Date : 13/07/2023
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function working :

%           stiffness.m
%           Eigen.m
%           eigen_vec
%           disp_curve.m


% Input :
%       Layer parameters ( vs, vp, rho, h)
%       w - frequency vector
%       d - Order of shape function polynomial
%       dh - layer thickness
% Output :

%       v - Dispersion curve
%       ev_mat - eigenvector matrix

function[v, ev_mat, dof] = Main_code(vs, vp, rho, h, w, d, dh, fm)

disp('Running HTLM_Lamb');
disp('Coded by Mrinal Bhaumik and Tarun Naskar, IIT Madras, India');


%% Checking for row and column's orientation (need all values in single row)
if ~isrow(vs)
    vs = vs';
end
if ~isrow(vp)
    vp = vp';
end
if ~isrow(rho)
    rho = rho';
end
if ~isrow(h)
    h = h';
end

%% Poisson's ratio to P-wave velocity  conversion

if vp < 0.5
    vp = vs.*sqrt(2.*(1-vp)/(1-2.*vp));
end
%%

h(isnan(h)) = [];                 % deleting the NaN value
nElem       = (h./dh);            % Calculating number of thin layers in each layer

nElem (1)      = ceil(nElem(1)); 
nElem(nElem<1) = 1;               % Making sure the minimum number of thin layer is 1
nElem          =   round(nElem);  % Rounding off the layer numbers
each_layer_div_in_m = h ./ nElem; % The thickness of thin layer in each individual layer
dpth = repelem(each_layer_div_in_m ./ d, nElem .* d); % Calculates the distance between each internal nodes

% nElem is number of thin layer in each layer, each thin layer has internal
% nodes depending upon the order of the interpolation used ( see Figure 2
% in the paper for better clarity)

%% Calculation and globalization of stiffness matrices of finite layers

[A, B, C, M] = stiffness(vs, vp, rho, nElem, h, d);

dof = length(M);  % Degree of freedom

%% Arranging the horizontal DOF first then vertical

Z  = zeros(size(A,1) / 2);
x  = 1 : 2 : size(A,1)-1;
z  = 2 : 2 : size(A,1);

% z- odd points (Horizontal); y-even points (Vertical); Refer to equation
% (18) in the paper
K2 = ([A(x,x) Z; -B(z,x) A(z,z)]); 
K0 = ([C(x,x) B(x,z);Z C(z,z)]);
M  = ([M(x,x) Z;Z M(z,z)]);

%% Eigen value and eigenvector(at frequency fm) calculation

[kz, evi] = Eigen(w, M, K0, K2, fm);

% sorting the eigenvalues

kz(kz==Inf)               = 0;
kz(real(kz) < 0)          = 0;   % eliminating -ve real part
kz(imag(kz) > 10^-4)      = 0;   % eliminating +ve complex part
kz(abs(imag(kz)) > 10^-4) = 0;   % eliminating -ve complex part

%% Dispersion curve plot

[v,v1] = disp_curve(w, kz, vp);

%% Calculation of mode shape at frequency fm
ev_mat = [];

if fm > 1
    
    dpth = [dpth dpth(end)];
    ev_mat = eigen_vec(evi, v1, w, fm, dpth);
    
end

end
