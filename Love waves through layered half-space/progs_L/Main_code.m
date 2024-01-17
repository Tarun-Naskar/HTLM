

% Date : 16/01/2024
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function working :

%           stiffness.m
%           PMDL.m
%           shape_fn.m


% Input :
%       Layer parameters ( vs, rho, h).
%       w   - frequency vector.
%       d   - Order of shape function.
%       dh  - layer thickness.
%       PML – Number of PMDL layer require to model the half-space.

% Output :

%       v     – is the matrix containing phase velocities of all the possible modes
%       EV    – is the matrix containing mode shapes
%      dh_vec – depth of the nodes
%       dof   - degree of freedom

function[v, EV, dh_vec, dof] = Main_code(vs, rho, h, w, d, dh, PML)

% Checking for row and column's orientation (need all values in single row)

if ~isrow(vs)
    vs = vs';
end
if ~isrow(rho)
    rho = rho';
end
if ~isrow(h)
    h = h';
end

disp('Running HTLM_L');
disp('Coded by Mrinal Bhaumik and Tarun Naskar, IIT Madras, India')


% Parameters of the finite layers
cs      =   vs(1:end-1);
roS     =   rho(1:end-1);


h(isnan(h))    = [];                  % deleting the NaN value
nElem          = (h./dh);             % Calculating number of thin layers in each layer

nElem (1)      = ceil(nElem(1));
nElem(nElem<1) = 1;                   % Making sure the minimum number of thin layer is 1
nElem          =   round(nElem);      % Rounding off the layer numbers
dh_n           = h./(nElem*d);
dh_vec         = repelem(dh_n,nElem*d); % Calculates the distance between each internal nodes

% nElem is number of thin layer in each layer, each thin layer has internal
% nodes depending upon the order of the interpolation used ( see Figure 2
% in the paper for better clarity)


% Half-space [using PMDLs]

cs_hs      =   vs(end) ;
rho_hs     =   rho(end);

%% Calculation and globalization of stiffness matrices of finite layers

[A, G, M]  = stiffness(cs, roS, nElem, h, d);

dof = length(M);

%% Adding Half-space stiffness matrix

[A, G, M, LPML]  = PMDL(A, G, M, cs_hs, rho_hs, PML);

dh_vec = [dh_vec LPML];
%% Solving the eigen value problem

kz = 1.*zeros(length(M), length(w));
EV = zeros(length(M), length(M), length(w));

for ii = 1 : length(w)
    
    f          = 2 * pi * w(ii);
    [ev, kzi]  = eig(G-f^2*M, -A);
    kzi        = sqrt(diag(kzi));
    kz(:,ii)   = kzi;            % storing the eigenvalues
    EV(:,:,ii) = ev;             % storing the eigenvectors
end

% Sorting the eigenvalues

kz(kz==Inf)               = 0;
kz(real(kz) < 0)          = 0;   % eliminating -ve real part
kz(imag(kz) > 10^-4)      = 0;   % eliminating +ve complex part
kz(abs(imag(kz)) > 10^-4) = 0;   % eliminating -ve complex part

%% Dispersion curve calculation

kz(kz==0) = NaN;
f         = repmat(w,size(kz,1),1);
f         = f.*2*pi;
v         = f./kz;

% limiting the phase velocities

vmin                   = 0.5 * min(cs);
v(v<vmin)              = NaN;   % neglecting waves with velocity lower than half of minimum S-wave velocity
vmax                   = max([cs cs_hs]);
v(v>vmax)              = NaN;   % neglecting waves with velocity higher than maximum S-wave velocity

end
