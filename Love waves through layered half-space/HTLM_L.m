%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program :: HTLM_L.m
%
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

% Last revision date:
% 16/01/2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "HTLM_L.m" is a Matlab script to generate dispersion spectra of
% Love wave (SH-wave or antiplane shear) using higher order thin layer method.
% Users are encouraged  to refer to the corresponding research paper for a 
% detailed explanation of the underlying concepts and methods.

% Manuscript title –Higher-order thin layer method as an efficient forward model for calculating 
% dispersion curves of surface and Lamb waves in layered media

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:

% data     - Contain the input model parameters (h, vs, vp).
% fmin     – minimum frequency (Hz).
% fmax     – maximum frequency (Hz) up to which the dispersion spectra are required.
% df       – the frequency resolution (Hz).
% d        – is the order of HTLM elements. The method can take from linear (d=1) to 15th order (d = 15).
% dh       – is the thickness (m) of the thin discretized element
% PML      – Number of PMDL layer require to model the half-space. 
%            For soil model (m scale) 10 PMDL layer is sufficient. For crustal model (km scale) 15 is recommended.


% Output:

% v        – is the matrix containing phase velocities of all the possible modes.
% EV       – is the matrix containing mode shapes. 
% dh_vec   – depth of the nodes.


%  READS

% 'Profile name file.xlsx'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________
%__________________________________________________________________________

clc
clear all
close all
format compact

%% %%%%%%%%%%%%%%%%%%%%%%%%%   Input Model   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Data Format
%    h | vs | rho


data = xlsread('Profile_I.xlsx');         % Input model parameters

h = data(1:end,1); vs = data(:,2);  rho = data(:,4);

%% %%%%%%%%%%%%%%%%%%%%%%%%%   Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin    = 0;                    % minimum frequency
fmax    = 50;                   % Maximum required frequency
df      = 0.5;                  % Frequency resolution

d       = 6;                    % Order of HTLM (1-15),(1-Linear, 2-Quadratic, 3-Cubic, 4-Quartic)
dh      = 10;                   % Thickness of the thin layers.

PML     = 10 ;                  % No of PMDL layer (default value 10)
w       = fmin : df : fmax;     % Frequency scale

%% %%%%%%%%%%%%%%%%%%%%%%%%%%   Main Code   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('progs_L');

tic

[v, EV, dh_vec, dof] = Main_code(vs, rho, h, w, d, dh, PML);

toc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dispersion curve_________________________________________________________

figure; plot(w, sort(v),'-');

xlabel ('Frequency (Hz)','FontSize', 16,'fontweight','bold');
ylabel ('Phase velocity (m/s)','FontSize', 16,'fontweight','bold');
set (gca,'fontname','times','TickDir','out');
title ('Love wave dispersion curves')


%% __________________________________________________________________________
%% __________________________________________________________________________
