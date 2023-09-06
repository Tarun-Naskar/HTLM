%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program :: HTLM.m
%
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

% Last revision date:
% 06/09/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "HTLM.m" is a Matlab script to generate dispersion spectra of
% Rayleigh and Lamb wave using higher order thin layer method.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:

    % fmax     – maximum frequency (Hz) up to which the dispersion spectra are required.
    % df       – the frequency resolution (Hz).
    % d        – is the order of HTLM elements. The method can take from linear (d=1) to 15th order (d = 15).
    % dh       – is the thickness (m) of the thin discretized element
    % type     - "type = 1" for plate, "type = 2" for soil
    % PML      – Number of PMDL layer require to model the half-space. For soil model (m scale) 10 PMDL layer is sufficient. For crustal model (km scale) 15 is recommended.
    % L_min    – Calculates the minimum possible wavelength presents in the model.
    % fm       – is the frequency (Hz) where mode shape is required. Other wise keep it empty [].
    % mode_n   - mode index.

% Output:

    % v        – is the matrix containing phase velocities of all the possible modes
    % ev_mat   – is the matrix containing mode shapes at frequency “fm”
    % dof      – provide the total degree of freedom (except half-space), involves in the calculation.

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
%    h | vs | vp/mu | rho


data = xlsread('Profile_I.xlsx');         % Input model parameters

h = data(1:end,1)'; vs = data(:,2)'; vp = data(:,3)'; rho = data(:,4)';

%% %%%%%%%%%%%%%%%%%%%%%%%%%   Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmax    = 50;                   % Maximum required frequency
df      = 1;                  % Frequency resolution
d       = 8;                   % Order of HTLM (1-15),(1-Linear, 2-Quadratic, 3-Cubic, 4-Quartic)
dh      = 10;                   % Thickness of the thin layers. 

type    = 2;                    % "type = 1" for plate, "type = 2" for soil
PML     = [10] ;                  % No of PMDL layer (default value 10)

w       = df : df : fmax;       % Frequency scale
L_min   = min(vs)/fmax;         % Minimum wavelength
fm      = [20];                 % Plot mode shape at frequency "fm"; else, fm = [], fm <= fmax
mode_n  = [9];                  % mode index. plots nth mode at frequency "fm"

%% %%%%%%%%%%%%%%%%%%%%%%%%   Main Code   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output : v - Phase velocities
%     ev_mat - eigenvector or mode shape of given frequency fm


addpath('progs');

tic

   [v,ev_mat,dof] = Main_code(vs,vp,rho,h,w,d,dh,PML,type,fmax,fm);

toc;





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dispersion curve_________________________________________________________

figure; plot(w(1:1:end),v(1:1:end,:),'-k');

xlabel('Frequency (Hz)','FontSize', 16,'fontweight','bold'); 
ylabel('Phase velocity (m/s)','FontSize', 16,'fontweight','bold');
set(gca,'fontname','times','TickDir','out'); 
box on
title('Dispersion curves')

%% Mode sape_______________________________________________________________
if fm > 0 
    mode_num = size(ev_mat,2)-1;
    if mode_n <= mode_num
        figure('Position',[586 446 323 420])
        plot(ev_mat(:,mode_n),ev_mat(:,end),'-k'); axis ij;
        set(gca, 'XAxisLocation', 'top','TickDir', 'out')
        h(isnan(h))=[];
        maxdepth = 2 * sum(h);
        ylim ([0 maxdepth])
        xline(0,'--k')
        title(['Mode: ',num2str(mode_n)]);
        ylabel('Depth (m)','FontSize', 16,'fontweight','bold');
        xlabel('Amplitude','FontSize', 16,'fontweight','bold');
        set(gca,'fontname','times','TickDir','out'); 
    else
        disp(['available number of modes: ', num2str(mode_num)]);
        error('the mode index "mode_n" exceeds the available number of modes at frequency "fm"')
    end
    
end

%%