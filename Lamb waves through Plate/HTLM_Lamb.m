%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program :: HTLM_Lamb.m
%
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

% Last revision date:
% 16/01/2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "HTLM_Lamb.m" is a Matlab script to generate dispersion curves of
% Lamb wave using higher order thin layer method. Users are encouraged 
% to refer to the corresponding research paper for a detailed explanation 
% of the underlying concepts and methods.

% Manuscript title –Higher-order thin layer method as an efficient forward model for calculating 
% dispersion curves of surface and Lamb waves in layered media

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:

    % fmin     – Starting frequency (Hz).
    % fmax     – maximum frequency (Hz) up to which the dispersion spectra are required.
    % df       – the frequency resolution (Hz).
    % d        – is the order of HTLM elements. The method can take from linear (d=1) to 15th order (d = 15).
    % dh       – is the thickness (m) of the thin discretized element
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


data = xlsread('Profile_Plate.xlsx');         % Input model parameters

h = data(1:end,1); vs = data(:,2); vp = data(:,3); rho = data(:,4);

%% %%%%%%%%%%%%%%%%%%%%%%%%%   Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin    = 0;                    % Starting frequency
fmax    = 30000;                % Maximum required frequency
df      = 100;                  % Frequency resolution
d       = 6;                    % Order of HTLM (1-15),(1-Linear, 2-Quadratic, 3-Cubic, 4-Quartic)
dh      = 0.05;                 % Thickness of the thin layers. 
w       = fmin : df : fmax;     % Frequency scale
fm      = [20000];              % Plot mode shape at frequency "fm"; else, fm = [], fm <= fmax
mode_n  = [3];                  % mode index. plots nth mode at frequency "fm"

%% %%%%%%%%%%%%%%%%%%%%%%%%   Main Code   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('progs_Lamb');

tic

   [v, ev_mat, dof] = Main_code(vs, vp, rho, h, w, d, dh, fm);

toc;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dispersion curve_________________________________________________________

figure; plot(w, v, '-');

xlabel('Frequency (Hz)','FontSize', 16,'fontweight','bold'); 
ylabel('Phase velocity (m/s)','FontSize', 16,'fontweight','bold');
set(gca,'fontname','times','TickDir','out'); 
box on
title('Dispersion curves')

%% Mode shape_______________________________________________________________
if fm > 0 
    mode_num   = size(ev_mat,2)-1;
    if mode_n <= mode_num
        figure('Position',[586 446 323 420])
        plot(ev_mat(:,mode_n),ev_mat(:,end),'-k'); axis ij;
        set(gca, 'XAxisLocation', 'top','TickDir', 'out')
        h(isnan(h)) = [];
        maxdepth    = sum(h);
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