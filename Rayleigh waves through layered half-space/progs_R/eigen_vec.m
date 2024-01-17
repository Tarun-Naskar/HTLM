% Date : 16/01/2024
% Coded by: Mrinal Bhaumik and Tarun naskar
% Indian Institute of Technology Madras, India

%% Sub-function of : Main_code.m
%  Sub functions   : None

% Description::
% This function calculates the mode shape and depth

% Input :
%       evi - eigenvectors
%       v1  - unsorted phase velocity
%        w  -  frequency vector
%        fm  - frequency at which the mode shape needs to be calculated 
%        Lpml - PMDL layers of half-space thicknesses

% Output :
%       ev_mat - mode shape and depth
%%
function [ev_mat] = eigen_vec(evi, v1, w, fm, dpth, Lpml)

    v1(isnan(v1))=0;
    pos = find(fm==w);                    % position index of the frequency fm
    
    dpth = [dpth Lpml];
    dpth = cumsum(dpth);                  % Calculating depth of each nodes
    dpth = [0 dpth(1:end-1)]'; 
    pp = find(v1(:,pos));                 % index of modes at frequency fm (location: pos)
    vec = evi(:,pp);                      % extracting the corresponding eigenvectors
    
    % shape = vec(size(v1,1)/2+1:end,:);  % for vertical motion
    shape = vec(1:size(v1,1)/2,:);        % for horizontal motion
    
    shape = shape./max(abs(shape));
    ev_mat = [shape dpth];

end

%%