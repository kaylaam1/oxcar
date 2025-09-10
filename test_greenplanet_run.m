%% Simplified Test Version for Greenplanet, Loading in our Group's model data and my own files
% 
%test_greenplanet_run.m
% Also want it to compute and store a .mat file with results that i can
% access

% this is a Driver that calls another function to solve the model
% (simplified right now for the test to just take input and provide the
% output)

% Ocean Biogeochemical Model Parameters

clear all; clear; clc;

% ------------------------------------------------------------------------


% OCIM (Ocean Circulation Inverse Model) output
load('CTL.mat');   
grd = output.grid;      % grid structure
msk = output.M3d;       % 3D ocean mask
iwet = find(msk(:));    % indices of wet (ocean) points

% Convert TR units to s^-1
spa = 365.25 * 24 * 60^2; 
TRdiv = -output.TR / spa;

% Store relevant fields in parameter struct
par = struct();
par.M3d = msk;
par.TRdiv = TRdiv;
par.iwet = iwet;
par.nwet = length(iwet);
par.grd = grd;
par.zw = grd.zw;

% Grid volumes and areas
dAt = grd.DXT3d .* grd.DYT3d; 
dVt = dAt .* grd.DZT3d; 
par.dVt = dVt;

% Environmental fields
temp = msk; temp(iwet) = output.T; 
salt = msk; salt(iwet) = output.S;
par.temp = temp; 
par.sal  = salt;

%----
% Load precomputed Bayesian data and model outputs
% This will come after (i will uncomment), first testing if it can load in
% research group's files that exist on greenplanet
% load('bayes_corrected_config3_50_OUR.mat','OUR_3D_vals');
% load('bayes_corrected_config3_50.mat','sumO2_vals','log_likelihood','pA_vals','b_vals');

% % Regridded WOA oxygen
% load('O2_WOA_regrid_mn.mat'); %another local file
%----


result = zeros(5,1);
par.z_ref   = -grd.zt(3);   % reference depth for fitting
test_vals = [1 2 3 4 5];
for k = 1:length(test_vals(:))

        par.A = test_vals(k)
        %solve model
        output = test_model_solver(par);

        result(k) = output * k;

end



%test function that it calls, but as its own script test_model_solver.m

% function output = test_model_solver(par)
% 
%     iwet = par.iwet;
%     msk = par.msk;
%     output_i = NaN(size(msk)); %initialize solution for O2
% 
%     % if strcmp(config, 'config1') 
%         % --- config1: no production, only respiration ---
%         output_i = par.A * sum(iwet(:));
% output = output_i;
% end
% 
 
