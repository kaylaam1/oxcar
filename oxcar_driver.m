% oxcar_driver: Runs oxygen model solver (solve_O2_model.m) given user
% input model configuration and max iterations
%
% Model solves for the 3D oxygen field and Oxygen Utilization Rate field
% Solver also returns iteration metrics if configuration includes the iterative production term

clear; clc;
% Load OCIM model output (Ocean Circulation Inverse Model)
load('CTL.mat');

% Extract grid and mask information
grd = output.grid;  % Grid structure
msk = output.M3d;   % 3D ocean mask
iwet = find(msk(:)); % Indices of ocean points
spa = 365.25 * 24 * 60^2; % Seconds in a year
par.M3d = output.M3d;
TRdiv = -output.TR / spa; % Transport matrix (s⁻¹)
% Set up parameters structure
par.TRdiv = TRdiv; par.iwet = iwet; par.nwet = length(iwet);
par.grd = grd; par.msk = msk; par.zw = grd.zw;

% Compute area-weighted depth volume
dAt = grd.DXT3d .* grd.DYT3d; 
dVt = dAt .* grd.DZT3d; par.dVt = dVt;

% Load Environmental Data (Temperature, Salinity, and Gas Transfer)
temp = msk; temp(iwet) = output.T; par.temp = temp; salt = msk; salt(iwet) = output.S; par.sal = salt;

% Load gas transfer velocity and atmospheric pressure
tempdata = load('tempPC_C13_C14_1850_frac_0.6');
par.kw = tempdata.par.kw; % Gas transfer velocity
par.atm_pres = tempdata.par.P; % Atmospheric pressure

% Load Optimized Model Parameters
load('optim_params_before_prod.mat');
par.b = optim_params_before_prod(1); par.A = optim_params_before_prod(2); par.z_ref = optim_params_before_prod(3);

% --- Run Model Solver ---
%   Configuration Descriptions:
%   config1 includes physical transport, air-sea gas exchange, and respiration
%   config2 adds production source term, proportionate to the integrated respired oxygen below the reference depth
%   config3 uses NPP to scale spatially-varying production 
%   config4 varies b (or A), tbd, based on NPP/temp/cell size/ etc
config = 'config2';  % Choose 'config1-5'
max_iters = 30; % could be Inf?
[O2_model, iter_metrics, OUR] = solve_O2_model(par, config, max_iters);

% --- Save or pass O2_model to diagnostics ---
save('O2_output_config2.mat', 'O2_model', 'iter_metrics','OUR');
