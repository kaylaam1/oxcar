% Optimization of Oxcar Parameters
clear all

% Load Optimized Model Parameters
load('optim_params_after_prod.mat'); % scalars for b, A, z_ref
par.b = optim_params_after_prod(1);
par.A = optim_params_after_prod(2);
par.z_ref = optim_params_after_prod(3);
% Load OCIM model output 
load('CTL.mat');
% Load regridded World Ocean Atlas (WOA) oxygen data
load('O2_WOA_regrid_mn.mat');
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

config = 'config2';
par.A = 2e-6; par.b = 1.4;
% Initial guess for parameters: using previous optimized values
initial_guess = log([par.b, par.A]);

% Bounds for parameters
% (removed after changing from fmincon to fminunc)
% lb = [initial_guess(1) * 0.1, initial_guess(2) * 0.1]; % Lower bounds for b, A
% ub = [initial_guess(1) * 10, initial_guess(2)*10]; % Upper bounds for b, A


options = optimoptions('fminunc', ...
    'Display', 'iter', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-9, ...
    'StepTolerance', 1e-9, ...
    'Algorithm', 'quasi-newton'); 

max_iters = 100;

% Objective function to minimize
obj_fun = @(params) oxygen_mse_with_zref(params, par, O2_WOA_regrid_mn, iwet, config, max_iters);

% Perform optimization -
% [opt_params, opt_mse] = fmincon(obj_fun, initial_guess, [], [], [], [], lb, ub, [], options);
[log_opt_params, opt_mse] = fminunc(obj_fun, initial_guess, options);

% Extract optimal values
b_opt = exp(log_opt_params(1));
A_opt = exp(log_opt_params(2));

% Update par structure with reoptimized values
par.A = A_opt;
par.b = b_opt;

% Optional: save new optimized parameters
optim_params_unc = [b_opt, A_opt, -grd.zt(3)]; % scalars for b, A, z_ref
%save('optim_params_unc.mat','optim_params_unc'); 

%%
function mse = oxygen_mse_with_zref(params, par, O2_WOA_regrid_mn, iwet,config, max_iters)
    % Unpack parameters
    b = exp(params(1)); %log_params
    A = exp(params(2));
    
    % Update par structure
    par.b = b;
    par.A = A;
    
    % Compute modeled oxygen field
    [O2_model, iter_metrics, OUR] = O2_model_solver(par, config, max_iters);


    % Compute MSE
    diff = O2_model(iwet) - O2_WOA_regrid_mn(iwet);
    %volume weighted difference
    diff = diff .* par.dVt(iwet);  diff = diff / sum(par.dVt(iwet));
    mse = mean(diff.^2, 'omitnan');
end

