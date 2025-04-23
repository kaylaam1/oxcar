% Optimization of Oxcar Parameters

% Load Optimized Model Parameters
load('optim_params_before_prod.mat'); % scalars for b, A, z_ref
par.b = optim_params_before_prod(1);
par.A = optim_params_before_prod(2);
par.z_ref = optim_params_before_prod(3);
% Initial guess for parameters: using previous optimized values
initial_guess = [par.b, par.A, par.z_ref];

% Bounds for parameters
lb = [inital_guess(1) * 0.3, inital_guess(2) * 0.3]; % Lower bounds for b, A
ub = [inital_guess(1) * 3, inital_guess(2)*3]; % Upper bounds for b, A

% Optimization options
options = optimset('Display', 'iter', 'TolFun', 1e-6);

% Objective function to minimize
obj_fun = @(params) oxygen_mse_with_zref(params, par, O2_WOA_regrid_mn, iwet);

% Perform optimization -- change to unc
[opt_params, opt_mse] = fmincon(obj_fun, initial_guess, [], [], [], [], lb, ub, [], options);

% Extract optimal values
b_opt = opt_params(1);
A_opt = opt_params(2);

% Update par structure with reoptimized values
par.A = A_opt;
par.b = b_opt;


function mse = oxygen_mse_with_zref(params, par, O2_WOA_regrid_mn, iwet)
    % Unpack parameters
    b = params(1);
    A = params(2);
    
    % Update par structure
    par.b_lat = b;
    par.A_lat = A;
    
    % Compute modeled oxygen field
    [f, ~, J] = O2simple_with_remin2(par);

    dO = par.msk + NaN;
    dO(iwet) = J \ f; % Solve linear system
    
    % Compute MSE
    diff = dO(iwet) - O2_WOA_regrid_mn(iwet);
    mse = mean(diff.^2, 'omitnan');
end

