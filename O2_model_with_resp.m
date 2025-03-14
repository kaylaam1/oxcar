function [f, f_remin, J] = O2_model_with_resp(par)
% O2_model_with_remin: The function constructs a system of equations for solving the 
% steady-state oxygen field in the ocean, including air-sea oxygen fluxes, physical transport, 
% and oxygen consumption due to the respiration of organic matter.
%
% Steady-state Oxygen Model Equation:
% [∇∙(u+K∙∇) + KO2 + f_remin] * [dO] = KO2 * O2_sat
%
% In Matrix form:
%     {         J           } * [dO] =    { f }
% 
% Inputs:
%   par - A structure containing grid information, the transport matrix, and model parameters.
%
% Outputs:
%   f       - Air-sea flux term (gas transfer velocity * oxygen saturation)
%   f_remin - Rate coefficient of oxygen utilization due to respiration
%   J       - Jacobian matrix incorporating transport, air-sea gas exchange, and oxygen consumption.
%
% Notes:
% - The model assumes oxygen utilization follows a depth-dependent power-law function below a reference depth (z_ref).
% - The function constructs a system of equations for solving the steady-state oxygen field.
% - Gas exchange at the surface acts as a boundary condition, restoring oxygen toward equilibrium.
% - Oxygen Utilization Rates (mmol/m^3 s) can be computed after solving the model by multiplying modeled oxygen (mmol/m^3) 
%   by f_remin, the rate coefficient of oxygen consumption (s^-1).

    % Compute air-sea gas exchange parameters (KO2: gas transfer velocity, o2sat: oxygen saturation)
    [KO2, o2sat] = Fsea2air_v2(par); 
    
    % Extract grid information
    M3d = par.M3d;
    dzt = par.grd.dzt; % Layer thickness at each depth
    zt = -par.grd.zt;  % Convert depth to negative values (below surface)
    zw = -par.zw;      % Convert interface depths to negative values
    iwet = par.iwet;   % Index of wet grid points (ocean cells)

    % Respiration parameters
    % R(z) = (-A * (zw / z_ref).^(-b)) dz  
    z_ref = par.z_ref; % Reference depth where remineralization starts
    b = par.b;         % Exponent for power-law remineralization
    A = par.A;         % Scaling constant for remineralization flux

    % Extend zw to include the bottom layer
    zw(end + 1) = zw(end) - dzt(end); 

    % Compute j_z, the power law representing organic carbon flux divided by an assumed O2:C ratio 
    % (negative because downward carbon flux)
    j_z = -A * (zw / z_ref).^(-b); 
    j_z(1) = [];  % Remove the first value (infinity at zw = 0) 
    j_z(end + 1) = 0; % Add zero at the bottom

    [ny, nx, nz] = size(par.msk); % Grid size
    % Copy flux profile horizontally and apply ocean mask
    j_z = reshape(j_z, [1 1 nz+1]);  
    j_z = repmat(j_z, [ny, nx, 1]); 
    M3d_bottom = M3d;
    M3d_bottom(:,:,end+1) = M3d(:,:,end);
    dzt = reshape(dzt, [1 1 nz]);
    dzt = repmat(dzt, [ny, nx, 1]);
    j_z(M3d_bottom == 0) = 0;
     
    % Take derivative of j_z to get the carbon flux divergence (representing respiration). 
    % The respiration rate coefficient, R_z, is assumed to be proportional to the carbon 
    % remineralization rate & oxygen availability

    R_z = -((j_z(:,:,1:end-1) - j_z(:,:,2:end)) ./ dzt);
    i_above_zo = find(zt > z_ref); % index for depths shallower than z_ref
    R_z(:,:,i_above_zo) = 0; % Apply respiration only below the z_ref

    J = par.TRdiv + KO2 + d0(R_z(iwet)); % Compute the Jacobian matrix J 
    f = KO2 * o2sat; % Compute the RHS of steady-state oxygen model equation (f)
    f_remin = R_z; % Store respiration rate coefficent (s^-1)
end



function mse = oxygen_model_optim(params, par, O2_WOA_regrid_mn, iwet)
% oxygen_model_optim: Computes the mean squared error (MSE) between modeled 
% and observed (World Ocean Atlas) oxygen concentrations.
%
% This function serves as an objective function for parameter optimization,
% adjusting remineralization parameters (A, b, z_ref) to minimize the mismatch 
% between the modeled oxygen field and observational data.
%
% Inputs:
%   params  - Vector of optimization parameters:
%             params(1) = b      (power-law exponent)
%             params(2) = A      (scaling factor)
%             params(3) = -z_ref (reference depth; negative for depth)
%   par     - Structure containing grid information and transport terms.
%   O2_WOA_regrid_mn - Observed oxygen concentrations (from WOA), regridded to model resolution.
%   iwet    - Indices of wet (ocean) grid cells.
%
% Output:
%   mse - Mean squared error (MSE) between modeled and observed oxygen concentrations.
%
% Notes:
% - Parameters are exponentiated to enforce positivity constraints.
% - The modeled oxygen field is computed using `O2_model_with_resp()`.
% - The MSE is computed by weighting squared errors by ocean volume.

    % Exponentiate optimization parameters
    b = exp(params(1)); 
    A = exp(params(2)); 
    z_ref = -exp(params(3)); % Ensure reference depth remains negative

    % Update parameter structure with new values
    par.b = b;
    par.A = A;
    par.z_ref = z_ref;

    % Solve for steady-state oxygen distribution using updated parameters
    [f, ~, J] = O2_model_with_resp(par);
    
    % Initialize modeled oxygen field
    dO = par.msk + NaN; 
    dO(iwet) = J \ f; % Solve the linear system for oxygen concentration

    % Compute residuals between modeled and observed oxygen
    diff = dO(iwet) - O2_WOA_regrid_mn(iwet);

    % Compute mean squared error and Weight squared differences by the 
    % volume of each ocean grid box
    mse = sum(diff.^2 .* dVt(iwet)) / sum(dVt(iwet));

end
