function [O2_solution, iter_metrics, OUR] = O2_model_solver(par, config, max_iters)
% Inputs:
%   par       - structure with all model parameters and grid info
%   config    - string, one of: 'config1', 'config2', 'config3'
%   max_iters - maximum # iterations for production term loop (not needed
%               for config1)
%  
%   Configuration Descriptions:
%   config1 includes physical transport, air-sea gas exchange, and respiration
%   config2 adds production source term, equal to the integrated respired oxygen below the reference depth
%   config3 uses NPP to scale spatially-varying production 
%   config4 varies b (or A), tbd, based on NPP/temp/cell size/ etc
%    
% Outputs:
%   O2_solution - final 3D oxygen field
%   iter_metrics - structure with iteration convergence info
%   OUR - Oxygen Utilization Rate: oxygen_solution .* Rz (respiration rate coefficent) 
%   Note: add output that combines positive source term for production
%   layers and negative loss in deeper layers?

% Compute outputs from model subfunction
[f, Rz, J] = O2simple_with_remin_before_prod(par); % See subfunction below
%   f       - Air-sea flux vector (gas transfer velocity * oxygen saturation) (mmol/m^3 s)
%   J       - Jacobian matrix incorporating transport, air-sea gas exchange, and oxygen consumption (s^-1)
%   Rz      - Rate coefficient of oxygen utilization due to respiration (s^-1)

FJ = decomposition(J); % decomposition of Jacobian for computational efficiency

iwet = par.iwet; msk = par.msk;
O2_solution = NaN(size(msk)); %initialize solution for O2
        if strcmp(config, 'config1') 
            % --- config1: no production, only respiration ---

            O2_solution(iwet) = FJ \ f; % solve for O2
            iter_metrics = []; % config 1 has no iteration metrics
        elseif strcmp(config, 'config2') || strcmp(config, 'config3') % config2 and config3 include production term
            A_prev = msk + NaN; % initialize solution for O2
            A_prev(iwet) = FJ \ f; % solve for O2
 
            tol = 1e-6; % tolerance for max difference in mmol/m^3 oxygen between A_prev and A_new
            iter_metrics = struct('max_diff', [], 'mean_diff', [], 'sum_diff', []); % initialize iter_metrics
            i_above_zo = find(-par.grd.zt > par.z_ref); % find depth layers above z_ref
            msk2 = msk; msk2(:,:,i_above_zo) = 0; % mask set to 0 for all layers above z_ref
            dVt2 = par.dVt .* msk2;  % 3D volumes set to 0 for all layers above z_ref
            dVt_weight2 = dVt2(iwet) / sum(dVt2(iwet)); % for computing volume-weighted mean OUR
            msk3 = msk; below_zref = find(-par.grd.zt < par.z_ref);
            msk3(:,:,min(below_zref):end) = 0; % Keep only surface layers
    
                dVt_prod = par.dVt .* msk3;
                vol_prod = sum(dVt_prod(iwet));  % Compute volume of surface layers (production layers)
                dVt_resp = par.dVt .* msk2;
                vol_resp = sum(dVt_resp(iwet)); % Compute volume of deep ocean (respiration layers)

            for iter = 1:max_iters
                % compute oxygen utilization rate
                OUR = A_prev .* Rz;  OUR(isnan(OUR)) = 0; %Note: check if i can i remove this (the reassign Nans to 0)? 
                % compute mean_OUR (volume-weighted)
                mean_OUR = sum(dVt_weight2 .* OUR(iwet));
              
    
            if strcmp(config, 'config2')
                % --- config2: uniform production ---
            
                % scale mean_OUR by ratio of respiration to production layers volume to get correct production rate
                scaled_mean_OUR = mean_OUR * (vol_resp / vol_prod);
                P = msk3 * scaled_mean_OUR;
                P3 = P(iwet);
            
            elseif strcmp(config, 'config3')
                % --- config3: spatially-varying production using NPP ---
                
                % expand NPP vertically (assumed uniform in depth in production layers)
                n_layers = min(below_zref) - 1;
                NPP_3D = repmat(par.NPP ./ n_layers, 1, 1, n_layers); % mmol/m²/s
            
                % normalize using dz
                dz = abs(diff(par.grd.zt));
                dz_prod = reshape(dz(1:n_layers), 1, 1, n_layers); 
                NPP_vol = NPP_3D ./ dz_prod;
            
                % mask onto ocean points
                NPP_vol(~logical(msk3(:,:,1:n_layers))) = 0;
            
                % compute total shape production
                shape_total = sum(NPP_vol(:) .* dVt_prod(:,:,1:n_layers), 'omitnan');
            
                % scale NPP_vol so that production = total respiration
                total_respiration = sum(dVt2(:) .* OUR(:), 'omitnan');
                P = (total_respiration / shape_total) * NPP_vol;
            
                % embed into full 3D array (fill deeper layers with 0)
                P_full = msk * 0;
                P_full(:,:,1:n_layers) = P;
                P = P_full;
                P3 = P_full(iwet); 
            end

            A_new = msk + NaN;
            A_new(iwet) = FJ \ (f + P3);

            % save iteration metrics
            iter_metrics.max_diff(iter) = max(abs(A_new(iwet) - A_prev(iwet)),[],'omitnan');
            % Optional other metrics, removed here for optimization to be fast
            % iter_metrics.mean_diff(iter) = mean(abs(A_new(:) - A_prev(:)), 'omitnan');
            % iter_metrics.sum_diff(iter) = sum(abs(A_new(:) - A_prev(:)), 'omitnan');

            A_prev = A_new;

            if iter_metrics.max_diff(iter) < tol % Stop iterating if max difference between solutions is below tolerance
                break;
            end
        end

        O2_solution = A_new; % Store solution after convergence
        OUR = O2_solution .* Rz; OUR(isnan(OUR)) = 0;

        else
        error('Unknown config. Choose ''config1'' or ''config2'' or ''config3''.');

    end
end


function [f, Rz, J] = O2simple_with_remin_before_prod(par)

% O2_model_with_remin: The function constructs a system of equations for solving the 
% steady-state oxygen field in the ocean, including air-sea oxygen fluxes,
% physical transport, and oxygen consumption due to the respiration of organic matter.
%
% Steady-state Oxygen Model Equation:
% [∇∙(u+K∙∇) + KO2 + Rz] * [dO] = KO2 * O2_sat
%
% In Matrix form:
%     [         J           ]  * [dO] =      f 
% 
% Inputs:
%   par - A structure containing grid information, the transport matrix, and model parameters (A, b, and z_ref).
%   The adjustable parameters are:
%   par.A      - Scaling constant for remineralization flux
%   par.b      - Exponent for remineralization power-law
%   The fixed parameters are:
%   par.z_ref  - Reference depth below which remineralization is applied (93.7 meters)
%
% Outputs:
%   f       - Air-sea flux vector (gas transfer velocity * oxygen saturation) (mmol/m^3 s)
%   J       - Jacobian matrix incorporating transport, air-sea gas exchange, and oxygen consumption (s^-1)
%   Rz      - Rate coefficient of oxygen utilization due to respiration (s^-1)
%
% Notes:
% - The model assumes oxygen utilization follows a first order loss term with the rate coefficient 
%   given by the derivative of a depth-dependent power-law function below a reference depth (z_ref).
% - The function constructs a linear system of equations for solving the steady-state oxygen field.
% - Gas exchange at the surface acts as a boundary condition, restoring oxygen toward equilibrium (assuming atmospheric concentration is fully saturated/O2_sat).
% - Oxygen Utilization Rates (mmol/m^3 s) can be computed after solving the model by multiplying modeled oxygen (mmol/m^3) 
%   by Rz, the rate coefficient of oxygen consumption (s^-1).
% - Oxygen production (Pxy) is assumed to be equal to oxygen respiration (Rz) integrated from the ocean bottom to z_ref

    % Compute air-sea gas exchange parameters (KO2: gas transfer velocity, o2sat: oxygen saturation)
    [KO2, o2sat] = Fsea2air_v2(par); 
    
    % Extract grid information
    M3d = par.M3d; % Wet-dry mask
    msk = par.msk; % ocean mask 
    dzt = par.grd.dzt; % Layer thickness at each depth
    zt = -par.grd.zt;  % Convert depth to negative values (below surface)
    zw = -par.zw;      % Convert interface depths to negative values
    iwet = par.iwet;   % Index of wet grid points (ocean cells)
    z_ref = par.z_ref;
    isdry = find(msk(:) == 0); % Index of dry grid points
    i_above_zo = find(zt > z_ref); % Index for depths shallower than z_ref
    % Note:move above out of function so doesn't have to repeat and make sure in all par structures
    
    % Respiration parameters
    % R(z) = (d/dz)(-A * (zw / z_ref).^(-b))  
    z_ref = par.z_ref; % Reference depth where remineralization starts
    b = par.b;         % Exponent for power-law remineralization
    A = par.A;         % Scaling constant for remineralization flux

    % Extend zw to include the bottom dry layer
    zw(end + 1) = zw(end) - dzt(end); 

    % Compute j_z, the power law representing organic carbon flux divided by an assumed O2:C ratio 
    % (negative because downward carbon flux)
    j_z = -A * (zw / z_ref).^(-b); 
    j_z(1) = [];  % Remove the first value (infinity at zw = 0) 
    j_z(end + 1) = 0; % Add zero at the bottom because the assumption is the flux that hits the sediments is respired in the deepest grid cell of the model

    [ny, nx, nz] = size(par.msk); % Grid size
    j_z = reshape(j_z, [1 1 nz+1]); % Copy flux profile horizontally
    j_z = repmat(j_z, [ny, nx, 1]); 
    M3d_bottom = M3d;
    M3d_bottom(:,:,end+1) = M3d(:,:,end); 
    dzt = reshape(dzt, [1 1 nz]);
    dzt = repmat(dzt, [ny, nx, 1]);
    j_z(M3d_bottom == 0) = 0; % Apply ocean mask
      
    % The respiration rate coefficient, Rz, is assumed to be equal to the carbon remineralization rate & oxygen availability
    Rz = -((j_z(:,:,1:end-1) - j_z(:,:,2:end)) ./ dzt); % Take derivative of j_z to get the carbon flux divergence (representing respiration)
    Rz(:,:,i_above_zo) = 0; % Apply respiration only below the z_ref
    Rz(isdry) = 0; % Apply respiration only at wet points
    d0 = @(x) spdiags(x(:),0,length(x(:)),length(x(:))); % sparse diagonal function
    J = par.TRdiv + KO2 + d0(Rz(iwet));  %the Jacobian matrix J including transport, gas exchange, and respiration
    % Compute the air-sea oxygen flux (f) as KO2 * oxygen saturation
    f = KO2 * o2sat; 

end


function [KO2, o2sat] = Fsea2air_v2(par)
    % Inputs: par structure with temperature, salinity, and other parameters
    grd  = par.grd;  % Retrieve grid structure
    msk  = par.msk;  % Retrieve 3D mask (ocean points)
    iwet = par.iwet; % Indices of wet grid boxes (ocean)
    T = par.temp;    % Temperature data
    S = par.sal;     % Salinity data

    % Initialize arrays for gas transfer velocity (kw) and atmospheric pressure
    kw = msk*0;
    kw(:,:,1) = par.kw; % Assign gas transfer velocity at surface layer only
    atm_pres = msk*0;
    atm_pres(:,:,1) = par.atm_pres; % Assign atmospheric pressure at surface layer

    KO2  = msk*0; % Initialize KO2 array for gas exchange operator

    % Compute Schmidt number (sco2) for oxygen in water, dependent on temperature
    sco2 = 1638.0 - 81.83*T + 1.483*T.^2 - 0.008004*T.^3;

    % Compute gas transfer velocity kw, scaled by Schmidt number
    kw = kw.*sqrt(660./sco2);
    KO2 = kw/grd.dzt(1); % Divide by surface layer thickness (depth)

    KO2(:,:,2:end) = 0; % Set gas transfer velocity to 0 below the surface

    KO2 = d0(KO2(iwet)); % Reshape KO2 into a 1D vector for wet points
    o2sat = 0*msk; % Initialize oxygen saturation array
    o2sat = 1000 * o2sato(T,S) .* atm_pres; % Compute oxygen saturation (mmol/m^3)
    o2sat = o2sat(iwet); % Extract values only for wet grid cells
end

