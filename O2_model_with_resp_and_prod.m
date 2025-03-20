function [f, J, Rz, Pxy] = O2_model_with_resp(par)
% O2_model_with_remin: The function constructs a system of equations for solving the 
% steady-state oxygen field in the ocean, including air-sea oxygen fluxes, physical transport, 
% oxygen consumption due to the respiration of organic matter, and oxygen production due to photosynthesis.
%
% Steady-state Oxygen Model Equation:
% [∇∙(u+K∙∇) + KO2 + Rz - Pxy] * [dO] = KO2 * O2_sat
%
% In Matrix form:
%     [         J           ]  * [dO] =      f 
% 
% Inputs:
%   par - A structure containing grid information, the transport matrix, and model parameters (A, b, and z_ref).
%   The adjustable parameters are:
%   par.A      - Scaling constant for remineralization flux
%   par.b      - Exponent for remineralization power-law
%   par.z_ref  - Reference depth where remineralization starts
%
% Outputs:
%   f       - Air-sea flux vector (gas transfer velocity * oxygen saturation) (mmol/m^3 s)
%   J       - Jacobian matrix incorporating transport, air-sea gas exchange, and oxygen consumption (s^-1)
%   Rz      - Rate coefficient of oxygen utilization due to respiration (s^-1)
%   Pxy     - Rate coefficient of oxygen production due to photosynthesis (s^-1)
%
% Notes:
% - The model assumes oxygen utilization follows a first order loss term with the rate coefficient 
%   given by the derivative of a depth-dependent power-law function below a reference depth (z_ref).
% - The function constructs a linear system of equations for solving the steady-state oxygen field.
% - Gas exchange at the surface acts as a boundary condition, restoring oxygen toward equilibrium.
% - Oxygen Utilization Rates (mmol/m^3 s) can be computed after solving the model by multiplying modeled oxygen (mmol/m^3) 
%   by Rz, the rate coefficient of oxygen consumption (s^-1).
% - Oxygen production (Pxy) is assumed to be proportional to oxygen respiration (Rz) integrated from the ocean bottom to z_ref

    % Compute air-sea gas exchange parameters (KO2: gas transfer velocity, o2sat: oxygen saturation)
    [KO2, o2sat] = Fsea2air_v2(par); 
    
    % Extract grid information
    M3d = par.M3d; % Wet-dry mask
    msk = par.msk;
    dzt = par.grd.dzt; % Layer thickness at each depth
    zt = -par.grd.zt;  % Convert depth to negative values (below surface)
    zw = -par.zw;      % Convert interface depths to negative values
    iwet = par.iwet;   % Index of wet grid points (ocean cells)
    isdry = find(msk(:) == 0); % Index of dry grid points

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
    % Copy flux profile horizontally and apply ocean mask
    j_z = reshape(j_z, [1 1 nz+1]);  
    j_z = repmat(j_z, [ny, nx, 1]); 
    M3d_bottom = M3d;
    M3d_bottom(:,:,end+1) = M3d(:,:,end);
    dzt = reshape(dzt, [1 1 nz]);
    dzt = repmat(dzt, [ny, nx, 1]);
    j_z(M3d_bottom == 0) = 0;
      
    % The respiration rate coefficient, Rz, is assumed to be proportional to the carbon 
    % remineralization rate & oxygen availability
    Rz = -((j_z(:,:,1:end-1) - j_z(:,:,2:end)) ./ dzt); % Take derivative of j_z to get the carbon flux divergence (representing respiration)
    i_above_zo = find(zt > z_ref); % index for depths shallower than z_ref
    Rz(:,:,i_above_zo) = 0; % Apply respiration only below the z_ref
    Rz(isdry) = 0;

    % Compute rate coefficient of Oxygen Production Pxy to incorporate into the Jacobian J
    % Production is assumed to be proportionate to Oxygen Respiration (Rz * [dO]) integrated from the ocean bottom to z_ref
    % Oxygen Production (mmol/m^3/s) = ∫(lb: z_bottom, ub: z_ref) Rz * [dO] dz  
    %                                = Pxy * [dO] = ∫(lb: z_bottom, ub: z_ref) Rz dz 

    P = sum(Rz,3); % Production is the sum of respiration rate coefficient in each water column 
    % from the deepest grid cell to z_ref (all other cells outside of this range in Rz are already zero)
    sumP = sum(P(:)); % Take sum of production rate coefficients in every wet grid box above z_ref
   
    % Only apply production to wet points above z_ref
    msk2 = msk; % Make a copy of msk
    below_zref = find(zt < z_ref);
    msk2(:,:,min(below_zref):end) = 0; % zero out points below z_ref
    above_zref_wet = find(msk2 == 1); % find wet points above z_ref
    meanP = sumP / (length(above_zref_wet)); % Distribute total production evenly across # wet boxes above z_ref
  
    % Create a production rate coefficient field Pxy
    Pxy = msk2;
    Pxy = Pxy * meanP; % multiply wet grid cells above z_ref (value = 1) by mean production rate coefficient
    for i = 1:length(dzt)
    Pxy(:,:,i) = Pxy(:,:,i) / dzt(i); % divide each horizontal layer above z_ref by thickness of grid cell 
    end
    % Note: make in-line function for d0
    % X = spdiags(v(:),0,length(v(:)),length(v(:)));
    J = par.TRdiv + KO2 + d0(Rz(iwet)) - d0(Pxy(iwet)); % Compute the Jacobian matrix J 
    f = KO2 * o2sat; % Compute the RHS of steady-state oxygen model equation (f)
  
end

