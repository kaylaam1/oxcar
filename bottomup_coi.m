
function [flux_model] = bottomup_coi(CRR_mmol_yr)
% bottomup_coi: Computes the depth-integrated carbon export flux through a
% water column
% 
%
% This function integrates the Carbon Remineralization Rate (CRR) from the 
% deepest ocean layer upward, resulting in the carbon export flux at each depth.
%
% Input:
%   CRR_mmol_yr - 3D array (m x n x p) of Carbon Remineralization Rate 
%                 in mmol C/m^3/year, where p represents depth levels.
%
% Output:
%   flux_model  - 3D array (m x n x p) of cumulative carbon export flux 
%                 (mmol C/m^3/year), computed from the bottom upward.
%
% Notes:
% - The flux at each depth level includes all remineralization occurring below it.
% - The flux at the deepest layer is equal to the local CRR value.
% - The summation proceeds upward, accumulating flux contributions.

    % Get the dimensions of the input CRR array
    [m, n, p] = size(CRR_mmol_yr);
    
    % Initialize flux_model array with zeros
    flux_model = zeros(m, n, p); 

    % Set the bottom boundary condition (deepest layer flux equals local CRR)
    flux_model(:,:,end) = CRR_mmol_yr(:,:,end);

    % Compute cumulative flux moving upward through the depth levels
    for k = p-1:-1:1 % Iterate from the second-to-last layer to the surface
        flux_model(:,:,k) = flux_model(:,:,k+1) + CRR_mmol_yr(:,:,k);
    end
end
