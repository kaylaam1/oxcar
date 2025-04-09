clear all; % Clear workspace

% Load full biogeochemical (BGC) model oxygen output
load('CTL_He_PCO_Gamma0_kl12h_O5_POC2DIC_GM15_Nowicki_npp1_aveTeu_diffSig_O2C_uniEta_DICrmAnthro_2L_Pnormal_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00.mat');
full_oxygen_model = data.O2;

% Load regridded World Ocean Atlas (WOA) oxygen data
load('O2_WOA_regrid_mn.mat');

% Load OCIM model output (Ocean Circulation Inverse Model)
load('CTL.mat');

% Extract grid and mask information
grd = output.grid;  % Grid structure
msk = output.M3d;   % 3D ocean mask
iwet = find(msk(:)); % Indices of ocean points
spa = 365.25 * 24 * 60^2; % Seconds in a year

% Transport matrix (s⁻¹)
TRdiv = -output.TR / spa;

% Set up parameters structure
par.TRdiv = TRdiv; 
par.iwet = iwet;
par.nwet = length(iwet);
par.grd = grd;
par.msk = msk;
par.zw = grd.zw;
% Compute area-weighted depth volume
dAt = grd.DXT3d .* grd.DYT3d; 
dVt = dAt .* grd.DZT3d;
par.dVt = dVt;
%%
% Load Environmental Data (Temperature, Salinity, and Gas Transfer)
temp = msk; temp(iwet) = output.T; par.temp = temp;
salt = msk; salt(iwet) = output.S; par.sal = salt;

% Load gas transfer velocity and atmospheric pressure
tempdata = load('tempPC_C13_C14_1850_frac_0.6');
par.kw = tempdata.par.kw; % Gas transfer velocity
par.atm_pres = tempdata.par.P; % Atmospheric pressure

% Load Optimized Model Parameters
load('b_opt_feb.mat');
load('A_opt_feb.mat');
load('z_ref_opt_feb.mat');

% % Assign respiration model parameters
% par.z_ref = -100; % Reference depth (m)
% par.b = b_opt_feb;
% par.A = A_opt_feb * 2.625; % Scaling factor

par.z_ref = -grd.zt(3);
par.b = b_opt_feb * 1;
par.A = A_opt_feb * 2.725;

% par.z_ref = -grd.zt(4);
% par.b = b_opt_feb * 0.801;
% par.A = A_opt_feb * 1.52;

% Load Net Primary Production (NPP) Data
NPP = ncread('biopump_model_output_Nowicki.nc','/NPP'); % NPP (mmolC/m²/yr)
par.NPP = NPP(:,:,1) / spa; % Convert to mmolC/m2/s (note-- need convert to m3 before using to scale production probably)

%% **Model Configuration 1: Oxygen Model with Gas Exchange, Transport, and Respiration Only**
[f, Rz, J] = O2simple_with_remin_before_prod(par);
dO_config1 = par.msk + NaN; % Initialize output
dO_config1(iwet) = J \ f;   % Solve for steady-state oxygen concentration

% Compute Oxygen Utilization Rate (OUR)
OUR = dO_config1 .* Rz; 
OUR(~iwet) = NaN;

%% **Model Configuration 2: Add Production Proportional to Integrated Respired Oxygen**

% Initialize iterative solver parameters
tol = 1e-6; 
max_iters = 20;
FJ = decomposition(J); % Use LU decomposition to speed up solving

A_prev = dO_config1; % Initial condition
mean_diff_all = NaN(1, max_iters);
sum_diff_all = NaN(1, max_iters);
max_diff_all = NaN(1, max_iters);

for iter = 1:max_iters
    % Compute mean remineralization below reference depth
    i_above_zo = find(-grd.zt > par.z_ref);
    msk2 = par.msk; msk2(:,:,i_above_zo) = 0;
    dVt2 = dVt .* msk2;
    
    OUR_vec = Rz(:) .* A_prev(:);
    mean_OUR = sum((dVt2(:) / sum(dVt2(:), 'omitnan')) .* OUR_vec, 'omitnan');

    % Define production term (P)
    msk3 = par.msk;
    below_zref = find(-grd.zt < par.z_ref);
    msk3(:,:,min(below_zref):end) = 0;
    P = msk3 * mean_OUR;
    P3 = P(iwet);

    % Solve for updated oxygen concentration
    A_new = par.msk + NaN;
    A_new(iwet) = FJ \ (f + P3);

    % Compute convergence metrics
    max_diff_all(iter) = max(abs(A_new(:) - A_prev(:)), [], 'omitnan');
    mean_diff_all(iter) = mean(abs(A_new(:) - A_prev(:)), 'omitnan');
    sum_diff_all(iter) = sum(abs(A_new(:) - A_prev(:)), 'omitnan');

    A_prev = A_new; % Update for next iteration

    fprintf('Iteration %d: max diff = %.5e, mean diff = %.5e\n', iter, max_diff_all(iter), mean_diff_all(iter));
end

dO_config2 = A_new; % Solution from Model Configuration 2

%% **Model Configuration 3: Distribute Production Based on NPP Field**
% Compute spatially varying production using satellite-derived NPP


A_prev = dO_config1; % Initial condition
mean_diff_all = NaN(1, max_iters);
sum_diff_all = NaN(1, max_iters);
max_diff_all = NaN(1, max_iters);
% Normalize NPP to maintain total production
NPP_masked = par.NPP .* par.msk(:,:,1);
sum_zt = sum(dAt .* msk3, 3, 'omitnan'); % reminder: msk3 is msk with 0s below z_ref
NPP_masked = NPP_masked .* sum_zt; %mmol/m^2/s to mmol/s oxygen produced
%which one? above or below
sum_zt = sum(grd.DZT3d .* msk3, 3, 'omitnan');
NPP_masked = NPP_masked ./ sum_zt;

% Define NPP scaling function
S = 1.34 * NPP_masked; % Redfield ratio applied (NPP is production of organic carbon and we want production of oxygen)
%can probably remove bc it's just a scalar?
S_norm = S / sum(S(:), 'omitnan'); % Normalize


for iter = 1:max_iters
% Compute distributed production

   % Compute mean remineralization below reference depth
    i_above_zo = find(-grd.zt > par.z_ref);
    msk2 = par.msk; msk2(:,:,i_above_zo) = 0;
    dVt2 = dVt .* msk2;
    
    OUR_vec = Rz(:) .* A_prev(:);
    mean_OUR = sum((dVt2(:) / sum(dVt2(:), 'omitnan')) .* OUR_vec, 'omitnan');


P_total = sum(msk3(:) * mean_OUR, 'omitnan');
P_scaled = P_total * S_norm;
num_wet_layers = sum(msk3, 3, 'omitnan');
num_wet_layers(num_wet_layers == 0) = NaN; % Avoid division by zero

% Distribute P_scaled evenly across wet layers
P_3D = msk3;
for k = 1:size(msk3, 3)
    P_3D(:,:,k) = (P_scaled ./ num_wet_layers) .* msk3(:,:,k);
end

P3 = P_3D(iwet);
A_new3 = par.msk + NaN;
A_new3(iwet) = FJ \ (f + P3);

    % Compute convergence metrics
    max_diff_all(iter) = max(abs(A_new(:) - A_prev(:)), [], 'omitnan');
    mean_diff_all(iter) = mean(abs(A_new(:) - A_prev(:)), 'omitnan');
    sum_diff_all(iter) = sum(abs(A_new(:) - A_prev(:)), 'omitnan');

    A_prev = A_new; % Update for next iteration

    % fprintf('Iteration %d: max diff = %.5e, mean diff = %.5e\n', iter, max_diff_all(iter), mean_diff_all(iter));
end
%second parameter set converged at 11
dO_config3 = A_new; % Solution from Model Configuration 3

%% Compute Model Performance Metrics
for i = 1
% Compute global volume-weighted statistics and print results
O2_WOA_vec = O2_WOA_regrid_mn(:);

dO_vec = dO_config1(:);
dO_vec2 = dO_config2(:);
dO_vec3 = dO_config3(:);
full_oxygen_vec = full_oxygen_model(:);
dVt_vec = dVt(:);

% Normalize volume weights (so they sum to 1)
dVt_norm = dVt_vec / sum(dVt_vec, 'omitnan');

% Compute weighted mean
mean_WOA = sum(dVt_norm .* O2_WOA_vec, 'omitnan');
mean_dO = sum(dVt_norm .* dO_vec, 'omitnan');
mean_dO2 = sum(dVt_norm .* dO_vec2, 'omitnan');
mean_dO3 = sum(dVt_norm .* dO_vec3, 'omitnan');
mean_full = sum(dVt_norm .* full_oxygen_vec, 'omitnan');

% Compute weighted standard deviation
std_WOA = sqrt(sum(dVt_norm .* (O2_WOA_vec - mean_WOA).^2, 'omitnan'));
std_dO = sqrt(sum(dVt_norm .* (dO_vec - mean_dO).^2, 'omitnan'));
std_dO2 = sqrt(sum(dVt_norm .* (dO_vec2 - mean_dO2).^2, 'omitnan'));
std_dO3 = sqrt(sum(dVt_norm .* (dO_vec3 - mean_dO3).^2, 'omitnan'));
std_full = sqrt(sum(dVt_norm .* (full_oxygen_vec - mean_full).^2, 'omitnan'));

% Compute weighted correlation
weighted_corr_optimized = sum(dVt_norm .* ((O2_WOA_vec - mean_WOA) .* (dO_vec - mean_dO)), 'omitnan') / (std_WOA * std_dO);
weighted_corr_optimized2 = sum(dVt_norm .* ((O2_WOA_vec - mean_WOA) .* (dO_vec2 - mean_dO2)), 'omitnan') / (std_WOA * std_dO2);
weighted_corr_optimized3 = sum(dVt_norm .* ((O2_WOA_vec - mean_WOA) .* (dO_vec3 - mean_dO3)), 'omitnan') / (std_WOA * std_dO3);
weighted_corr_fullBGC = sum(dVt_norm .* ((O2_WOA_vec - mean_WOA) .* (full_oxygen_vec - mean_full)), 'omitnan') / (std_WOA * std_full);

% Compute Global Oxygen Inventory (Total Oxygen Content in mmol)
total_oxygen_WOA = sum(O2_WOA_regrid_mn(:) .* dVt(:), 'omitnan');
total_oxygen_optimized = sum(dO_config1(:) .* dVt(:), 'omitnan');
total_oxygen_optimized2 = sum(dO_vec2(:) .* dVt(:), 'omitnan');
total_oxygen_optimized3 = sum(dO_vec3(:) .* dVt(:), 'omitnan');
total_oxygen_fullBGC = sum(full_oxygen_model(:) .* dVt(:), 'omitnan');

% Compute Volume-Weighted MSE
mse_fullBGC = sum(dVt_norm .* (full_oxygen_vec - O2_WOA_vec).^2, 'omitnan');
mse_optimized = sum(dVt_norm .* (dO_vec - O2_WOA_vec).^2, 'omitnan');
mse_optimized2 = sum(dVt_norm .* (dO_vec2 - O2_WOA_vec).^2, 'omitnan');
mse_optimized3 = sum(dVt_norm .* (dO_vec3 - O2_WOA_vec).^2, 'omitnan');

end
%% CRR and Carbon Export (Model, Full Model & WOA)

for i = 1
        below_zref = find(-grd.zt < par.z_ref);
        OUR = dO_config1 .* Rz;
        OUR2 = dO_config2 .* Rz;
        OUR3 = dO_config3 .* Rz;
load('OUR_WOA_tf.mat');% OUR is Modeled Oxygen Utilization Rate 
OUR_WOA = OUR_WOA_tf; % (mmol/m^3/yr)

isdry = find(msk(:) == 0); OUR_WOA(isdry) = NaN;
OUR_WOA(OUR_WOA < 0) = 0; %save only oxygen utilization (consumption)

TF = isoutlier(OUR_WOA(:)); 
OUR_WOA(TF) = NaN;

ratio_O_C = 1.34; % Redfield ratio of -O2: C organic C (mol O2 per mol C)
mol_C_per_mol_O2 = 1 / ratio_O_C; % mol C per mol O2

% Compute Carbon Remineralization Rate
CRR = OUR * mol_C_per_mol_O2 * spa; % mmol C / m^3 year
CRR = CRR / 1000;

CRR2 = OUR2 * mol_C_per_mol_O2 * spa; % mmol C / m^3 year
CRR2 = CRR2 / 1000; %mol C/ m^3 yr

CRR3 = OUR3 * mol_C_per_mol_O2 * spa; % mmol C / m^3 year
CRR3 = CRR3 / 1000;
% Compute CRR for WOA
CRR_WOA = OUR_WOA * mol_C_per_mol_O2; % mmol C / m^3 year

% Global Integrated Carbon Export 
dVt(isdry) = NaN; %try setting to 0 instead?
CRR(:,:,1:min(below_zref)) = zeros(91,180,min(below_zref));
CRR_sum = sum(CRR(iwet).* dVt(iwet),'omitnan'); 
CRR_sum = CRR_sum  * 12 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C

CRR2(:,:,1:min(below_zref)) = zeros(91,180,min(below_zref));
CRR_sum2 = sum(CRR2(iwet).* dVt(iwet),'omitnan'); 
CRR_sum2 = CRR_sum2  * 12 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C

CRR3(:,:,1:min(below_zref)) = zeros(91,180,min(below_zref));
CRR_sum3 = sum(CRR3(iwet).* dVt(iwet),'omitnan'); 
CRR_sum3 = CRR_sum3  * 12 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C

CRR_WOA(isnan(CRR_WOA)) = 0;
CRR_WOA_sum = sum(CRR_WOA(iwet).* dVt(iwet),'omitnan'); 
CRR_WOA_sum = CRR_WOA_sum * 12 / 1000 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C

[CRR_sum CRR_sum2 CRR_sum3 CRR_WOA_sum]

% Compute carbon flux by depth level
CRR_mmol_yr = CRR .* dVt; % mol/ yr 
CRR_mmol_yr = CRR_mmol_yr * 1000; % mmol/ yr 
CRR_mmol_yr(isnan(CRR_mmol_yr)) = 0; % Set NaNs to zeros

flux_model = bottomup_coi(CRR_mmol_yr);

CRR_mmol_yr2 = CRR2 .* dVt; % mol/ yr 
CRR_mmol_yr2 = CRR_mmol_yr2 * 1000; % mmol/ yr 
CRR_mmol_yr2(isnan(CRR_mmol_yr2)) = 0; % Set NaNs to zeros

flux_model2 = bottomup_coi(CRR_mmol_yr2);

% Compute carbon flux by depth level
CRR_mmol_yr3 = CRR3 .* dVt; % mol/ yr 
CRR_mmol_yr3 = CRR_mmol_yr3 * 1000; % mmol/ yr 
CRR_mmol_yr3(isnan(CRR_mmol_yr3)) = 0; % Set NaNs to zeros

flux_model3 = bottomup_coi(CRR_mmol_yr3);
end
       below_zref = find(-grd.zt < par.z_ref);
  %      OUR = dO_store .* Rz;
        OUR = dO_config1 .* Rz;
load('OUR_WOA_tf.mat');% OUR is Modeled Oxygen Utilization Rate 
OUR_WOA = OUR_WOA_tf; % (mmol/m^3/yr)
%try computing AOU from WOA using oxygen not AOU (like use i guess our
%models temp and salinity fields? if not theirs?)
isdry = find(msk(:) == 0); OUR_WOA(isdry) = NaN;
OUR_WOA(OUR_WOA < 0) = 0; %save only oxygen utilization (consumption)

TF = isoutlier(OUR_WOA(:)); 
OUR_WOA(TF) = NaN;

ratio_O_C = 1.34; % Redfield ratio of -O2: C organic C (mol O2 per mol C)
mol_C_per_mol_O2 = 1 / ratio_O_C; % mol C per mol O2

% Compute Carbon Remineralization Rate
CRR = OUR * mol_C_per_mol_O2 * spa; % mmol C / m^3 year

CRR = CRR / 1000;
% Compute CRR for WOA
CRR_WOA = OUR_WOA * mol_C_per_mol_O2; % mmol C / m^3 year

% Global Integrated Carbon Export 
dVt(isdry) = NaN; %try setting to 0 instead?
CRR(:,:,1:min(below_zref)) = zeros(91,180,min(below_zref));
CRR_sum = sum(CRR(iwet).* dVt(iwet),'omitnan'); 
CRR_sum = CRR_sum  * 12 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C


CRR_WOA(isnan(CRR_WOA)) = 0;
CRR_WOA_sum = sum(CRR_WOA(iwet).* dVt(iwet),'omitnan'); 
CRR_WOA_sum = CRR_WOA_sum * 12 / 1000 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C

[CRR_sum CRR_WOA_sum]

% Compute carbon flux by depth level
CRR_mmol_yr = CRR .* dVt; % mol/ yr 
CRR_mmol_yr(isnan(CRR_mmol_yr)) = 0; % Set NaNs to zeros

% version 1: my function
flux_model = bottomup_coi(CRR_mmol_yr);

% version 2: matlab functions
CRR_mmol_yr = flip(CRR_mmol_yr,3);
flux_Pg_C_model = cumsum(CRR_mmol_yr,[3],'omitnan');
flux_Pg_C_model = flip(flux_Pg_C_model,3);

%mol C / m^2 / yr for comparison with reccap
% flux_model = flux_model / 1000; 
flux_model = flux_model ./ dAt;
flux_model(isdry) = NaN;


% Define depth levels to be plotted
depth_indices = [1, 2, 3, 4, 5, 6]; % Indices of depth levels to plot
depths = grd.zt(depth_indices); % Depths in meters corresponding to these indices
num_depths = length(depth_indices);

flux_model(isdry) = NaN;
flux_model = flux_model ./ dAt; %mmol C/m^2/yr

% Store results in a table
CRR_sum2 = 0; CRR_sum3 = 0; CRR_sum_full = 0; 

results = table(["Optimized with Resp"; "Optimized with Resp and Prod"; "Optimized with Resp and NPP scaled Prod"; "Full BGC"; "World Ocean Atlas"], ...
                [total_oxygen_optimized; total_oxygen_optimized2; total_oxygen_optimized3; total_oxygen_fullBGC; total_oxygen_WOA], ... % Added comma here
                [weighted_corr_optimized; weighted_corr_optimized2; weighted_corr_optimized3; weighted_corr_fullBGC; 0], ...
                [mse_optimized; mse_optimized2; mse_optimized3; mse_fullBGC; 0], ...
                [CRR_sum; CRR_sum2; CRR_sum3; CRR_sum_full; CRR_WOA_sum],'VariableNames', {'Model', 'Total Oxygen Inventory', 'Correlation', 'MSE','Total Carbon Export'});

disp(results)

% modify this one to take second table comparing different z_refs
% CRR_sum2 = 0; CRR_sum3 = 0; CRR_sum_full = 0; 
% total_oxygen_WOA = 0;
% results2 = table(["Optimized with Resp"; "Optimized with Resp and Prod"; "Optimized with Resp and NPP scaled Prod"; "Full BGC"; "World Ocean Atlas"], ...
%                 [total_oxygen_optimized; total_oxygen_optimized2; total_oxygen_optimized3; total_oxygen_fullBGC; total_oxygen_WOA], ... % Added comma here
%                 [weighted_corr_optimized; weighted_corr_optimized2; weighted_corr_optimized3; weighted_corr_fullBGC; 0], ...
%                 [mse_optimized; mse_optimized2; mse_optimized3; mse_fullBGC; 0], ...
%                 [CRR_sum; CRR_sum2; CRR_sum3; CRR_sum_full; CRR_WOA_sum],'VariableNames', {'Reference Depth', 'A', 'b', 'Total Oxygen Inventory', 'Correlation', 'MSE','Total Carbon Export'});
% 

%%
% flux_model2(isdry) = NaN;
% flux_model2 = flux_model2 ./ dAt; %mmol C/m^2/yr
% flux_model3(isdry) = NaN;
% flux_model3 = flux_model3 ./ dAt;
% figure
% for i = 1:num_depths
%     % Column 1: WOA Observational Data
%     subplot(num_depths, 2, (i-1)*2 + 1);
%     contourf(flux_model(:, :, depth_indices(i)));
%     title(['Carbon Export Resp only (', num2str(depths(i)), ' m)']);
%     colorbar; 
%     % clim([0 350]);
% 
%     % subplot(num_depths, 2, (i-1)*2 + 2);
%     % contourf(flux_model2(:, :, depth_indices(i)));   
%     % title(['Carbon Export Prod (', num2str(depths(i)), ' m)']);
%     % colorbar; 
%     % % clim([0 350]);
% end

% figure
% for i = 1:num_depths
%     % Column 1: WOA Observational Data
%     subplot(num_depths, 2, (i-1)*2 + 1);
%     contourf(flux_model2(:, :, depth_indices(i)));
%     title(['Carbon Export Resp only (', num2str(depths(i)), ' m)']);
%     colorbar; 
%     % clim([0 350]);
% 
%     subplot(num_depths, 2, (i-1)*2 + 2);
%     contourf(flux_model3(:, :, depth_indices(i)));   
%     title(['Carbon Export NPP-scaled Prod (', num2str(depths(i)), ' m)']);
%     colorbar; 
%     % clim([0 350]);
% end

% Define depth levels to be plotted
depth_indices = [2, 5, 8, 11, 14, 17]; % Indices of depth levels to plot
depths = grd.zt(depth_indices); % Depths in meters corresponding to these indices
num_depths = length(depth_indices);

% Create oxygen maps by depth with 2 columns (WOA vs. Optimized Model)

% figure
% for i = 1:num_depths
%     % Column 1: WOA Observational Data
%     subplot(num_depths, 2, (i-1)*2 + 1);
%     contourf(O2_WOA_regrid_mn(:, :, depth_indices(i)));
%     title(['WOA Observed O2 (', num2str(depths(i)), ' m)']);
%     colorbar; 
%     clim([0 350]);
% 
%     % Column 2: Optimized Model with Respiration Only
%     subplot(num_depths, 2, (i-1)*2 + 2);
%     contourf(dO(:, :, depth_indices(i)));   
%     title(['Optimized Modeled O2 (', num2str(depths(i)), ' m)']);
%     colorbar; 
%     clim([0 350]);
% end
% sgtitle('Comparison of Observed and Modeled Oxygen Concentrations Across Depths (mmol/m^3)');


% Define depth levels to be plotted
% depth_indices = [1, 2, 3, 4, 5, 6]; % Indices of depth levels to plot
depths = grd.zt(depth_indices); % Depths in meters corresponding to these indices
num_depths = length(depth_indices);
%%
figure
for i = 1:num_depths
    % Column 1: Optimized Model with Production
    subplot(num_depths, 2, (i-1)*2 + 1);
    contourf(O2_WOA_regrid_mn(:, :, depth_indices(i)));
    title(['World Ocean Atlas (', num2str(depths(i)), ' m)']);
    colorbar; 
    clim([0 350]);
    
    % Column 2: Optimized Model with NPP-scaled Production
    subplot(num_depths, 2, (i-1)*2 + 2);
    contourf(dO_config3(:, :, depth_indices(i)));
   
    title(['Optimized Model with NPP-scaled Production(', num2str(depths(i)), ' m)']);
    colorbar; 
    clim([0 350]);
end
sgtitle('Comparison of Observed and Modeled Oxygen Concentrations (mmol/m^3)');


%% Functions

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
%   par.z_ref  - Reference depth where remineralization starts
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

    % Compute the Jacobian matrix J including transport, gas exchange,
    % respiration
     J = par.TRdiv + KO2 + d0(Rz(iwet));
 
    % Compute the air-sea oxygen flux (f) as KO2 * oxygen saturation
    f = KO2 * o2sat; 

end

function X = d0(v) 
X = spdiags(v(:),0,length(v(:)),length(v(:)));

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
