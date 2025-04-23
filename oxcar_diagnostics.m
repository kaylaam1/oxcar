clear all

% Load oxcar model output
load('O2_output_config2.mat', 'O2_model'); % load oxygen field
load('O2_output_config2.mat', 'OUR'); % load Oxygen Utilization Rate

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
par.M3d = output.M3d;% 3D ocean mask
TRdiv = -output.TR / spa; % Transport matrix (s⁻¹)
dAt = grd.DXT3d .* grd.DYT3d; 
dVt = dAt .* grd.DZT3d;

% Set up parameters structure
par.TRdiv = TRdiv; 
par.iwet = iwet;
par.nwet = length(iwet);
par.grd = grd;
par.msk = msk;
par.zw = grd.zw;
par.dVt = dVt;

% Load Environmental Data (Temperature, Salinity, and Gas Transfer)
temp = msk; temp(iwet) = output.T; par.temp = temp;
salt = msk; salt(iwet) = output.S; par.sal = salt;

% Load gas transfer velocity and atmospheric pressure
tempdata = load('tempPC_C13_C14_1850_frac_0.6');
par.kw = tempdata.par.kw; % Gas transfer velocity
par.atm_pres = tempdata.par.P; % Atmospheric pressure

% Load Optimized Model Parameters
load('optim_params_before_prod.mat'); % scalars for b, A, z_ref
par.b = optim_params_before_prod(1);
par.A = optim_params_before_prod(2);
par.z_ref = optim_params_before_prod(3);

% Compute Model Performance Metrics
for i = 1
% Compute global volume-weighted statistics and print results
O2_WOA_vec = O2_WOA_regrid_mn(:);

dO_vec = O2_model(:);
dO_vec2 = O2_model(:);
dO_vec3 = O2_model(:);
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
total_oxygen_optimized = sum(O2_model(:) .* dVt(:), 'omitnan');
total_oxygen_optimized2 = sum(dO_vec2(:) .* dVt(:), 'omitnan');
total_oxygen_optimized3 = sum(dO_vec3(:) .* dVt(:), 'omitnan');
total_oxygen_fullBGC = sum(full_oxygen_model(:) .* dVt(:), 'omitnan');

% Compute Volume-Weighted MSE
mse_fullBGC = sum(dVt_norm .* (full_oxygen_vec - O2_WOA_vec).^2, 'omitnan');
mse_optimized = sum(dVt_norm .* (dO_vec - O2_WOA_vec).^2, 'omitnan');
mse_optimized2 = sum(dVt_norm .* (dO_vec2 - O2_WOA_vec).^2, 'omitnan');
mse_optimized3 = sum(dVt_norm .* (dO_vec3 - O2_WOA_vec).^2, 'omitnan');

% zonal_bias = squeeze(mean(O2_model - O2_WOA_regrid_mn, 2, 'omitnan'));
% figure; contourf(grd.yt, grd.zt, zonal_bias'); 
% set(gca, 'YDir','reverse'); 
% xlabel('Latitude'); ylabel('Depth (m)');
% title('Zonal Mean Oxygen Bias (Model - WOA)');
% colorbar;

end

% CRR and Carbon Export (Model, Full Model & WOA)

for i = 1
        below_zref = find(-grd.zt < par.z_ref);
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

CRR2 = OUR * mol_C_per_mol_O2 * spa; % mmol C / m^3 year
CRR2 = CRR2 / 1000; %mol C/ m^3 yr

CRR3 = OUR * mol_C_per_mol_O2 * spa; % mmol C / m^3 year
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
% 
% load('OUR_WOA_tf.mat');% OUR is Modeled Oxygen Utilization Rate 
% OUR_WOA = OUR_WOA_tf; % (mmol/m^3/yr)
% %try computing AOU from WOA using oxygen not AOU (like use i guess our
% %models temp and salinity fields? if not theirs?)
% isdry = find(msk(:) == 0); OUR_WOA(isdry) = NaN;
% OUR_WOA(OUR_WOA < 0) = 0; %save only oxygen utilization (consumption)
% 
% TF = isoutlier(OUR_WOA(:)); 
% OUR_WOA(TF) = NaN;
% 
% ratio_O_C = 1.34; % Redfield ratio of -O2: C organic C (mol O2 per mol C)
% mol_C_per_mol_O2 = 1 / ratio_O_C; % mol C per mol O2
% 
% % Compute Carbon Remineralization Rate
% CRR = OUR * mol_C_per_mol_O2 * spa; % mmol C / m^3 year
% 
% CRR = CRR / 1000;
% % Compute CRR for WOA
% CRR_WOA = OUR_WOA * mol_C_per_mol_O2; % mmol C / m^3 year
% 
% % Global Integrated Carbon Export 
% dVt(isdry) = NaN; %try setting to 0 instead?
% CRR(:,:,1:min(below_zref)) = zeros(91,180,min(below_zref));
% CRR_sum = sum(CRR(iwet).* dVt(iwet),'omitnan'); 
% CRR_sum = CRR_sum  * 12 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C
% 
% 
% CRR_WOA(isnan(CRR_WOA)) = 0;
% CRR_WOA_sum = sum(CRR_WOA(iwet).* dVt(iwet),'omitnan'); 
% CRR_WOA_sum = CRR_WOA_sum * 12 / 1000 / (10^15); % mmol/yr to Pg C /yr, 1 mol C: 12 g C
% 
% [CRR_sum CRR_WOA_sum]

% Compute carbon flux by depth level
CRR_mmol_yr = CRR .* dVt; % mol/ yr 
CRR_mmol_yr(isnan(CRR_mmol_yr)) = 0; % Set NaNs to zeros

flux_model = bottomup_coi(CRR_mmol_yr);
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

% Define depth levels to be plotted
depth_indices = [2, 5, 8, 11, 14, 17]; % Indices of depth levels to plot
depths = grd.zt(depth_indices); % Depths in meters corresponding to these indices
num_depths = length(depth_indices);

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
    contourf(O2_model(:, :, depth_indices(i)));
   
    title(['Optimized Model with NPP-scaled Production(', num2str(depths(i)), ' m)']);
    colorbar; 
    clim([0 350]);
end
sgtitle('Comparison of Observed and Modeled Oxygen Concentrations (mmol/m^3)');


% 
% NEED to include volume weighting before doing some of this new bias stuff and residual
% analysis
% zonal_bias = squeeze(mean(O2_model - O2_WOA_regrid_mn, 2, 'omitnan'));
% figure; contourf(grd.yt, grd.zt, zonal_bias'); 
% set(gca, 'YDir','reverse'); 
% xlabel('Latitude'); ylabel('Depth (m)');
% title('Zonal Mean Oxygen Bias (Model - WOA)');
% colorbar;
% % 
% % %download fxn from matlab file exchange
% % taylordiagram([O2_model(:), O2_model(:), O2_model(:)], O2_WOA_regrid_mn(:), ...
% %     'Labels', {'Resp Only', 'Resp+Prod', 'Resp+NPP Prod'});
% 
% CRR_integrated = squeeze(sum(CRR, 3, 'omitnan')) .* dAt(:,:,1);
% figure; contourf(CRR_integrated / 1e6); % Pg C / m²
% title('Depth-Integrated Carbon Remineralization (Pg C/m²/yr)');
% colorbar;
% 
% figure;
% histogram(O2_model(:) - O2_WOA_regrid_mn(:), 50);
% title('Histogram of Model - WOA Residuals');
% xlabel('Residual (mmol/m^3)'); ylabel('Frequency');
