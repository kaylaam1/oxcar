function output = test_model_solver(par)

    iwet = par.iwet;
    msk = par.msk;
    output_i = NaN(size(msk)); %initialize solution for O2

    % if strcmp(config, 'config1') 
        % --- config1: no production, only respiration ---
        output_i = par.A * sum(iwet(:));
output = output_i;
end
         