function finalResults = ViscosityPredictor(moleFractions, LowT, StepT, HighT, useVFT, exportPath)
% ViscosityPredictor Runs viscosity prediction without App UI.
%
% Inputs:
%   moleFractions : 1x17 numeric vector representing mole fractions of solutes.
%   useVFT        : Boolean (true/false). If true, applies VFT fitting.
%   exportPath    : String representing the output Excel filename (e.g., 'Results.xlsx').
%
% Output:
%   finalResults  : A MATLAB table containing the prediction data.

    %% 1. Input Validation and Setup
    if length(moleFractions) ~= 17
        error('Input moleFractions must be a vector of length 17.');
    end
    
    % Solute Names (Hardcoded from original App)
    soluteNames = {
        'Water'
        'Methanol'
        'Ethylene Glycol'
        '1,2-Propanediol'
        '1,3-Butanediol'
        'Diethylene Glycol'
        'N,N-Dimethylethanolamine'
        '2-Methyl-2,4-Pentanediol'
        'Triethylene glycol (TER)'
        'Diethylene Glycol Monobutyl Ether'
        'Glycerol'
        'DMSO'
        '3-Methoxy-1,2-Propanediol'
        '1,2,4 butanetriol'
        '1,3-Propanediol'
        '2,3-Butanediol'
        '1,4-Butanediol'
    };

    % Temperature Settings (Defaulted from original App)

    Temp = (HighT:-StepT:LowT)';

    %% 2. Load Data
    % Ensure these .mat files are in your MATLAB Path
    try
        data_mixtures = load('mixture_data.mat', 'mixtures3');
        mixtures3 = data_mixtures.mixtures3;
        
        data_net = load('PINN.mat', 'net', 'inputMin', 'inputMax');
        netLR = data_net.net;
        inputMeans = data_net.inputMin;
        inputStds = data_net.inputMax - data_net.inputMin;
        
        data_Tg = load('Tg_fitted_VFT.mat','Tg_fit_all');
        Tg_fit_all_K = data_Tg.Tg_fit_all(:) + 273.15;
    catch ME
        error('Error loading .mat files. Ensure mixture_data.mat, PINN.mat, and Tg_fitted_VFT.mat are in the path. Details: %s', ME.message);
    end

    %% 3. Calculate Mass Fractions
    Atomic = zeros(1,17);
    for j = 1:17
        [~, Atomic_M, ~] = PureComponentProperties(mixtures3, j);
        Atomic(j) = Atomic_M;
    end
    
    % Calculate Binary Mass Fraction based on Mole Fraction and Atomic Mass
    % Note: If sum is 0 (all zeros), handle division by zero
    numerator = molefraction_vec(moleFractions) .* Atomic;
    total_mass = sum(numerator);
    
    if total_mass == 0
        Binary_MassFraction = zeros(1,17);
    else
        Binary_MassFraction = numerator ./ total_mass;
    end

    %% 4. Prepare Neural Network Inputs
    T_eval = double(Temp);
    num_repeats = length(T_eval);
    
    % Create input matrix: [MoleFractions, MassFractions, Temperature]
    newfraction_variable = repmat([molefraction_vec(moleFractions) Binary_MassFraction], num_repeats, 1);
    final_data_mixture = [newfraction_variable T_eval];
    
    % Normalize Inputs
    inputSize = size(final_data_mixture, 2);
    final_data_mixture_normalized = (final_data_mixture' - inputMeans(1:inputSize,:)) ./ inputStds(1:inputSize,:);

    %% 5. Run Prediction (GPU/Deep Learning Toolbox)
    % Use dlarray as in original code
    prediction_normalized = double(extractdata(forward(netLR, dlarray(final_data_mixture_normalized, "CB"))));
    
    % Denormalize Output
    prediction_log = (prediction_normalized .* inputStds(end,:)) + inputMeans(end,:);
    raw_ann_prediction = exp(prediction_log)';

    %% 6. VFT Calculation (Optional)
    denom = 0;
    for j = 1:17
        Tg_k_j = Tg_fit_all_K(j);
        if Binary_MassFraction(j) > 0 % Avoid division issues if mass fraction is 0
             denom = denom + (Binary_MassFraction(j) / Tg_k_j);
        end
    end
    
    % Calculate Tg Mix
    if denom == 0
        Tg_mix_C = -273.15; % Fallback if denom is 0
    else
        Tg_mix_K = 1 / denom;
        Tg_mix_C = Tg_mix_K - 273.15;
    end
    
    % Apply VFT fit if requested
    if useVFT
        [mu_fit, ~, ~] = VFT2(T_eval, raw_ann_prediction, Tg_mix_C);
        prediction_mixture = mu_fit(:);
    else
        prediction_mixture = raw_ann_prediction(:);
    end
    
    Temp = Temp(:);
    prediction_mixture = prediction_mixture(:);

    %% 7. Construct Output Table
    T_table = table(Temp, prediction_mixture, 'VariableNames', {'Temperature','Viscosity'});
    
    for j = 1:17
        colname = matlab.lang.makeValidName(soluteNames{j});
        T_table.(sprintf('MF_%s', colname)) = repmat(moleFractions(j), length(Temp), 1);
    end
    
    for j = 1:17
        colname = matlab.lang.makeValidName(soluteNames{j}); 
        T_table.(sprintf('MassF_%s', colname)) = repmat(Binary_MassFraction(j), length(Temp), 1);
    end
    
    finalResults = T_table;

    %% 8. Export to Excel
    try
        writetable(finalResults, exportPath);
        fprintf('Successfully exported results to: %s\n', exportPath);
    catch ME
        warning('Failed to export Excel file. Check path permissions. Error: %s', ME.message);
    end

end

%% Helper Functions

function vec = molefraction_vec(in)
    % Ensures input is a row vector
    vec = in(:)';
end

function [Mol_V, Atomic_M, dens] = PureComponentProperties(mixtures, i)
    mixtureField = sprintf('Mixture_%d', i);
    Mol_V = mixtures.(mixtureField).MolarVolume;
    Atomic_M = mixtures.(mixtureField).AtomicMass;
    dens = mixtures.(mixtureField).Density;
end

function [mu_fit, xFit, r2] = VFT2(T_exp_C, mu_exp, iTg)
    T_exp_K = T_exp_C + 273.15;
    
    valid_indices = mu_exp > 0 & ~isnan(mu_exp) & ~isnan(T_exp_K);
    T_data = T_exp_K(valid_indices);
    mu_data = mu_exp(valid_indices);
    
    log_vft_fun = @(p, T) p(1) + p(2) ./ (T - p(3));
    
    log_mu_data = log(mu_data);
    
    iTg_K = iTg + 273.15;
    p0 = [-5, 1000, iTg_K - 50]; 
    
    min_T = min(T_data);
    lb = [-50,   0,  0];       
    ub = [ 20, 50000, min_T - 5]; 
    
    options = optimoptions('lsqcurvefit', 'Display', 'off', ...
                           'MaxFunctionEvaluations', 5000, ...
                           'FunctionTolerance', 1e-8, ...
                           'StepTolerance', 1e-8);
    
    try
        [pFit, ~] = lsqcurvefit(log_vft_fun, p0, T_data, log_mu_data, lb, ub, options);
        
        A_final = exp(pFit(1));
        B_final = pFit(2);
        T0_C_final = pFit(3) - 273.15;
        
        xFit = [A_final, B_final, T0_C_final];
        
        mu_fit = A_final .* exp(B_final ./ (T_exp_K - pFit(3)));
        
        log_mu_fit = log(mu_fit(valid_indices));
        SS_tot = sum((log_mu_data - mean(log_mu_data)).^2);
        SS_res = sum((log_mu_data - log_mu_fit).^2);
        
        r2 = 1 - (SS_res / SS_tot);
        
    catch ME
        warning('VFT Fit failed: %s', ME.message);
        xFit = [NaN, NaN, NaN];
        mu_fit = nan(size(T_exp_C));
        r2 = NaN;
    end
end