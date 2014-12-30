function [best_model_BIC,BIC_best,best_model_nll,nll_best] = master_function(start_vals,file_name,fid_BIC,fid_unix_params,fid_uniy_params,fid_GLC_params,fid_responses_accounted,fid_unix_learn_params,fid_uniy_learn_params)

% This function fits all the current models to a data set in the form
% data = [response x y].


global z_limit A_indices B_indices data
global fit_guessing fit_biased_guessing fit_unix fit_uniy fit_unequal_GCC
global fit_equal_GCC fit_GLC fit_GQC fit_GCC_guessing fit_GCC_eq2
global fit_unix_learn fit_uniy_learn

z_limit = 3;
noise = 5; 

% We print this value to 'fitting_results.dat'
fid2 = fopen([cd '/~output/' file_name '_results.dat'],'w');
fprintf(fid2, 'Fitting results for file: %s\n',file_name);
fprintf(fid2, 'Z limit: %10.5f\n\n', z_limit);

% We write the initial guesses to the file 'fitting_results.dat'
uni_xc = start_vals(1);
uni_yc = start_vals(2);
gcc_xc = start_vals(3);
gcc_yc = start_vals(4);
hybrid_center_point = start_vals(5:6);
hybrid_lin_point = start_vals(7:8);

fprintf(fid2, 'Starting Uni X Bound: %10.5f\n', uni_xc);
fprintf(fid2, 'Starting Uni Y Bound: %10.5f\n', uni_yc);
fprintf(fid2, 'Starting GCC X Bound: %10.5f\n', gcc_xc);
fprintf(fid2, 'Starting GCC Y Bound: %10.5f\n', gcc_yc);
fprintf(fid2, 'Starting Hybrid Center: %10.5f %10.5f\n', hybrid_center_point);
fprintf(fid2, 'Starting Hybrid Lin Point: %10.5f %10.5f\n', hybrid_lin_point);
fprintf(fid2, 'Starting Noise: %10.5f\n\n', noise);

% Write the number of data points to function_info
[nobs col] = size(data);

% We arrange so the bounds on the parameters vary according to the data.

x = data(:,2);
y = data(:,3);

xrange = max(x)-min(x);
xub = max(x)+0.1*xrange;
xlb = min(x)-0.1*xrange;

yrange = max(y)-min(y);
yub = max(y)+0.1*yrange;
ylb = min(y)-0.1*yrange;

bub = 2*max(yub,-ylb); % This is the upper bound on the b parameter in the linear models
blb = -bub;

nub = max(xrange,yrange)/2;
nlb = 0.001;

BIC_best = 1000;
best_model_BIC = fill_out('default',20);
nll_best = 1000;
best_model_nll = fill_out('default',20);

%% Guessing
if fit_guessing == 1

    fprintf(fid2,'Guessing Model:\n\n'); 

    negloglike_guess = -nobs*log(0.5);
    BIC_guess = 2*negloglike_guess;

    BIC_best = BIC_guess;
    best_model_BIC = fill_out('Guessing',20);
    nll_best = negloglike_guess;
    best_model_nll = fill_out('Guessing',20);

    % Report results
    fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_guess);
    fprintf(fid2,'BIC for Guessing Model:  %10.5f\n\n', BIC_guess);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_guess);

end

% Next we fit the biased guessing model.
if fit_biased_guessing == 1

    fprintf(fid2,'Biased Guessing Model:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) ~= 1);

    num_A = length(A_indices);
    num_B = length(B_indices);

    p1 = num_A/nobs;

    negloglike_biased1 = -num_A*log(p1)-num_B*log(1-p1);

    A_indices = find(data(:,1) ~= 1);
    B_indices = find(data(:,1) == 1);

    num_A = length(A_indices);
    num_B = length(B_indices);

    p2 = num_A/nobs;

    negloglike_biased2 = -num_A*log(p2)-num_B*log(1-p2);

    if (negloglike_biased1 <= negloglike_biased2) 

        p = p1;
        negloglike_biased = negloglike_biased1;
        switch_tag = 1;

    else

        p = p2;
        negloglike_biased = negloglike_biased2;
        switch_tag = 2;

    end

    BIC_biased = 2*negloglike_biased + log(nobs);

    if (BIC_biased < BIC_best)

        BIC_best = BIC_biased;
        best_model_BIC = fill_out('Biased Guessing',20);

    end

    if (negloglike_biased < nll_best)

        nll_best = negloglike_biased;
        best_model_nll = fill_out('Biased Guessing',20);

    end

    % Report results
    fprintf(fid2,'Switch tag: %i\n', switch_tag);
    fprintf(fid2,'p: %10.5f\n', p);
    fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_biased);
    fprintf(fid2,'BIC for Biased Guessing Model: %10.5f\n\n', BIC_biased);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_biased); 

end

%%
% Now we fit a whole series of models, starting with the unidimensional
% rules.

% Set search options
options = optimset(...
     'Display', 'off',...  % Spit out progress info
 'TolFun', .001);      % Termination tolerance on f

%% unix learning
%%%% Fitting Unidimensional X learning rule %%%%%
if fit_unix_learn == 1
    
    noise = 5;
    alpha_pos = 0.1;
    alpha_neg = 0.1;
    
    % We create our starting guess for the parameters.
    in_params = [noise alpha_pos alpha_neg];

    % Set contraints on parameters
    vlb = [xlb nlb];
    vub = [xub nub];

    % Run model fitting program
    fprintf('\nFitting Unidimensional X Learn\n\n')
    fprintf(fid2,'Unidimensional X Learn:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_unix1_learn] = fmincon('nll_unix_learn',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_unix2_learn] = fmincon('nll_unix_learn',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_unix1_learn <= negloglike_unix2_learn) 

        out_params = out_params1;
        negloglike_unix_learn = negloglike_unix1_learn;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_unix_learn = negloglike_unix2_learn;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=3;
    BIC_unix_learn=2*negloglike_unix_learn+(log(nobs)*fp);

    if (BIC_unix_learn < BIC_best)

        BIC_best = BIC_unix_learn;
        best_model_BIC = fill_out('Unidimensional X learn',20);

    end

    if (negloglike_unix_learn < nll_best)

        nll_best = negloglike_unix_learn;
        best_model_nll = fill_out('Unidimensional X learn',20);

    end

    % Report results
    fprintf(fid2,'Switch tag: %i\n', switch_tag);
    fprintf(fid2,'Noise: %10.5f\n', out_params(1));
    fprintf(fid2,'alpha pos: %10.5f\n', out_params(2));
    fprintf(fid2,'alpha neg: %10.5f\n', out_params(3));
    fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_unix_learn);
    fprintf(fid2,'BIC for Unidimensional X: %10.5f\n\n', BIC_unix_learn);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_unix_learn);

    % Print bound parameters to file
    fprintf(fid_unix_learn_params, '%10.5f %10.5f %10.5f %10.5f ', out_params);
    fprintf(fid_unix_learn_params, '\n');

    % compute percent results accounted for and print result to file
%     P_unix = compute_percent_accounted(data, [1 0 -out_params(1)]);
%     fprintf(fid_responses_accounted, '%10.5f %10.5f \t', P_unix);

end

%% uniy learning
%%%% Fitting Unidimensional Y learning rule %%%%%
if fit_uniy_learn == 1

    uni_xc = max(data(:,2)) - (max(data(:,2))-min(data(:,2)))/2; % Seems like these initial guesses could make a bug difference
    uni_yc = max(data(:,3)) - (max(data(:,3))-min(data(:,3)))/2;
    noise = 5;
    alpha_pos = 0.1;
    alpha_neg = 0.1;
    
    % We create our starting guess for the parameters.
    in_params = [uni_yc noise alpha_pos alpha_neg];

    % Set contraints on parameters
    vlb = [xlb nlb];
    vub = [xub nub];

    % Run model fitting program
    fprintf('\nFitting Unidimensional Y Learn\n\n')
    fprintf(fid2,'Unidimensional Y Learn:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_uniy1_learn] = fmincon('nll_uniy_learn',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_uniy2_learn] = fmincon('nll_uniy_learn',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_uniy1_learn <= negloglike_uniy2_learn) 

        out_params = out_params1;
        negloglike_uniy_learn = negloglike_uniy1_learn;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_uniy_learn = negloglike_uniy2_learn;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=3;
    BIC_uniy_learn=2*negloglike_uniy_learn+(log(nobs)*fp);

    if (BIC_uniy_learn < BIC_best)

        BIC_best = BIC_uniy_learn;
        best_model_BIC = fill_out('Unidimensional Y learn',20);

    end

    if (negloglike_uniy_learn < nll_best)

        nll_best = negloglike_uniy_learn;
        best_model_nll = fill_out('Unidimensional Y learn',20);

    end

    % Report results
    fprintf(fid2,'Switch tag: %i\n', switch_tag);
    fprintf(fid2,'final yc: %10.5f\n', out_params(1));
    fprintf(fid2,'Noise: %10.5f\n', out_params(2));
    fprintf(fid2,'alpha pos: %10.5f\n', out_params(3));
    fprintf(fid2,'alpha neg: %10.5f\n', out_params(4));
    fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_uniy_learn);
    fprintf(fid2,'BIC for Unidimensional Y: %10.5f\n\n', BIC_uniy_learn);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_uniy_learn);

    % Print bound parameters to file
    fprintf(fid_uniy_learn_params, '%10.5f %10.5f %10.5f %10.5f ', out_params);
    fprintf(fid_uniy_learn_params, '\n');

    % compute percent results accounted for and print result to file
%     P_unix = compute_percent_accounted(data, [1 0 -out_params(1)]);
%     fprintf(fid_responses_accounted, '%10.5f %10.5f \t', P_unix);

end

%% unix
%%%% Fitting Unidimensional X rule %%%%%
if fit_unix == 1

    % We create our starting guess for the parameters.
    in_params = [uni_xc noise];

    % Set contraints on parameters
    vlb = [xlb nlb]
    vub = [xub nub]

    % Run model fitting program
    fprintf('\nFitting Unidimensional X\n\n')
    fprintf(fid2,'Unidimensional X:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_unix1] = fmincon('nll_unix_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_unix2] = fmincon('nll_unix_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_unix1 <= negloglike_unix2) 

        out_params = out_params1;
        negloglike_unix = negloglike_unix1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_unix = negloglike_unix2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=2;
    BIC_unix=2*negloglike_unix+(log(nobs)*fp);

    if (BIC_unix < BIC_best)

        BIC_best = BIC_unix;
        best_model_BIC = fill_out('Unidimensional X',20);

    end

    if (negloglike_unix < nll_best)

        nll_best = negloglike_unix;
        best_model_nll = fill_out('Unidimensional X',20);

    end


    % Report results
    fprintf(fid2,'Switch tag: %i\n', switch_tag);
    fprintf(fid2,'X Bound: %10.5f\n', out_params(1));
    fprintf(fid2,'Noise: %10.5f\n', out_params(2));
    fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_unix);
    fprintf(fid2,'BIC for Unidimensional X: %10.5f\n\n', BIC_unix);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_unix);

    % Print bound parameters to file
    fprintf(fid_unix_params, '%10.5f ', [1 0 -out_params(1)]);
    fprintf(fid_unix_params, '\n');

    % compute percent results accounted for and print result to file
    P_unix = compute_percent_accounted(data, [1 0 -out_params(1)]);
    fprintf(fid_responses_accounted, '%10.5f %10.5f \t', P_unix);

end

%% uniy
%%%% Fitting Unidimensional Y rule %%%%%
if fit_uniy == 1

    % We create our starting guess for the parameters.
    in_params = [uni_yc noise];

    % Set contraints on parameters
    vlb = [ylb nlb];
    vub = [yub nub];

    % Run model fitting program
    fprintf('\n\nFitting Unidimensional Y\n\n')
    fprintf(fid2,'Unidimensional Y:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_uniy1] = fmincon('nll_uniy_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_uniy2] = fmincon('nll_uniy_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_uniy1 <= negloglike_uniy2) 

        out_params = out_params1;
        negloglike_uniy = negloglike_uniy1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_uniy = negloglike_uniy2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=2;
    BIC_uniy=2*negloglike_uniy+(log(nobs)*fp);

    if (BIC_uniy < BIC_best)

        BIC_best = BIC_uniy;
        best_model_BIC = fill_out('Unidimensional Y',20);

    end

    if (negloglike_uniy < nll_best)

        nll_best = negloglike_uniy;
        best_model_nll = fill_out('Unidimensional Y',20);

    end

    % Report results
    fprintf(fid2,'Switch tag: %i\n', switch_tag);
    fprintf(fid2,'Y Bound: %10.5f\n', out_params(1));
    fprintf(fid2,'Noise: %10.5f\n', out_params(2));
    fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_uniy);
    fprintf(fid2,'BIC for Unidimensional Y: %10.5f\n\n', BIC_uniy);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_uniy);

    % Print bound parameters to file
    fprintf(fid_uniy_params, '%10.5f ', [0 1 -out_params(1)]);
    fprintf(fid_uniy_params, '\n');

    % compute percent results accounted for and print result to file
    P_uniy = compute_percent_accounted(data, [0 1 -out_params(1)]);
    fprintf(fid_responses_accounted, '%10.5f %10.5f \t', P_uniy);

end

%% GCC Unequal
%%%%%%  Fitting Unequal Variance Conjunctive Rule %%%%%%
if fit_unequal_GCC == 1

    % We create our starting guess for the parameters.
    in_params = [gcc_xc gcc_yc noise noise];

    % Set contraints on parameters
    vlb = [xlb ylb nlb nlb];
    vub = [xub yub nub nub];

    % Run model fitting program
    fprintf('\n\nFitting Unequal GCC\n\n')
    fprintf(fid2,'Unequal GCC:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) ~= 1);

    [out_params1,negloglike_gcc_uneq1] = fmincon('nll_gcc_uneq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_gcc_uneq2] = fmincon('nll_gcc_uneq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_gcc_uneq1 <= negloglike_gcc_uneq2) 

        out_params = out_params1;
        negloglike_gcc_uneq = negloglike_gcc_uneq1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_gcc_uneq = negloglike_gcc_uneq2;
        switch_tag = 2;    

    end

    % Compute BIC score
    fp=4;
    BIC_gcc_uneq=2*negloglike_gcc_uneq+(log(nobs)*fp);

    if (BIC_gcc_uneq < BIC_best)

        BIC_best = BIC_gcc_uneq;
        best_model_BIC = fill_out('Unequal GCC',20);

    end

    if (negloglike_gcc_uneq < nll_best)

        nll_best = negloglike_gcc_uneq;
        best_model_nll = fill_out('Unequal GCC',20);

    end

    % Report results
    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
    fprintf(fid2, 'X Bound: %10.5f\n', out_params(1));
    fprintf(fid2, 'Y Bound: %10.5f\n', out_params(2));
    fprintf(fid2, 'Noise X: %10.5f\n', out_params(3));
    fprintf(fid2, 'Noise Y: %10.5f\n', out_params(4));
    fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_gcc_uneq);
    fprintf(fid2, 'BIC for GCC (Unequal Variance Case): %10.5f\n\n', BIC_gcc_uneq);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_gcc_uneq);

end

%% GCC Equal
%%%%%%  Fitting Equal Variance Conjunctive Rule %%%%%%
if fit_equal_GCC == 1

    % We create our starting guess for the parameters.
    in_params = [gcc_xc gcc_yc noise];

    % Set contraints on parameters
    vlb = [xlb ylb nlb];
    vub = [xub yub nub];

    % Run model fitting program
    fprintf('\n\nFitting Equal GCC\n\n')
    fprintf(fid2,'Equal GCC:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_gcc_eq1] = fmincon('nll_gcc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_gcc_eq2] = fmincon('nll_gcc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_gcc_eq1 <= negloglike_gcc_eq2) 

        out_params = out_params1;
        negloglike_gcc_eq = negloglike_gcc_eq1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_gcc_eq = negloglike_gcc_eq2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=3;
    BIC_gcc_eq=2*negloglike_gcc_eq+(log(nobs)*fp);

    if (BIC_gcc_eq < BIC_best)

        BIC_best = BIC_gcc_eq;
        best_model_BIC = fill_out('Equal GCC',20);

    end

    if (negloglike_gcc_eq < nll_best)

        nll_best = negloglike_gcc_eq;
        best_model_nll = fill_out('Equal GCC',20);

    end

    % Report results
    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
    fprintf(fid2, 'X Bound: %10.5f\n', out_params(1));
    fprintf(fid2, 'Y Bound: %10.5f\n', out_params(2));
    fprintf(fid2, 'Noise: %10.5f\n', out_params(3));
    fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_gcc_eq);
    fprintf(fid2, 'BIC for GCC (Equal Variance Case): %10.5f\n\n', BIC_gcc_eq);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_gcc_eq);

end

%% GCC Guess
%%%%%%% Fitting GCC guessing model %%%%%%%%%%
if fit_GCC_guessing == 1

    % We create our starting guess for the parameters.
    in_params = [gcc_xc gcc_yc noise];

    % Set contraints on parameters
    vlb = [xlb ylb nlb];
    vub = [xub yub nub];

    % Run model fitting program
    fprintf('\n\nFitting GCC guessing\n\n')
    fprintf(fid2,'GCC guessing:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_gcc_guessing1] = fmincon('nll_gcc_guess_vMJC',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_gcc_guessing2] = fmincon('nll_gcc_guess_vMJC',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_gcc_guessing1 <= negloglike_gcc_guessing2) 

        out_params = out_params1;
        negloglike_gcc_guessing = negloglike_gcc_guessing1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_gcc_guessing = negloglike_gcc_guessing2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=3;
    BIC_gcc_guessing=2*negloglike_gcc_guessing+(log(nobs)*fp);

    if (BIC_gcc_guessing < BIC_best)

        BIC_best = BIC_gcc_guessing;
        best_model_BIC = fill_out('GCC Guessing',20);

    end

    if (negloglike_gcc_guessing < nll_best)

        nll_best = negloglike_gcc_guessing;
        best_model_nll = fill_out('GCC Guessing',20);

    end

    % Report results
    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
    fprintf(fid2, 'X Bound: %10.5f\n', out_params(1));
    fprintf(fid2, 'Y Bound: %10.5f\n', out_params(2));
    fprintf(fid2, 'Noise: %10.5f\n', out_params(3));
    fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_gcc_guessing);
    fprintf(fid2, 'BIC for GCC guessing: %10.5f\n\n', BIC_gcc_guessing);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_gcc_guessing);

end

%% GCC Equal 2
%%%%%%% Fitting GCC equal variance 2 model %%%%%%%%%%

if fit_GCC_eq2 == 1

    % We create our starting guess for the parameters.
    in_params = [gcc_xc gcc_yc noise];

    % Set contraints on parameters
    vlb = [xlb ylb nlb];
    vub = [xub yub nub];

    % Run model fitting program
    fprintf('\n\nFitting Equal Variance GCC 2\n\n')
    fprintf(fid2,'GCC eq2:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_gcc_eq2_1] = fmincon('nll_gcc_eq2_vMJC',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_gcc_eq2_2] = fmincon('nll_gcc_eq2_vMJC',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_gcc_eq1 <= negloglike_gcc_eq2) 

        out_params = out_params1;
        negloglike_gcc_eq2 = negloglike_gcc_eq2_1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_gcc_eq2 = negloglike_gcc_eq2_2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=3;
    BIC_gcc_eq2=2*negloglike_gcc_eq2+(log(nobs)*fp);

    if (BIC_gcc_eq2 < BIC_best)

        BIC_best = BIC_gcc_eq2;
        best_model_BIC = fill_out('Equal GCC 2',20);

    end

    if (negloglike_gcc_eq2 < nll_best)

        nll_best = negloglike_gcc_eq2;
        best_model_nll = fill_out('Equal GCC 2',20);

    end

    % Report results
    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
    fprintf(fid2, 'X Bound: %10.5f\n', out_params(1));
    fprintf(fid2, 'Y Bound: %10.5f\n', out_params(2));
    fprintf(fid2, 'Noise: %10.5f\n', out_params(3));
    fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_gcc_eq2);
    fprintf(fid2, 'BIC for Equal GCC 2: %10.5f\n\n', BIC_gcc_eq2);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_gcc_eq2);

end

%% GLC
%%%%% Fitting Equal Variance GLC %%%%%%
if fit_GLC == 1

    % Determine the best fitting linear bound based on the data to be analyzed (either response or category depending on use_resp)
    lin_bound=fisherdiscrim2d([data(:,1:3) data(:,1)],5);

    % The definition of in_params depends on a1^2 + a2^2 = 1 and on a2 >= 0, so
    % we normalize just to be safe.
    a1 = lin_bound(1)/sqrt(lin_bound(1)^2+lin_bound(2)^2);
    b = lin_bound(3)/sqrt(lin_bound(1)^2+lin_bound(2)^2);

    if (lin_bound(2) < 0) 

        a1 = -a1;
        b = -b;

    end

    in_params = [a1 b noise];

    % In addition, we redefine data for the sake of nll_glc_eq_vJ.  
    global data_info

    data_info = data;
    data_info(:,4) = ones(length(data(:,1)),1);

    % Set contraints on parameters (Note that a1 varies between -1 and 1)
    vlb = [-1 blb nlb];
    vub = [1 bub nub];

    % Run model fitting program
    fprintf('\n\nFitting Equal GLC\n'); 
    fprintf(fid2,'Equal GLC:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_glc_eq1] = fmincon('nll_glc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_glc_eq2] = fmincon('nll_glc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_glc_eq1 <= negloglike_glc_eq2) 

        out_params = out_params1;
        negloglike_glc_eq = negloglike_glc_eq1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_glc_eq = negloglike_glc_eq2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=3;
    BIC_glc_eq=2*negloglike_glc_eq+(log(nobs)*fp);

    if (BIC_glc_eq < BIC_best)

        BIC_best = BIC_glc_eq;
        best_model_BIC = fill_out('Equal GLC',20);

    end

    if (negloglike_glc_eq < nll_best)

        nll_best = negloglike_glc_eq;
        best_model_nll = fill_out('Equal GLC',20);

    end

    % Report results
    a1 = out_params(1);
    a2 = sqrt(1-a1^2);
    b = out_params(2);

    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
    fprintf(fid2, 'a1: %10.5f\n', a1);
    fprintf(fid2, 'a2: %10.5f\n', a2);
    fprintf(fid2, 'b: %10.5f\n', b);
    fprintf(fid2, 'Noise: %10.5f\n', out_params(3));
    fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_glc_eq);
    fprintf(fid2, 'BIC for GLC (Equal Variance Case): %10.5f\n\n', BIC_glc_eq);

    % Print results to global BIC file
    fprintf(fid_BIC,'%10.5f\t', BIC_glc_eq);

    % Print bound parameters to file
    fprintf(fid_GLC_params, '%10.5f ', [a1 a2 b]);
    fprintf(fid_GLC_params, '\n');

    % compute percent results accounted for and print result to file
    P_glc = compute_percent_accounted(data, [a1 a2 b]);
    fprintf(fid_responses_accounted, '%10.5f %10.5f\t', P_glc);

end

%% GQC
%%%%%%% Fitting GQC model %%%%%%%%%%

if fit_GQC == 1
  
    % Determine the best fitting linear bound based on the data to be analyzed (either response or category depending on use_resp)
    lin_bound=fisherdiscrim2d([data(:,1:3) data(:,1)],5);

    % The definition of in_params depends on b1^2 + b2^2 = 1 and on b2 >= 0, so
    % we normalize just to be safe.
    b1 = lin_bound(1)/sqrt(lin_bound(1)^2+lin_bound(2)^2);
    c = lin_bound(3)/sqrt(lin_bound(1)^2+lin_bound(2)^2);

    if (lin_bound(2) < 0) 

        b1 = -b1;
        c = -c;

    end

    % We define our initial parameters
    in_params = [0 0 0 b1 c noise noise];

    % Set contraints on parameters (Note that b1 varies between -1 and 1)
    vlb = [blb blb blb -1 blb nlb nlb];
    vub = [bub bub bub 1 bub nub nub];

    fprintf('\n\nFitting GQC\n'); 
    fprintf(fid2,'GQC:\n\n'); 

    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) == 2);

    [out_params1,negloglike_gqc_1] = fmincon('nll_gqc_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    A_indices = find(data(:,1) == 2);
    B_indices = find(data(:,1) == 1);

    [out_params2,negloglike_gqc_2] = fmincon('nll_gqc_vJ',in_params,[],[],[],[],vlb,vub,[],options);

    if (negloglike_gqc_1 <= negloglike_gqc_2) 

        out_params = out_params1;
        negloglike_gqc = negloglike_gqc_1;
        switch_tag = 1;

    else

        out_params = out_params2;
        negloglike_gqc = negloglike_gqc_2;
        switch_tag = 2;

    end

    % Compute BIC score
    fp=5;
    BIC_gqc = 2*negloglike_gqc+(log(nobs)*fp);

    if (BIC_gqc < BIC_best)

        BIC_best = BIC_gqc;
        best_model_BIC = fill_out('GQC',20);

    end

    if (negloglike_gqc < nll_best)

        nll_best = negloglike_gqc;
        best_model_nll = fill_out('GQC',20);

    end

    % Report results
    a1 = out_params(1);
    a2 = out_params(2);
    a3 = out_params(3);

    b1 = out_params(4);
    b2 = sqrt(1-b1^2);

    c = out_params(5);

    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
    fprintf(fid2, 'a1: %10.5f\n', a1);
    fprintf(fid2, 'a2: %10.5f\n', a2);
    fprintf(fid2, 'a3: %10.5f\n', a3);
    fprintf(fid2, 'b1: %10.5f\n', b1);
    fprintf(fid2, 'b2: %10.5f\n', b2);
    fprintf(fid2, 'c: %10.5f\n', c);
    fprintf(fid2, 'Perceptual Noise: %10.5f\n', out_params(6));
    fprintf(fid2, 'Criterial Noise: %10.5f\n', out_params(7));
    fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_gqc);
    fprintf(fid2, 'BIC for GQC: %10.5f\n\n', BIC_gqc);
    fprintf(fid_BIC,'%10.5f\t', BIC_gqc);

end

%%
fprintf(fid_responses_accounted, '\n');
fprintf(fid_BIC,'\n');
fprintf(fid2, 'Best fitting model by BIC:\n');
fprintf(fid2, '%s\n', best_model_BIC);
fprintf(fid2, 'BIC for that model:\n');
fprintf(fid2, '%10.5f\n', BIC_best);
fprintf(fid2, 'Best fitting model by Negative Log Likelihood:\n');
fprintf(fid2, '%s\n', best_model_nll);
fprintf(fid2, 'Negative Log Likelihood for that model:\n');
fprintf(fid2, '%10.5f\n', nll_best);
    
fclose(fid2);