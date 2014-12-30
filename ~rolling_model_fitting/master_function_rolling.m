function [best_model_BIC,BIC_best,best_model_nll,nll_best,out_vals] = master_function_rolling(start_vals,file_name,fit_hybrid)

% [best_model_BIC,BIC_best,best_model_nll,nll_best] = master_function(start_vals,file_name,fit_hybrid)

% This function fits all the current models to a data set in the form
% data = [response x y].  Also:
% 
% start_vals = [uni_xc uni_yc gcc_xc gcc_yc hybrid_center_point hybrid_lin_point] contains the starting
% bounds for the model fitting

% The user has a choice whether or not to fit one or both of the hybrid models.  If
% fit_hybrid = 0, neither of the hybrids are fit (default).  If fit_hybrid = 1, only 
% the equal variance hybrid model is fit.  If fit_hybrid = 2 then both
% hybrid models are fit.

if (nargin < 3) 
    fit_hybrid = 0;
end

global z_limit A_indices B_indices data

z_limit = 3;
noise = 5; 

% We print this value to 'fitting_results.dat'

fid2 = fopen([cd '\output\' file_name '_results.dat'],'w');
fprintf(fid2, 'Fitting results for file: %s\n',file_name);
fprintf(fid2, 'Z limit: %10.5f\n\n', z_limit);

% We write the initial guesses to the file 'fitting_results.dat'

uni_xc = start_vals(1);
uni_yc = start_vals(2);
gcc_xc = start_vals(3);
gcc_yc = start_vals(4);
hybrid_center_point = start_vals(5:6);
hybrid_lin_point = start_vals(7:8);

out_vals(1:8) = 0;

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

% We now fit a large number of models.  Let's start by fitting the guessing
% models.

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

% Next we fit the biased guessing model.

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

% Now we fit a whole series of models, starting with the unidimensional
% rules.

% Set search options

options = optimset(...
     'Display', 'off',...  % Spit out progress info
 'TolFun', .001);      % Termination tolerance on f    

%%%% Fitting Unidimensional X rule %%%%%

% We create our starting guess for the parameters.

in_params = [uni_xc noise];

% Set contraints on parameters

vlb = [xlb nlb];
vub = [xub nub];

% Run model fitting program

fprintf('\nFitting Unidimensional X\n\n')
fprintf(fid2,'Unidimensional X:\n\n'); 

A_indices = find(data(:,1) == 1);
B_indices = find(data(:,1) ~= 1);

[out_params1,negloglike_unix1] = fmincon('nll_unix_vJ',in_params,[],[],[],[],vlb,vub,[],options);

A_indices = find(data(:,1) ~= 1);
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

out_vals(1) = out_params(1); % Record best unidimensional x bound for next fit

% Report results

fprintf(fid2,'Switch tag: %i\n', switch_tag);
fprintf(fid2,'X Bound: %10.5f\n', out_params(1));
fprintf(fid2,'Noise: %10.5f\n', out_params(2));
fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_unix);
fprintf(fid2,'BIC for Unidimensional X: %10.5f\n\n', BIC_unix);

%%%% Fitting Unidimensional Y rule %%%%%

% We create our starting guess for the parameters.

in_params = [uni_yc noise];

% Set contraints on parameters

vlb = [ylb nlb];
vub = [yub nub];

% Run model fitting program
fprintf('\n\nFitting Unidimensional Y\n\n')
fprintf(fid2,'Unidimensional Y:\n\n'); 

A_indices = find(data(:,1) == 1);
B_indices = find(data(:,1) ~= 1);

[out_params1,negloglike_uniy1] = fmincon('nll_uniy_vJ',in_params,[],[],[],[],vlb,vub,[],options);

A_indices = find(data(:,1) ~= 1);
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

out_vals(2) = out_params(1); % Record best unidimensional y bound for next fit

% Report results

fprintf(fid2,'Switch tag: %i\n', switch_tag);
fprintf(fid2,'Y Bound: %10.5f\n', out_params(1));
fprintf(fid2,'Noise: %10.5f\n', out_params(2));
fprintf(fid2,'Negative Log Likelihood: %10.5f\n', negloglike_uniy);
fprintf(fid2,'BIC for Unidimensional Y: %10.5f\n\n', BIC_uniy);

% Next we fit the conjunctive rules.

%%%%%%  Fitting Unequal Variance Conjunctive Rule %%%%%%

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

A_indices = find(data(:,1) ~= 1);
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

%%%%%%  Fitting Equal Variance Conjunctive Rule %%%%%%

% We create our starting guess for the parameters.

in_params = [gcc_xc gcc_yc noise];

% Set contraints on parameters

vlb = [xlb ylb nlb];
vub = [xub yub nub];

% Run model fitting program

fprintf('\n\nFitting Equal GCC\n\n')
fprintf(fid2,'Equal GCC:\n\n'); 

A_indices = find(data(:,1) == 1);
B_indices = find(data(:,1) ~= 1);

[out_params1,negloglike_gcc_eq1] = fmincon('nll_gcc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

A_indices = find(data(:,1) ~= 1);
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

gcc_switch_tag = switch_tag;

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

out_vals(3) = out_params(1); % Record best gcc x bound for next fit
out_vals(4) = out_params(2); % Record best gcc y bound for next fit

% Report results

fprintf(fid2, 'Switch tag: %i\n', switch_tag);
fprintf(fid2, 'X Bound: %10.5f\n', out_params(1));
fprintf(fid2, 'Y Bound: %10.5f\n', out_params(2));
fprintf(fid2, 'Noise: %10.5f\n', out_params(3));
fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_gcc_eq);
fprintf(fid2, 'BIC for GCC (Equal Variance Case): %10.5f\n\n', BIC_gcc_eq);

% Next we fit the GLC model.

%%%%% Fitting Equal Variance GLC %%%%%%

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
B_indices = find(data(:,1) ~= 1);

[out_params1,negloglike_glc_eq1] = fmincon('nll_glc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

A_indices = find(data(:,1) ~= 1);
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

% Now we fit the GQC model.

%%%%%%% Fitting GQC model %%%%%%%%%%
  
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

% Run model fitting program

fprintf('\n\nFitting GQC\n'); 
fprintf(fid2,'GQC:\n\n'); 

A_indices = find(data(:,1) == 1);
B_indices = find(data(:,1) ~= 1);

[out_params1,negloglike_gqc_1] = fmincon('nll_gqc_vJ',in_params,[],[],[],[],vlb,vub,[],options);

A_indices = find(data(:,1) ~= 1);
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

if (fit_hybrid == 1)|(fit_hybrid == 2)

	% Finally we fit the hybrid models
	
	%%%%% Fitting Unequal Variance Hybrid %%%%%%  % We are currently not
	%%%%% fitting the unequal variance hybrid, but we still need the
	%%%%% introductory lines of codes so that the equal variance hybrid will
	%%%%% fit properly.
	
	% We declare how large we'd like the z-grid (maximum width = 1000), and the
	% smallest step size used to compute numerical derivatives.  We write these
	% values to fitting_results.dat.
	
	grid_width = 100;
	dmc = 0.005;
	
	fprintf(fid2, 'Grid Width: %i\n', grid_width);
	fprintf(fid2, 'DiffMinChange: %10.5f\n\n', dmc);
	
	% We create a starting guess for the linear bound using the points in
	% start_vals.  
	
	[a1 a2 b] = line_from_points(hybrid_center_point,hybrid_lin_point);
	
	in_params = [a1 b hybrid_center_point(1) noise noise];
	
	% Set contraints on parameters (Note that a1 is constrained between -1 and 1)
	
	vlb = [-1 blb xlb nlb nlb];
	vub = [1 bub xub nub nub];
	
	% We set the search options.
	
	options = optimset(...
         'Display', 'off',...  % Spit out progress info
     'DiffMinChange', dmc,...    
     'TolFun', .001);      % Termination tolerance on f    
	
    % We run the model fitting based on the gcc switch tag.
    
    switch_tag = gcc_switch_tag;
 
     if (fit_hybrid == 2)
		% Run model fitting program
		
		fprintf('\n\nFitting Unequal Hybrid\n'); 
		fprintf(fid2,'Unequal Hybrid:\n\n'); 
		
		% The hybrid model fitting programs access the function_info.dat
		% file to retreive the desired size of the z-grid, the number of
		% observations, the indices of the A and B responses, and the coordinates of the data.
		
		fid=fopen('function_info.dat','w');
		fprintf(fid, '%i\n', grid_width);

        % We run the model fitting based on the gcc switch tag.
        
        if (switch_tag == 1)
            
    		A_indices = find(data(:,1) == 1);
			B_indices = find(data(:,1) ~= 1);
    
        else
            
            A_indices = find(data(:,1) ~= 1);
			B_indices = find(data(:,1) == 1);
	
        end
        
		fprintf(fid, '%i %i %i\n',nobs, num_A, num_B);
		
		% Write the A and B indices to function_info
		
		A_indices_long = A_indices;
		B_indices_long = B_indices;
		
		A_indices_long((num_A+1:1000))=0;
		B_indices_long((num_B+1:1000))=0;
		
		fprintf(fid, '%i ',A_indices_long);
		fprintf(fid, '\n');
		fprintf(fid, '%i ',B_indices_long);
		fprintf(fid, '\n');
		
		% Write the coordinates of the stimuli to function_info
		
		fprintf(fid, '%6.2f %6.2f\n', data(:,2:3)');
		fclose(fid);
		
		[out_params,negloglike_hybrid_uneq] = fmincon('nll_hybrid_uneq_vJ',in_params,[],[],[],[],vlb,vub,[],options);
				
		% Compute BIC score
		
		fp=5;
		BIC_hybrid_uneq=2*negloglike_hybrid_uneq+(log(nobs)*fp);
		
		if (BIC_hybrid_uneq < BIC_best)
		    
		    BIC_best = BIC_hybrid_uneq;
		    best_model_BIC = fill_out('Unequal Hybrid',20);
		    
		end
		
		if (negloglike_hybrid_uneq < nll_best)
		    
		    nll_best = negloglike_hybrid_uneq;
		    best_model_nll = fill_out('Unequal Hybrid',20);
		    
		end
			
		% Report results
		
		a1 = out_params(1);
		a2 = sqrt(1-a1^2);
		b = out_params(2);
		
		fprintf(fid2, 'Switch tag: %i\n', switch_tag);
		fprintf(fid2, 'a1: %10.5f\n', a1);
		fprintf(fid2, 'a2: %10.5f\n', a2);
		fprintf(fid2, 'b: %10.5f\n', b);
		fprintf(fid2, 'Unidimensional X Bound: %10.5f\n', out_params(3));
		fprintf(fid2, 'Noise X: %10.5f\n', out_params(4));
		fprintf(fid2, 'Noise Y: %10.5f\n', out_params(5));
		fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_hybrid_uneq);
		fprintf(fid2, 'BIC for Hybrid (Unequal Variance Case): %10.5f\n\n', BIC_hybrid_uneq);
		
	end
    
	% %%%%% Fitting Equal Variance Hybrid %%%%%%
	
	% We use the same starting parameters as in the unequal variance
	% hybrid, but with one fewer noise parameter.
	
	in_params = in_params(1:4);
	
	% Set contraints on parameters (Note that a1 is constrained between -1 and 1)
	
	vlb = [-1 blb xlb nlb];
	vub = [1 bub xub nub];
	
	% Run model fitting program
	
	fprintf('\n\nFitting Equal Hybrid\n'); 
	fprintf(fid2,'Equal Hybrid:\n\n'); 
	
	% The hybrid model fitting programs access the function_info.dat
	% file to retreive the desired size of the z-grid, the number of
	% observations, the indices of the A and B responses, and the coordinates of the data.
	
	fid=fopen('function_info.dat','w');
	fprintf(fid, '%i\n', grid_width);

    % We run the model fitting based on the gcc switch tag.
    
    if (switch_tag == 1)
        
		A_indices = find(data(:,1) == 1);
		B_indices = find(data(:,1) ~= 1);

    else
        
        A_indices = find(data(:,1) ~= 1);
		B_indices = find(data(:,1) == 1);

    end
       
	num_A = length(A_indices);
	num_B = length(B_indices);
	
	fprintf(fid, '%i %i %i\n',nobs, num_A, num_B);
	
	% Write the A and B indices to function_info
	
	A_indices_long = A_indices;
	B_indices_long = B_indices;
	
	A_indices_long((num_A+1:1000))=0;
	B_indices_long((num_B+1:1000))=0;
	
	fprintf(fid, '%i ',A_indices_long);
	fprintf(fid, '\n');
	fprintf(fid, '%i ',B_indices_long);
	fprintf(fid, '\n');
	
	% Write the coordinates of the stimuli to function_info
	
	fprintf(fid, '%6.2f %6.2f\n', data(:,2:3)');
	fclose(fid);
	
	[out_params,negloglike_hybrid_eq] = fmincon('nll_hybrid_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);
	
	% Compute BIC score
	
	fp=4;
	BIC_hybrid_eq=2*negloglike_hybrid_eq+(log(nobs)*fp);
	
	if (BIC_hybrid_eq < BIC_best)
        
        BIC_best = BIC_hybrid_eq;
        best_model_BIC = fill_out('Equal Hybrid',20);
        
	end
	
	if (negloglike_hybrid_eq < nll_best)
        
        nll_best = negloglike_hybrid_eq;
        best_model_nll = fill_out('Equal Hybrid',20);
        
	end
	
	% Report results
	
	a1 = out_params(1);
	a2 = sqrt(1-a1^2);
	b = out_params(2);

    if (a2 ~= 0) 
        
        out_vals(5) = out_params(3);
        out_vals(6) = -(a1*out_params(3)+b)/a2;
        out_vals(7) = 0;
        out_vals(8) = -b/a2;

    else
        
        out_vals(5) = out_params(3);
        out_vals(6) = 50;
        out_vals(7) = out_params(3);
        out_vals(8) = 0;
        
    end
        
    fprintf(fid2, 'Switch tag: %i\n', switch_tag);
	fprintf(fid2, 'a1: %10.5f\n', a1);
	fprintf(fid2, 'a2: %10.5f\n', a2);
	fprintf(fid2, 'b: %10.5f\n', b);
	fprintf(fid2, 'Unidimensional X Bound: %10.5f\n', out_params(3));
	fprintf(fid2, 'Noise: %10.5f\n', out_params(4));
	fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_hybrid_eq);
	fprintf(fid2, 'BIC for Hybrid (Equal Variance Case): %10.5f\n\n', BIC_hybrid_eq);
    
  	% %%%%% Fitting Lin Bound Only Hybrid %%%%%%
	
    % We use the same starting parameters as in the equal hybrid, but this time with no
    % unidimensional x bound as it is fixed to be 50
    
   	in_params = [in_params(1:2) in_params(4)];
	
	% Set contraints on parameters (Note that a1 is constrained between -1 and 1)
	
	vlb = [-1 blb nlb];
	vub = [1 bub nub];
	
	% Run model fitting program
	
	fprintf('\n\nFitting Lin Only Hybrid\n'); 
	fprintf(fid2,'Lin Only Hybrid:\n\n'); 
	
	% The hybrid model fitting programs access the function_info.dat
	% file to retreive the desired size of the z-grid, the number of
	% observations, the indices of the A and B responses, and the coordinates of the data.
	
	fid=fopen('function_info.dat','w');
	fprintf(fid, '%i\n', grid_width);

    % We run the model fitting based on the gcc switch tag.
    
    if (switch_tag == 1)
        
		A_indices = find(data(:,1) == 1);
		B_indices = find(data(:,1) ~= 1);

    else
        
        A_indices = find(data(:,1) ~= 1);
		B_indices = find(data(:,1) == 1);

    end
       
	num_A = length(A_indices);
	num_B = length(B_indices);
	
	fprintf(fid, '%i %i %i\n',nobs, num_A, num_B);
	
	% Write the A and B indices to function_info
	
	A_indices_long = A_indices;
	B_indices_long = B_indices;
	
	A_indices_long((num_A+1:1000))=0;
	B_indices_long((num_B+1:1000))=0;
	
	fprintf(fid, '%i ',A_indices_long);
	fprintf(fid, '\n');
	fprintf(fid, '%i ',B_indices_long);
	fprintf(fid, '\n');
	
	% Write the coordinates of the stimuli to function_info
	
	fprintf(fid, '%6.2f %6.2f\n', data(:,2:3)');
	fclose(fid);
	
	[out_params,negloglike_hybrid_lin_only] = fmincon('nll_hybrid_lin_only_vJ',in_params,[],[],[],[],vlb,vub,[],options);
	
	% Compute BIC score
	
	fp=3;
	BIC_hybrid_lin_only=2*negloglike_hybrid_lin_only+(log(nobs)*fp);
	
	if (BIC_hybrid_lin_only < BIC_best)
        
        BIC_best = BIC_hybrid_lin_only;
        best_model_BIC = fill_out('Lin Only Hybrid',20);
        
	end
	
	if (negloglike_hybrid_lin_only < nll_best)
        
        nll_best = negloglike_hybrid_lin_only;
        best_model_nll = fill_out('Lin Only Hybrid',20);
        
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
	fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_hybrid_lin_only);
	fprintf(fid2, 'BIC for Hybrid (Lin Only Case): %10.5f\n\n', BIC_hybrid_lin_only);

    % %%%%% Fitting Optimal Hybrid %%%%%%
	
	% We use the same starting parameters as in the unequal variance hybrid 
	
	in_params = [in_params(3)];
	
	% Set contraints on parameters (Note that a1 is constrained between -1 and 1)
	
	vlb = [nlb];
	vub = [nub];
	
	% Run model fitting program
	
	fprintf('\n\nFitting Optimal Hybrid\n'); 
	fprintf(fid2,'Optimal Hybrid:\n\n'); 
	
	% The hybrid model fitting programs access the function_info.dat
	% file to retreive the desired size of the z-grid, the number of
	% observations, the indices of the A and B responses, and the coordinates of the data.
	
	fid=fopen('function_info.dat','w');
	fprintf(fid, '%i\n', grid_width);
	
    % We run the model fitting based on the gcc switch tag.
    
    if (switch_tag == 1)
        
		A_indices = find(data(:,1) == 1);
		B_indices = find(data(:,1) ~= 1);
	
	else
        
        A_indices = find(data(:,1) ~= 1);
		B_indices = find(data(:,1) == 1);
	
	end
      
	num_A = length(A_indices);
	num_B = length(B_indices);
	
	fprintf(fid, '%i %i %i\n',nobs, num_A, num_B);
	
	% Write the A and B indices to function_info
	
	A_indices_long = A_indices;
	B_indices_long = B_indices;
	
	A_indices_long((num_A+1:1000))=0;
	B_indices_long((num_B+1:1000))=0;
	
	fprintf(fid, '%i ',A_indices_long);
	fprintf(fid, '\n');
	fprintf(fid, '%i ',B_indices_long);
	fprintf(fid, '\n');
	
	% Write the coordinates of the stimuli to function_info
	
	fprintf(fid, '%6.2f %6.2f\n', data(:,2:3)');
	fclose(fid);
	
	[out_params,negloglike_hybrid_optimal] = fmincon('nll_hybrid_optimal_vJ',in_params,[],[],[],[],vlb,vub,[],options);
	
	% Compute BIC score
	
	fp=1;
	BIC_hybrid_optimal=2*negloglike_hybrid_optimal+(log(nobs)*fp);
	
	if (BIC_hybrid_optimal < BIC_best)
        
        BIC_best = BIC_hybrid_optimal;
        best_model_BIC = fill_out('Optimal Hybrid',20);
        
	end
	
	if (negloglike_hybrid_optimal < nll_best)
        
        nll_best = negloglike_hybrid_optimal;
        best_model_nll = fill_out('Optimal Hybrid',20);
        
	end
	
	% Report results
	
	fprintf(fid2, 'Switch tag: %i\n', switch_tag);
	fprintf(fid2, 'Noise: %10.5f\n', out_params(1));
	fprintf(fid2, 'Negative Log Likelihood: %10.5f\n', negloglike_hybrid_optimal);
	fprintf(fid2, 'BIC for Hybrid (Optimal Case): %10.5f\n\n', BIC_hybrid_optimal);
    
    
end
    
% Now we write which model fit best to file.

fprintf(fid2, 'Best fitting model by BIC:\n');
fprintf(fid2, '%s\n', best_model_BIC);
fprintf(fid2, 'BIC for that model:\n');
fprintf(fid2, '%10.5f\n', BIC_best);
fprintf(fid2, 'Best fitting model by Negative Log Likelihood:\n');
fprintf(fid2, '%s\n', best_model_nll);
fprintf(fid2, 'Negative Log Likelihood for that model:\n');
fprintf(fid2, '%10.5f\n', nll_best);

fclose(fid2);