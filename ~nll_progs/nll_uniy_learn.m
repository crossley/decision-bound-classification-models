function negloglike = nll_uniy_learn(in_params)

%  negloglike = nll_uniy_learn(in_params)
%  returns the negative loglikelihood of the unidimensional X
%  - criterial learning variant - bound fit

%  Parameters:
%    params format:  [xc_init noise alpha_pos alpha_neg]
%    data row format:  [subject_response x y correct_response]
%    z_limit is the z-score value beyond which one should truncate

%  We assume that B category is on the left.  

global A_indices B_indices z_limit data

A_ind = A_indices;
B_ind = B_indices;
num_trials = size(data,1);

xc_init = in_params(1);
noise = in_params(2);
alpha_pos = in_params(3);
alpha_neg = in_params(4);

x = data(:,3);
fb = data(:,1) == data(:,4);

xc = zeros(num_trials,1); xc(1) = xc_init;
delta_pos = zeros(num_trials,1);
delta_neg = zeros(num_trials,1);

% init xc (we want to transform into [0 100] or something)
xc(1) = (max(x) - min(x))*rand();
for t = 2:num_trials-1;

	% update xc
	if fb(t) == 1; xc(t+1) = xc(t) + alpha_pos*randn; end % NOTE: the drift is from the previous xc
	if fb(t) == 0; xc(t+1) = x(t) + noise + alpha_neg*randn; end % NOTE: the drift is from the percept

end

zscoresX = (x-xc)./noise; 
zscoresX = min(zscoresX,z_limit);
zscoresX = max(zscoresX,-z_limit);

pXA = normcdf(zscoresX, 0.0, 1.0);

prA = pXA;
log_A_probs = log(prA(A_ind));

prB = 1-prA;
log_B_probs = log(prB(B_ind));

negloglike = -(sum(log_A_probs)+sum(log_B_probs));
