function negloglike = nll_gcc_eq2_vMJC(in_params)

global A_indices B_indices data z_limit

%  negloglike = nll_gcc_eq_vJ(in_params)
%  returns the negative loglikelihood of the 2d data for the
%  General Conjuctive Classifier with unequal variance in the two dimensions.

%  Parameters:
%    params format: [biasX biasY noise] (so x = biasX and
%    y = biasY make boundary)
%    data row format:  [subject_response x y correct_response]
%    z_limit is the z-score value beyond which one should truncate

% The A category is assumed to be in the lower right.

xc = in_params(1);
yc = in_params(2);
noise = in_params(3);

zscoresX = (data(:,2)-xc)./noise; 
zscoresY = (yc-data(:,3))./noise;
zscoresX = min(zscoresX,z_limit);
zscoresX = max(zscoresX,-z_limit);
zscoresY = min(zscoresY,z_limit);
zscoresY = max(zscoresY,-z_limit);

pXA=normalcdf(zscoresX);
pYA=normalcdf(zscoresY);

prA = (pXA).*pYA;
prB = 1-prA;

log_A_probs = log(prA(A_indices));
log_B_probs = log(prB(B_indices));

% Sum up loglikelihoods and return the negative

negloglike = -(sum(log_A_probs)+sum(log_B_probs));