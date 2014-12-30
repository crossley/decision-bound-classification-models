function negloglike = nll_gcc_uneq_vJ(in_params)

global A_indices B_indices data z_limit

%  negloglike = nll_gcc_uneq_vJ(in_params)
%  returns the negative loglikelihood of the 2d data for the
%  General Conjuctive Classifier with unequal variance in the two dimensions.

%  Parameters:
%    params format: [biasX biasY noiseX noiseY] (so x = biasX and
%    y = biasY make boundary)
%    data row format:  [subject_response x y correct_response]
%    z_limit is the z-score value beyond which one should truncate

% The B category is assumed to be in the upper left.

xc = in_params(1);
yc = in_params(2);
noiseX = in_params(3);
noiseY = in_params(4);

zscoresX = (xc-data(:,2))./noiseX; 
zscoresY = (data(:,3)-yc)./noiseY;
zscoresX = min(zscoresX,z_limit);
zscoresX = max(zscoresX,-z_limit);
zscoresY = min(zscoresY,z_limit);
zscoresY = max(zscoresY,-z_limit);

pXB=normalcdf(zscoresX);
pYB=normalcdf(zscoresY);

prB = (pXB).*pYB;
prA = 1-prB;

log_A_probs = log(prA(A_indices));
log_B_probs = log(prB(B_indices));

% Sum up loglikelihoods and return the negative

negloglike = -(sum(log_A_probs)+sum(log_B_probs));
