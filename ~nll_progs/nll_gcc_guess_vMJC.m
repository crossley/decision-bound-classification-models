function negloglike = nll_gcc_guess_vMJC(in_params)

global A_indices B_indices data z_limit

%  negloglike = nll_gcc_eq_vJ(in_params)
%  returns the negative loglikelihood of the 2d data for the
%  General Conjuctive Classifier with unequal variance in the two dimensions.

%  Parameters:
%    params format: [biasX biasY noise] (so x = biasX and
%    y = biasY make boundary)
%    data row format:  [subject_response x y correct_response]
%    z_limit is the z-score value beyond which one should truncate

% The B category is assumed to be in the upper left. The A category is
% assumed to be in the lower right. The model guesses in the remaining 
% two quadrants.

xc = in_params(1);
yc = in_params(2);
noise = in_params(3);

zscoresX1 = (data(:,2)-xc)./noise; 
zscoresY1 = (data(:,3)-yc)./noise;
zscoresX1 = min(zscoresX1,z_limit);
zscoresX1 = max(zscoresX1,-z_limit);
zscoresY1 = min(zscoresY1,z_limit);
zscoresY1 = max(zscoresY1,-z_limit);

zscoresX2 = (xc-data(:,2))./noise; 
zscoresY2 = (yc-data(:,3))./noise;
zscoresX2 = min(zscoresX2,z_limit);
zscoresX2 = max(zscoresX2,-z_limit);
zscoresY2 = min(zscoresY2,z_limit);
zscoresY2 = max(zscoresY2,-z_limit);

pXB1=normalcdf(zscoresX1);
pYB1=normalcdf(zscoresY1);

pXB2=normalcdf(zscoresX2);
pYB2=normalcdf(zscoresY2);

prB = (pXB2).*(pYB1) + 0.5*(pXB2).*(pYB2) + 0.5*(pXB1).*(pYB1);
prA = 1-prB;

log_A_probs = log(prA(A_indices));
log_B_probs = log(prB(B_indices));

% Sum up loglikelihoods and return the negative

negloglike = -(sum(log_A_probs)+sum(log_B_probs));