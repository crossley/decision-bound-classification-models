function negloglike = nll_unix_vJ(in_params)

%  negloglike = nll_unix_vJ(in_params)
%  returns the negative loglikelihood of the unidimensional X
%  bound fit

%  Parameters:
%    params format:  [bias noise] (so x=bias is boundary)
%    data row format:  [subject_response x y correct_response]
%    z_limit is the z-score value beyond which one should truncate

%  We assume the B category is on the left.  

global A_indices B_indices z_limit data

xc = in_params(1);
noise = in_params(2);

zscoresX = (data(:,2)-xc)./noise; 
zscoresX = min(zscoresX,z_limit);
zscoresX = max(zscoresX,-z_limit);

pXA = normcdf(zscoresX, 0.0, 1.0);

prA = pXA;
log_A_probs = log(prA(A_indices));

prB = 1-prA;
log_B_probs = log(prB(B_indices));

% Sum up loglikelihoods and return the negative

negloglike = -(sum(log_A_probs)+sum(log_B_probs));
