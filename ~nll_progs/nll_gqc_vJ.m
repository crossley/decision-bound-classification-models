function negloglike = nll_gqc_vJ(in_params)

global A_indices B_indices data_info z_limit

%  negloglike = nll_gqc_vJ(params)
%  returns the negative loglikelihood of the 2d data for the
%  General Quadratic Classifier.

%  Note:  These computations are based on the likelihood model

%  Parameters:
%    params format:  [Amat_entries(1:3) b1 c_bias pnoise cnoise]
%    data_info row format:  [correct_response x y 1]
%    z_limit is the z-score value beyond which one should truncate

pnoise = in_params(6);
cnoise = in_params(7);

% Combine individual noises into single noise parameter

noise = pnoise+cnoise;

% Define the A matrix, the B vector, and c constant

A = [   in_params(1) .5*in_params(3)
     .5*in_params(3)    in_params(2)];
b = [in_params(4) sqrt(1-in_params(4)^2)]';
c = in_params(5);

% Compute mean of h(x) for each x (as shown in Eqn. 11, Ashby p.462).

temp1 = trace(A*noise);
xAxs = in_params(1)*data_info(:,2).^2 + in_params(2)*data_info(:,3).^2 + (in_params(3)*prod(data_info(:,2:3)'))';
meanhxs = temp1 + xAxs + data_info(:,2:3)*b + c;

% Compute variance of h(x) for each x (as shown in Eqn. 12, Ashby p.462).
% Take advantage of fact that perceptual noise matrix is diagonal with
% pnoise the same in all three dimensions

ndatapts = size(data_info,1);
bvec = [b(1)*ones(1,ndatapts);b(2)*ones(1,ndatapts)];
temp2 = bvec+2*A*data_info(:,2:3)';
varhxs = 2*temp1*temp1 + noise*(sum(temp2.^2))';

% Compute z-scores for each data point

zscores = meanhxs./(sqrt(varhxs));   
zscores = min(zscores,z_limit);
zscores = max(zscores,-z_limit);
   
% Find log of cumulative probability

prB = normcdf(zscores);

prA = 1-prB;

log_A_probs = log(prA(A_indices));
log_B_probs = log(prB(B_indices));
   
% Sum them up and return the negative

negloglike = -(sum(log_A_probs)+sum(log_B_probs));




   
   

