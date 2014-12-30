function negloglike = nll_hybrid_eq_si(in_params)

% negloglike = nll_hybrid_eq_vJ(in_params) computes the negative log
% likelihood for a pattern of responses given a hybrid decision rule in the
% form in_params = [a1 b x_criterion noise1 noise2] where a1*x + a2*y + b =
% 0 is the equation of the linear bound, with a1^2 + a2^2 = 1 and a2 >=0.

global A_indices B_indices z_limit data

a1 = in_params(1);
a2 = sqrt(1-a1^2);
b = in_params(2);
B = in_params(3); % B is the x rule criterion
sigma_x = in_params(4);
sigma_y = in_params(4); % We use the same value for the noise in each dimension

C = -a1/a2;
D = -b/a2;

% For each data point we must compute the probability of that data point
% coming from the bound specified by in_params.

for i = 1:length(data(:,1))
     prA(i) = compute_hybrid_prob(data(i,2),data(i,3),sigma_x,sigma_y,B,C,D);
end

prB = 1-prA;

log_A_probs = log(prA(A_indices));
log_B_probs = log(prB(B_indices));

% Sum up loglikelihoods and return the negative

negloglike = -(sum(log_A_probs)+sum(log_B_probs));

