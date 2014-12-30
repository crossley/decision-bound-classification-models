
% Master program for fitting GLC model to 2 dimensional stimuli

close all
clear all

global z_limit data_info

z_limit = 3;




% Load the data 

use_resp = 1;
sub_num = 501
block_num = 1;
trials_per_block = 100;

% Load the data for this subject

cd .. % We move to the model_fitting directory so we can access the data.

    all_data = load_sub_data(sub_num,use_resp); % Load the data for this subject

cd GUI % We move back to the GUI directory

rel_trials = (trials_per_block*(1-block_num)+1):(trials_per_block*block_num);

global z_limit A_indices B_indices data

z_limit = 3;
noise = 5; 
% 
% load_b_stim
% data = data4;

 data = all_data(rel_trials,:);

x = data(:,2); 
y = data(:,3);

% We give the user the chance to switch the category labels.

cat_org = 1; % input('Switch Categories (1 if yes, 0 if no):');

[nobs col] = size(data);

global A_indices B_indices

% We switch the responses if necessary

if cat_org == 0 
    
    A_indices = find(data(:,1) == 1);
    B_indices = find(data(:,1) ~= 1);

else

    A_indices = find(data(:,1) ~= 1);
    B_indices = find(data(:,1) == 1);

end
    
% We'll start with the best bound from the actual categories using fisherdiscrim2d.
% At this point bound is in the form [a1 a2 b], where a1*x + a2*y + b = 0
% is the boundary.

lin_bound=fisherdiscrim2d([data(:,1:3) data(:,1)],5);

% The definition of in_params depends on a1^2 + a2^2 = 1 and on a2 >= 0, so
% we normalize just to be safe.

a1 = lin_bound(1)/sqrt(lin_bound(1)^2+lin_bound(2)^2);
b = lin_bound(3)/sqrt(lin_bound(1)^2+lin_bound(2)^2);

if (lin_bound(2) < 0) 

    a1 = -a1;
    b = -b;
    
end

% We plot this starting bound

figure;
xyaxes=[min(x) max(x) min(y) max(y)];
plot2dstim(data,xyaxes,0);
hold on;

title(['Subjects scores Start Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');
plot2dlinbnd(lin_bound,'k-'); 

% We define our initial parameters

noise=5;

in_params = [a1 b noise];

% In addition, we redefine data for the sake of negloglike_2dGLC.  

data_info = data;
data_info(:,4) = ones(length(data(:,4)),1);

% Set contraints on parameters (Note that a1 varies between -1 and 1)

vlb = [-1 -100 0.001];
vub = [1 100 10];

% Set search options
options = optimset(...
     'Display', 'iter',...  % Spit out progress info
    'TolFun', .001');      % Termination tolerance on f    

% Run model fitting program

fprintf('\n%s\n\n',['Running Equal Variance GLC Fits...']); 

[out_params,negloglike_glc_eq] = fmincon('nll_glc_eq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

% Compute BIC score

fp=3;
BIC_glc_eq=2*negloglike_glc_eq+(log(nobs)*fp);

% Plot results

a1 = out_params(1);
a2 = sqrt(1-a1^2);
b = out_params(2);

plot_bound = [a1 a2 b];

figure;
plot2dstim(data,xyaxes,0);
hold on;

title(['Subjects scores Calculated Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');

plot2dlinbnd(plot_bound,'r-'); 

% Report results

fprintf('a1: %10.5f\n', a1);
fprintf('a2: %10.5f\n', a2);
fprintf('b: %10.5f\n', b);
fprintf('Noise: %10.5f\n', out_params(3));
fprintf('Negative Log Likelihood: %10.5f\n', negloglike_glc_eq);
fprintf('BIC for GLC (Equal Variance Case): %10.5f\n\n', BIC_glc_eq);
