
% Framework for fitting unidimensional Y rule to 2 dimensional stimuli

close all
clear all

global z_limit data

z_limit = 3;

% Load the data 

% sub_num = 501;
% use_resp = 1;
% 
% trials_per_block = 100;
% block_num = 1;
% 
% cd .. % We move to the model_fitting directory so we can access the data.
% 
%     all_data = load_sub_data(sub_num,use_resp); % Load the data for this subject
% 
% cd GUI % We move back to the GUI directory
% 
% rel_trials = (trials_per_block*(1-block_num)+1):(trials_per_block*block_num);
% 
% data = all_data(rel_trials,:);

load_b_stim
data = data2;

x = data(:,2); 
y = data(:,3);

% We plot the data so we can make our original guesses for the X and Y
% bounds.

xyaxes=[min(x) max(x) min(y) max(y)];

figure;
plot2dstim(data,xyaxes,0);
hold on;

title(['Subject Responses'])
xlabel('Orientation');
ylabel('Spatial Frequency');
legend('A','B');

% We read in our initial guesses for xc, yc and the noise standard
% deviation.  We also give the user the chance to switch the category
% labels.

% cat_org = 0;
% xc = 80;
% yc = 55;
% noise = 5;

cat_org = input('Switch Categories (1 if yes, 0 if no):');
xc = input('Starting X Bound:');
yc = input('Starting Y Bound:');
noise = input('Starting Noise:');

close all

% Write the number of data points to function_info

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

% Plot subject's responses and the starting bound

figure; 
xyaxes=[min(x) max(x) min(y) max(y)];
plot2dstim(data,xyaxes,0);
hold on;

title(['Subject scores Start Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');

% Y bound (plots y = yc)

y_bnd = [0 1 -yc];
plot2dlinbnd(y_bnd,'k-'); 

legend('A','B','y')

% We create our starting guess for the parameters.

in_params = [yc noise];

% Set contraints on parameters

vlb = [0 0.001];
vub = [100 10];

% Set search options
options = optimset(...
     'Display', 'off',...  % Spit out progress info
    'TolFun', .001');      % Termination tolerance on f    

% Run model fitting program

fprintf('\n%s\n\n',['Running Unidimensional Y Fits...']); 

[out_params,negloglike_uniy] = fmincon('nll_uniy_vJ',in_params,[],[],[],[],vlb,vub,[],options);

% Compute BIC score

fp=2;
BIC_uniy=2*negloglike_uniy+(log(nobs)*fp);

% Report results

fprintf('Y Bound: %10.5f\n', out_params(1));
fprintf('Noise: %10.5f\n', out_params(2));
fprintf('Negative Log Likelihood: %10.5f\n', negloglike_uniy);
fprintf('BIC for Unidimensional Y: %10.5f\n\n', BIC_uniy);

% Plot results

figure;
plot2dstim(data,xyaxes,0);
hold on;

title(['Subjects scores Calculated Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');

y_bnd = [0 1 -out_params(1)];
plot2dlinbnd(y_bnd,'r-'); 

legend('A','B','y')

