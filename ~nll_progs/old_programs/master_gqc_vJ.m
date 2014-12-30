
% Master program for fitting GQC model to 2 dimensional stimuli

close all
clear all

global z_limit data_info

z_limit = 3;

% Load the data 

%load_b_stim
%data = data1;

load_cmd = ['load  ' cd '\data\subject501.dat;'];
eval(load_cmd);	

data = subject501;

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

lin_bound=fisherdiscrim2d(data(:,1:4),5);

% The definition of in_params depends on b1^2 + b2^2 = 1 and on b2 >= 0, so
% we normalize just to be safe.

b1 = lin_bound(1)/sqrt(lin_bound(1)^2+lin_bound(2)^2);
c = lin_bound(3)/sqrt(lin_bound(1)^2+lin_bound(2)^2);

if (lin_bound(2) < 0) 

    b1 = -b1;
    c = -c;
    
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

in_params = [0 0 0 b1 c noise noise];

% In addition, we redefine data for the sake of negloglike_2dGLC.  

data_info = data;
data_info(:,4) = ones(length(data(:,4)),1);

% Set contraints on parameters (Note that b1 varies between -1 and 1)

vlb = [-100 -100 -100 -1 -100 0.001 0.001];
vub = [100 100 100 1 100 10 10];

dmc = 0.005;

% We set the search options.

options = optimset(...
     'Display', 'off',...  % Spit out progress info
 'DiffMinChange', dmc,...    
 'TolFun', .001);      % Termination tolerance on f    

% Run model fitting program

fprintf('\n%s\n\n',['Running GQC Fits...']); 

[out_params,negloglike_gqc] = fmincon('nll_gqc_vJ',in_params,[],[],[],[],vlb,vub,[],options);

% Compute BIC score

fp=6;
BIC_gqc=2*negloglike_gqc+(log(nobs)*fp);

% Plot results

% a1 =    0.00328;
% a2 =   -0.00536;
% a3 =   -0.00094;
% b1 =   -0.62144;
% b2 =    0.78346;
% c =   -1.60052;

a1 = out_params(1);
a2 = out_params(2);
a3 = out_params(3);

b1 = out_params(4);
b2 = sqrt(1-b1^2);

c = out_params(5);

plot_bound = [b1 b2 c];

figure;
plot2dstim(data,xyaxes,0);
hold on;

title(['Subjects scores Calculated Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');

bnd = [a1 a2 a3 b1 b2 c];
range = [-200 200]; 
plotstr = 'r-';

plot2dquadbnd(bnd, range, plotstr)

% Report results

fprintf('a1: %10.5f\n', a1);
fprintf('a2: %10.5f\n', a2);
fprintf('a3: %10.5f\n', a3);
fprintf('b1: %10.5f\n', b1);
fprintf('b2: %10.5f\n', b2);
fprintf('c: %10.5f\n', c);
fprintf('Perceptual Noise: %10.5f\n', out_params(6));
fprintf('Criterial Noise: %10.5f\n', out_params(7));
fprintf('Negative Log Likelihood: %10.5f\n', negloglike_gqc);
fprintf('BIC for GQC: %10.5f\n\n', BIC_gqc);
