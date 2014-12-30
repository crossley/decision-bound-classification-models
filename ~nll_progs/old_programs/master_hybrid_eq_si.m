
% This program uses a z-grid to fit the hybrid model to a data set in the form
% data = [response x y correct_category].  For the purposes of this
% program, the B category must be in the upper left.  

close all
clear all

% Load the data 

global data

load_b_stim
data = data4;

x = data(:,2); 
y = data(:,3);

% The user visually selects the best x criterial bound, then select two points that determine
% the linear part of the hybrid bound.

global pos % This is the position of the graphical interface

pos = [3 33 1279 921];

[cen_pt(1) cen_pt(2)] = load_hybrid_center_point(data,1,1);
[lin_pt(1) lin_pt(2)] = load_hybrid_lin_point(cen_pt,data,1,1);

xc = cen_pt(1);

cat_org = 0; %input('Switch Categories (1 if yes, 0 if no):');
noise = 5; %input('Starting Noise:');

close all

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

% We create the linear part of the hybrid bound.

[a1 a2 b] = line_from_points(cen_pt,lin_pt);

lin_bound = [a1 a2 b];

in_params = [a1 b xc noise];

% We plot this starting bound

figure;
xyaxes=[min(x) max(x) min(y) max(y)];
plot2dstim(data,xyaxes,0);
hold on;

title(['Subjects scores Start Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');

plot2dlinbnd(lin_bound,'k-'); 
x_bound = [1 0 -xc];
plot2dlinbnd(x_bound,'k-'); 

% Set contraints on parameters (Note that a1 is constrained between -1 and 1)

vlb = [-1 -100 0 0.001];
vub = [1 100 100 10];

% We set the search options.  Note that dmc is the
% smallest step size used to compute numerical derivatives.

dmc = .005;

options = optimset(...
     'Display', 'on',...  % Spit out progress info
 'DiffMinChange', dmc,...    
 'TolFun', .001);      % Termination tolerance on f    

% Run model fitting program

fprintf('\n\nFitting Equal Hybrid\n'); 

warning off MATLAB:quadl:ImproperFcnValue

[out_params,negloglike_hybrid_eq] = fmincon('nll_hybrid_eq_si',in_params,[],[],[],[],vlb,vub,[],options);

% Compute BIC score

fp=4;
BIC_hybrid_eq=2*negloglike_hybrid_eq+(log(nobs)*fp);

% Report results

a1 = out_params(1);
a2 = sqrt(1-a1^2);
b = out_params(2);

fprintf('a1: %10.5f\n', a1);
fprintf('a2: %10.5f\n', a2);
fprintf('b: %10.5f\n', b);
fprintf('Unidimensional X Bound: %10.5f\n', out_params(3));
fprintf('Noise X: %10.5f\n', out_params(4));
fprintf('Negative Log Likelihood: %10.5f\n', negloglike_hybrid_eq);
fprintf('BIC for Hybrid (Equal Variance Case): %10.5f\n\n', BIC_hybrid_eq);

% Plot final bound

figure;
xyaxes=[min(x) max(x) min(y) max(y)];
plot2dstim(data,xyaxes,0);
hold on;

plot_bound = [a1 a2 b];

title(['Subjects scores Calculated Bound'])
xlabel('Orientation');
ylabel('Spatial Frequency');
plot2dlinbnd(plot_bound,'r-'); 
x_bound = [1 0 -out_params(3)];
plot2dlinbnd(x_bound,'r-'); 
