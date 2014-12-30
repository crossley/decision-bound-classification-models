
% This program uses a z-grid to fit the hybrid model to a data set in the form
% data = [response x y correct_category].  For the purposes of this
% program, the B category must be in the upper left.  

close all
clear all

% We declare how large we'd like the z-grid (maximum width = 1000), and the
% smallest step size used to compute numerical derivatives.

grid_width = 100;
dmc = 0.005;

% All of the model fitting programs will access the same function_info.dat
% file to retreive the desired size of the z-grid, the number of
% observations, the indices of the A and B responses, and the coordinates of the data.

fid=fopen('function_info.dat','w');
fprintf(fid, '%i\n', grid_width);

% Load the data 

load_b_stim
data = data4;

x = data(:,2); 
y = data(:,3);

% The user visually selects the best x criterial bound, then select two points that determine
% the linear part of the hybrid bound.

[xc yc] = load_bounds(data);
[point1(1) point1(2)] = load_first_point(data);
[point2(1) point2(2)] = load_second_point(data);

cat_org = 0; %input('Switch Categories (1 if yes, 0 if no):');
noise = 5; %input('Starting Noise:');

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
    
num_A = length(A_indices);
num_B = length(B_indices);

fprintf(fid, '%i %i %i\n',nobs, num_A, num_B);

% Write the A and B indices to function_info

A_indices_long = A_indices;
B_indices_long = B_indices;

A_indices_long((num_A+1:1000))=0;
B_indices_long((num_B+1:1000))=0;

fprintf(fid, '%i ',A_indices_long);
fprintf(fid, '\n');
fprintf(fid, '%i ',B_indices_long);
fprintf(fid, '\n');

% Write the coordinates of the stimuli to function_info

fprintf(fid, '%6.2f %6.2f\n', data(:,2:3)');
fclose(fid);

% We creat the linear part of the hybrid bound.

[a1 a2 b] = line_from_points(point1,point2);

lin_bound = [a1 a2 b];

in_params = [a1 b xc noise noise];

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

vlb = [-1 -100 0 0.001 0.001];
vub = [1 100 100 10 10];

% We set the search options.

options = optimset(...
     'Display', 'off',...  % Spit out progress info
 'DiffMinChange', dmc,...    
 'TolFun', .001);      % Termination tolerance on f    

% Run model fitting program

fprintf('\n\nFitting Unequal Hybrid\n'); 

[out_params,negloglike_hybrid_uneq] = fmincon('nll_hybrid_uneq_vJ',in_params,[],[],[],[],vlb,vub,[],options);

% Compute BIC score

fp=5;
BIC_hybrid_uneq=2*negloglike_hybrid_uneq+(log(nobs)*fp);

% Report results

a1 = out_params(1);
a2 = sqrt(1-a1^2);
b = out_params(2);

fprintf('a1: %10.5f\n', a1);
fprintf('a2: %10.5f\n', a2);
fprintf('b: %10.5f\n', b);
fprintf('Unidimensional X Bound: %10.5f\n', out_params(3));
fprintf('Noise X: %10.5f\n', out_params(4));
fprintf('Noise Y: %10.5f\n', out_params(5));
fprintf('Negative Log Likelihood: %10.5f\n', negloglike_hybrid_uneq);
fprintf('BIC for Hybrid (Unequal Variance Case): %10.5f\n\n', BIC_hybrid_uneq);

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
