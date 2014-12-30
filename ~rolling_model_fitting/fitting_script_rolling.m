% Script to load data, reformat it, allow user to input starting guesses
% for the bounds and run model fit.

clear all
home

global data subjects use_resp fit_hybrid num_trials window_size step_size
 
% User specified variables

subjects = [100];

% subjects = [1010 2040 3010 3020 3040 3050 4010]; [10101 10102 20101 20102 20201 20202 20301 20302 20402 20402 30101 30102 30201 30202 30301 30302 30401 30402 30501 30502 40101 40102 40201 40202 40301 40302 40401 40402];
use_resp = 1; % use_resp = 1 means use the actual subject responses, 0 means use true categories
fit_hybrid = 2; % 2 means fit all hybrids, 1 means fit all hybrids but unequal, 0 means fit no hybrid models
num_trials = 600;
window_size = 50;
step_size = 25;

num_windows = (num_trials-window_size)/step_size + 1

for window_num = 1:num_windows

    window_num
    
((window_num-1)*step_size + 1):(window_size+(window_num-1)*step_size)

end

% We transform the output of hybrid_experiment.m and write the results to
% subject###.dat files in the data directory

% transform_data;

% We run the script that records the starting bounds

% record_starting_bounds;

% We run the script that runs the actual model fits

run_model_fits_rolling;

