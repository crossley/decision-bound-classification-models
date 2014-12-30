% Script to load data, reformat it, allow user to input starting guesses
% for the bounds and run model fit.

clear all
% home
cd(fileparts(mfilename('fullpath')))

global data subjects use_resp fit_hybrid num_trials blocks block_len
global fit_guessing fit_biased_guessing fit_unix fit_uniy fit_unequal_GCC
global fit_equal_GCC fit_GLC fit_GQC fit_GCC_guessing fit_GCC_eq2 
global fit_unix_learn fit_uniy_learn
global by_session session_block_len


% User specified variables
subs_delay = [201:220];
subs_immed = [101:121];
subs_maip = [301:317];
subs_rb = [1:30];

subjects = [subs_immed subs_delay subs_maip subs_rb];
use_resp = 1; %use_resp = 1 means use the actual subject responses, 0 means use true categories
num_trials = 0;
block_len = 0;
blocks = num_trials / block_len;

by_session = 1;
session_block_len = 0; % set to zero to use all trials of each session

fit_unix_learn = 1;
fit_uniy_learn = 0;

fit_guessing = 0;
fit_biased_guessing = 0;
fit_unix = 0;
fit_uniy = 0;
fit_unequal_GCC = 0;
fit_equal_GCC = 0;
fit_GCC_guessing = 0;
fit_GCC_eq2 = 0;
fit_GLC = 0;
fit_GQC = 0;

% We run the script that runs the actual model fits
run_model_fits;