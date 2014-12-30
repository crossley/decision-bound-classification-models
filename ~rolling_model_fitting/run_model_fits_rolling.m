% This script runs the actual model fits.  This script expects the data to
% be in the cd\data directory, in the form subject###.dat.  For each
% subject, it also expects a file called start_vals_sub_###.dat that
% contains one line for each block, and contains the starting guesses for
% the unidimensional rules, the gcc, and the hybrid center and linear
% points (See record_starting_bounds.m for more info).

global data subjects use_resp fit_hybrid num_trials window_size step_size

% We run the model fits and output all the of best fits to file.

fid = fopen([cd '\output\all_results.dat'],'w');

fprintf(fid, '%20s %20s %20s %20s %20s %20s %20s %20s\n', fill_out('Subject Number',20),fill_out('Starting Trial',20),fill_out('Ending Trial',20),fill_out('Accuracy',20),fill_out('Best Model by BIC',20),fill_out('BIC Score',20),fill_out('Best Model by NLL',20),fill_out('NLL Score',20));

num_windows = (num_trials-window_size)/step_size + 1;

for sub_ind = subjects

    % Load the data for this subject
    
    all_data = load_sub_data(sub_ind,use_resp);

    % We load the starting bounds for this subject
    
    start_vals = [50 50 50 50 50 50 0 0]; % [uni_bound gcc_bound hybrid_center_point hybrid_lin_point];

    for window_num = 1:num_windows

        start_trial = ((window_num-1)*step_size + 1);
        end_trial = (window_size+(window_num-1)*step_size);
        data = all_data(start_trial:end_trial,:);
	
		% Name the file for the current subject on the current block
		
		file_name = ['Sub' num2str(sub_ind) 'Window' num2str(window_num)];
        
%         uni_bound =  [30 30]; %[start_bounds(block_num,1), start_bounds(block_num,2)];
%         gcc_bound = [40 30]; %[start_bounds(block_num,3), start_bounds(block_num,4)];
%         hybrid_center_point = [50 50]; %[start_bounds(block_num,5), start_bounds(block_num,6)];
%         hybrid_lin_point = [0 0]; %[start_bounds(block_num,7), start_bounds(block_num,8)];
    
		% Run the model fits
		
        fprintf('\nSubject: %4i\n', sub_ind);
        fprintf('Window:   %4i\n', window_num);
		[best_model_BIC,BIC_best,best_model_nll,nll_best,out_vals] = master_function_rolling(start_vals,file_name,fit_hybrid);

        accuracy = compute_accuracy(data);
        
        accuracy = max(accuracy, 1-accuracy);
        
% 		We create a fixed width output file.

        fprintf(fid,'%-20i %-20i %-20i %-20.10f %20s %-20.10f %20s %-20.10f\n',sub_ind,start_trial,end_trial,accuracy,best_model_BIC,BIC_best,best_model_nll,nll_best);

%       We use the bounds from the previous trial as the starting bounds
%       for the next

        start_vals = out_vals;
        
        
    end
    
end
   
fclose(fid);
