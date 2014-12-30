% This script runs the actual model fits.  This script expects the data to
% be in the cd\data directory, in the form subject###.dat.  For each
% subject, it also expects a file called start_vals_sub_###.dat that
% contains one line for each block, and contains the starting guesses for
% the unidimensional rules, the gcc, and the hybrid center and linear
% points (See record_starting_bounds.m for more info).

global subjects use_resp num_trials blocks block_len 
global by_session session_block_len

ignore_trials = num_trials-block_len*blocks;

% We run the model fits and output all the of best fits to file.

fid = fopen([cd '/~output/all_results.dat'],'w');

fprintf(fid, '%20s %20s %20s %20s %20s %20s %20s %20s\n', fill_out('Subject Number',20),fill_out('Block Number',20),fill_out('Trials per Block',20),fill_out('Accuracy',20),fill_out('Best Model by BIC',20),fill_out('BIC Score',20),fill_out('Best Model by NLL',20),fill_out('NLL Score',20));

for sub_ind = subjects
    
    % Open a new file to hold onto to raw BIC scores 
    % (with no descriptive text added) - MJC 9/17/07
    label=['/~output/raw_BIC_sub' num2str(sub_ind) '.dat'];
    fid_BIC = fopen([cd label], 'w');
    
    % We will also open a file to store raw accuracy
    label=['/~output/raw_accuracy_sub' num2str(sub_ind) '.dat'];
    fid_accuracy = fopen([cd label], 'w');
    
    % We will also open several files to store best fitting bound params
    label=['/~output/unix_learn_params_sub' num2str(sub_ind) '.dat'];
    fid_unix_learn_params = fopen([cd label], 'w');
    
    label=['/~output/uniy_learn_params_sub' num2str(sub_ind) '.dat'];
    fid_uniy_learn_params = fopen([cd label], 'w');
    
    label=['/~output/unix_params_sub' num2str(sub_ind) '.dat'];
    fid_unix_params = fopen([cd label], 'w');
    
    label=['/~output/uniy_params_sub' num2str(sub_ind) '.dat'];
    fid_uniy_params = fopen([cd label], 'w');
    
    label=['/~output/GLC_params_sub' num2str(sub_ind) '.dat'];
    fid_GLC_params = fopen([cd label], 'w');
    
    % Open a file to store percent responses accounted for: various bounds
    label=['/~output/percent_responses_accounted' num2str(sub_ind) '.dat'];
    fid_responses_accounted = fopen([cd label], 'w');

    % Load the data for this subject
    all_data = load_sub_data(sub_ind,use_resp);
    if (length(all_data) < num_trials) && by_session ~= 1
       num_missing = num_trials - length(all_data);
       all_data = [all_data(1:num_missing,:); all_data];
    end
    
    if by_session == 1
        
        sessions = max(all_data(:,5));
        
        for session_num = 1:sessions
            
            rel_trials = find(all_data(:,5)==session_num);
            data = all_data(rel_trials,:);
            
            if session_block_len ~= 0
                num_skip = length(data) - session_block_len;
                data = data((num_skip+1):length(data),:);
            end

            % Name the file for the current subject on the current block
            file_name = ['Sub' num2str(sub_ind) 'Block' num2str(session_num)];

            % [uni_bound gcc_bound hybrid_center_point hybrid_lin_point];
            start_vals = [50 50 50 50 50 50 0 0];

            % Run the model fits

            fprintf('\nSubject: %4i\n', sub_ind);
            fprintf('Block:   %4i\n', session_num);
            [best_model_BIC,BIC_best,best_model_nll,nll_best] = master_function(start_vals,file_name,fid_BIC,fid_unix_params,fid_uniy_params,fid_GLC_params,fid_responses_accounted,fid_unix_learn_params,fid_uniy_learn_params);

            accuracy = compute_accuracy(data);

            % We create a fixed width output file
            fprintf(fid,'%-20i %-20i %-20i %-20.10f %20s %-20.10f %20s %-20.10f\n',sub_ind,session_num,block_len,accuracy,best_model_BIC,BIC_best,best_model_nll,nll_best);

            % Write the raw accuracy information to a file
            fprintf(fid_accuracy,'%10.5f',accuracy);
        
        end
        
    else
        
        for block_num = 1:blocks
        
            rel_trials = ((block_num-1)*block_len + 1):block_num*block_len;
            rel_trials = rel_trials + ignore_trials;
            data = all_data(rel_trials,:);

            % Name the file for the current subject on the current block
            file_name = ['Sub' num2str(sub_ind) 'Block' num2str(block_num)];

            % [uni_bound gcc_bound hybrid_center_point hybrid_lin_point];
            start_vals = [50 50 50 50 50 50 0 0];

            % Run the model fits

            fprintf('\nSubject: %4i\n', sub_ind);
            fprintf('Block:   %4i\n', block_num);
            [best_model_BIC,BIC_best,best_model_nll,nll_best] = master_function(start_vals,file_name,fid_BIC,fid_unix_params,fid_uniy_params,fid_GLC_params,fid_responses_accounted);

            accuracy = compute_accuracy(data);

            % We create a fixed width output file
            fprintf(fid,'%-20i %-20i %-20i %-20.10f %20s %-20.10f %20s %-20.10f\n',sub_ind,block_num,block_len,accuracy,best_model_BIC,BIC_best,best_model_nll,nll_best);

            % Write the raw accuracy information to a file
            fprintf(fid_accuracy,'%10.5f',accuracy);
        
        end
        
    end
    
    % close BIC file
    fclose(fid_BIC);
    fclose(fid_accuracy);
    fclose(fid_unix_learn_params);
    fclose(fid_uniy_learn_params);
    fclose(fid_unix_params);
    fclose(fid_uniy_params);
    fclose(fid_GLC_params);
    fclose(fid_responses_accounted);
    
end
   
fclose(fid);