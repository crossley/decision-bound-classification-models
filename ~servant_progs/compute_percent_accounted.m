function [P] = compute_percent_accounted(block_data, bound_params)

% This function assumes bound params is in the format:
% [a1 a2 b]

a1 = bound_params(1);
a2 = bound_params(2);
b = bound_params(3);

if a1 ~= 0 && a2 ~= 0 % GLC

    m = -a1/a2;
    y_int = -b/a2;
    
    respond_A_ind = find(block_data(:,1)==1);
    x = block_data(respond_A_ind,2);
    y = block_data(respond_A_ind,3);
    num_consistent_A = length(find(y > (m*x + y_int)));

    respond_B_ind = find(block_data(:,1)==2);
    x = block_data(respond_B_ind,2);
    y = block_data(respond_B_ind,3);
    num_consistent_B = length(find(y <= (m*x + y_int)));

elseif a2 == 0 % unix
    
    respond_A_ind = find(block_data(:,1)==1);
    x = block_data(respond_A_ind,2);
    y = block_data(respond_A_ind,3);
    num_consistent_A = length(find(x < -b));

    respond_B_ind = find(block_data(:,1)==2);
    x = block_data(respond_B_ind,2);
    y = block_data(respond_B_ind,3);
    num_consistent_B = length(find(x >= -b));
    
elseif a1 == 0 % uniy
    
    respond_A_ind = find(block_data(:,1)==1);
    x = block_data(respond_A_ind,2);
    y = block_data(respond_A_ind,3);
    num_consistent_A = length(find(y < -b));

    respond_B_ind = find(block_data(:,1)==2);
    x = block_data(respond_B_ind,2);
    y = block_data(respond_B_ind,3);
    num_consistent_B = length(find(y >= -b));
    
end

total_consistent = num_consistent_A + num_consistent_B;
P = total_consistent / size(block_data,1);