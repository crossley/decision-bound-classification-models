function [P] = compute_percent_accounted_old(block_data, bound_params)

if size(bound_params,2) < 4

    a1 = bound_params(1);
    a2 = bound_params(2);
    b = bound_params(3);

    if -a1/a2 >= 0 % positve slope bound

        respond_A_ind = find(block_data(:,1)==1);
        x = block_data(respond_A_ind,2);
        y = block_data(respond_A_ind,3);
        num_consistent_A = length(find(a1*x+a2*y+b >= 0));

        respond_B_ind = find(block_data(:,1)==2);
        x = block_data(respond_B_ind,2);
        y = block_data(respond_B_ind,3);
        num_consistent_B = length(find(a1*x+a2*y+b < 0));

    else

        respond_A_ind = find(block_data(:,1)==1);
        x = block_data(respond_A_ind,2);
        y = block_data(respond_A_ind,3);
        num_consistent_A = length(find(a1*x+a2*y+b < 0));

        respond_B_ind = find(block_data(:,1)==2);
        x = block_data(respond_B_ind,2);
        y = block_data(respond_B_ind,3);
        num_consistent_B = length(find(a1*x+a2*y+b >= 0));

    end

elseif size(bound_params,2) == 4
    
    a1 = bound_params(1);
    a2 = bound_params(2);
    b = bound_params(3);
    x_bound = bound_params(4);
    
    num_consistent_A = 0;
    num_consistent_B = 0;
    
    % Anything greater than x_bound is a B
    num_consistent_B = num_consistent_B + length(find(block_data(:,1)==2 & block_data(:,2) > x_bound));
    
    % Less than x_bound and below the line a1*x + a2*y + b = 0 is also a B
    respond_B_ind = find(block_data(:,1)==2 & block_data(:,2)<=x_bound);
    x = block_data(respond_B_ind,2);
    y = block_data(respond_B_ind,3);
    num_consistent_B = num_consistent_B + length(find(a1*x+a2*y+b < 0));

    % Less than x_bound and above the line a1*x + a2*y + b = 0 is an A
    respond_A_ind = find(block_data(:,1)==1 & block_data(:,2)<=x_bound);
    x = block_data(respond_A_ind,2);
    y = block_data(respond_A_ind,3);
    num_consistent_A = num_consistent_A + length(find(a1*x+a2*y+b >= 0));

end

total_consistent = num_consistent_A + num_consistent_B;

P = total_consistent / size(block_data,1);