function accuracy = compute_accuracy(data)

corr_count = length(find(data(:,1) == data(:,4)));
accuracy = corr_count/length(data(:,1));