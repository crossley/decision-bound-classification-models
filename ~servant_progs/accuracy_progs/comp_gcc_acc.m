% This script computes the accuracy of a given unidimensional x rule when
% used to determine category membership as given by subject###.dat

% We start by loading the subject file

load_cmd = ['load  ' cd '\data\subject501.dat;'];
eval(load_cmd);	

data = subject501;

% We define the bound we want to use to compute accuracy

xc = 33.7;
yc = 20;

% We compute the accuracy

corr_1 = length(find( (data(:,1) == 1) & ((data(:,2) < xc) & (data(:,3) > xc)) ));
corr_2 = length(find( (data(:,1) == 2) & ~((data(:,2) < xc) & (data(:,3) > xc)) ));

accuracy = (corr_1 + corr_2)/length(data(:,1))
   