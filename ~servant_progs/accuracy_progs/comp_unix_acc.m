% This script computes the accuracy of a given unidimensional x rule when
% used to determine category membership as given by subject###.dat

% We start by loading the subject file

sub_ind = 708;

% We define the bound we want to compute accuracy for

xc = 0;

label = ['subject' num2str(sub_ind)];

load_cmd = ['load  ' cd '\data\' label '.dat;'];
eval(load_cmd);	

def_cmd = ['data = ' label ';'];
eval(def_cmd);

% We compute the accuracy

corr_1 = length(find((data(:,1) == 1) & (data(:,2) < xc)));
corr_2 = length(find((data(:,1) == 2) & (data(:,2) > xc)));

accuracy = (corr_1 + corr_2)/length(data(:,1))
   