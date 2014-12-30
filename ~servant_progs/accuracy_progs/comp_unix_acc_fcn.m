function [accuracy, max_accuracy, xc_max]  = comp_unix_acc_fcn(xc_vector)

% This script computes the accuracy of a given unidimensional x rule when
% used to determine category membership as given by subject###.dat

% We start by loading the subject file

sub_ind = 501;

label = ['subject' num2str(sub_ind)];

load_cmd = ['load  ' cd '\data\' label '.dat;'];
eval(load_cmd);	

def_cmd = ['data = ' label ';'];
eval(def_cmd);

% We compute the accuracy for various values of xc_vector

for i = 1:length(xc_vector)
	
	corr_1 = length(find((data(:,1) == 1) & (data(:,2) < xc_vector(i))));
	corr_2 = length(find((data(:,1) == 2) & (data(:,2) > xc_vector(i))));
	
	accuracy(i) = (corr_1 + corr_2)/length(data(:,1));
       
end

max_accuracy = 0;
xc_max = 0;

           