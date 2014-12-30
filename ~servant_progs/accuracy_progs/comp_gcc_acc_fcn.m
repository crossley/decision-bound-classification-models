function [accuracy, max_accuracy, xc_max, yc_max] = comp_gcc_acc_fcn(xc_vector, yc_vector)

% This script computes the accuracy of a given gcc rule when
% used to determine category membership as given by subject###.dat

% We start by loading the file

sub_ind = 501;

label = ['subject' num2str(sub_ind)];

load_cmd = ['load  ' cd '\data\' label '.dat;'];
eval(load_cmd);	

def_cmd = ['data = ' label ';'];
eval(def_cmd);

clear accuracy 

for i = 1:length(xc_vector)
    for j = 1:length(yc_vector)
    		
		% We compute the accuracy
		
		corr_1 = length(find( (data(:,1) == 1) & ((data(:,2) < xc_vector(i)) & (data(:,3) > yc_vector(j))) ));
		corr_2 = length(find( (data(:,1) == 2) & ~((data(:,2) < xc_vector(i)) & (data(:,3) > yc_vector(j))) ));
		
		accuracy(i,j) = (corr_1 + corr_2)/length(data(:,1));
        
%         if (accuracy(i,j) > max(max(accuracy)))
%             max_accuracy = accuracy(i,j)
%             xc_max = xc_vector(i)
%             yc_max = yc_vector(j)
%         end

    end
end

max_accuracy = 0;
xc_max = 0;
yc_max = 0;
           