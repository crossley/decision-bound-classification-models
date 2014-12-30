function out_data = load_sub_data(sub_num,use_resp)

% out_data = load_sub_data(sub_num,use_resp) loads in data for subject sub_num
 
% This script expects a file entitled "subject###.dat" in the form of
% [correct_category x y response].  Load_sub_data then returns a matrix
% with 4 columns, the first of which is either the actual category label or the
% response and the last of which is always the correct categories

label=['subject' num2str(sub_num)];

load_cmd = ['load  ' cd '/~data/' label '.txt;'];
eval(load_cmd);	

define_cmd = ['temp_data = ' label ';'];
eval(define_cmd);

corr_cat = temp_data(:,1);
x = temp_data(:,2);
y = temp_data(:,3);
resp = temp_data(:,4);
session = [];

if size(temp_data,2) == 5
    session = temp_data(:,5);
end

if (use_resp == 1)

    % To fit responses put responses in column 1 of all_data
    	
	out_data = [resp x y corr_cat session];

else
       
    % To fit the actual categories put category labels in column 1
    
	out_data = [corr_cat x y resp session];

end
