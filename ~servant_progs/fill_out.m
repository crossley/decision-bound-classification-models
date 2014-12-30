function new_string = fill_out(old_string,n)

% new_string = fill_out(old_string,n) fills out or truncates a string using blanks so it
% has length n

new_string = [old_string blanks(n - length(old_string))];
new_string = new_string(1:n);