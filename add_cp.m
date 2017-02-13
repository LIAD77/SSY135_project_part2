function output = add_cp(input, prefix_length)
%ADD_CP Add cyclic prefix to a signal.
% input - The input signal.
% prefix_length - The length of the prefix .

if prefix_length > length(input)
    error('Cyclic prefix must be shorter than the signal length.')
end

output = [input(end-prefix_length+1:end) ; input];
