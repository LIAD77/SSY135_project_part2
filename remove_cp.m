function output = remove_cp(input, prefix_length)
%REMOVE_CP Remove cyclic prefix of a signal.
% input - The input signal.
% prefix_length - The length of the prefix .
output = input(prefix_length+1:end);
