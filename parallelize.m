function parallel = parallelize(serial, channels)
%PARALLELIZE Parallelize a vector by reshaping it into a vector
%with channels rows. The matrix is filled column-wise.
% serial - Input vector.
% channels - The number of channels (rows of the matrix).

if channels > numel(serial)
    error('There must be fewer channels than the length of the input.')
end

if numel(serial) ~= size(serial, 1)
    error('Input must be a column vector.')
end

% Append zeroes if the numel of the input is not a multiple of the
% number of channels.
if mod(numel(serial) / channels, 1) ~= 0
    serial = [serial ; zeros(channels - mod(numel(serial), channels), 1)];
end
parallel = reshape(serial, channels, numel(serial) / channels);

end